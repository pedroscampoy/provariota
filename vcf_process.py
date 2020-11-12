#!/usr/bin/env python

import os
import logging
import sys
import argparse
import pandas as pd
import numpy as np
import re
import gzip
import subprocess
from misc import check_file_exists, check_create_dir, execute_subprocess, check_remove_file, \
list_to_bed, count_lines


logger = logging.getLogger()

def filter_tsv_variants(tsv_file, output_filtered, min_frequency=0.8, min_total_depth=20, min_alt_dp=4, is_pass=True, only_snp=True):
    input_file_name = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    output_file = os.path.join(output_filtered, input_file_name)

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")
    filtered_df = df[(df.PASS == is_pass) &
                    (df.TOTAL_DP >= min_total_depth) &
                    (df.ALT_DP >= min_alt_dp) &
                    (df.ALT_FREQ >= min_frequency)]
    if only_snp == True:
        final_df = filtered_df[~(filtered_df.ALT.str.startswith('+') | filtered_df.ALT.str.startswith('-'))]
        final_df.to_csv(output_file, sep='\t', index=False)
    else:
        filtered_df.to_csv(output_file, sep='\t', index=False)
    

def define_var_type(row):
    len_ref = len(row.REF)
    len_alt = len(row.ALT)
    
    if len_ref == len_alt == 1:
        return "SNP"
    else:
        return "INDEL"


def import_VCF_to_pandas(vcf_file):
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            #print(next_line)
            next_line = f.readline()

    if first_line.endswith('VCFv4.1'):
        df = pd.read_csv(vcf_file, sep='\t', skiprows=[header_lines], header=header_lines)
        
        for index, _ in df.iterrows():
            info_fields = re.findall(r';*([a-zA-Z]{1,20})=', df.loc[index,'INFO'])
            info_values = re.sub(r'([a-zA-Z]{1,20})=', '', df.loc[index,'INFO']).split(";") #Remove fields and split the remaining
        
            for ifield, ivalue in zip(info_fields,info_values):
                df.loc[index,ifield] = ivalue
        
        #df = df[(~df['RES'].str.startswith("phylo"))] #Remove phylo(lineage) markers
        df['ALT']=df['ALT'].str.upper()
        df['REF']=df['REF'].str.upper()

        return df
    else:
        print("This vcf file is not v4.1")
        sys.exit(1)



def import_VCF42_to_pandas_legacy(vcf_file, sep='\t'):
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            #print(next_line)
            next_line = f.readline()
    
    if first_line.endswith('VCFv4.2'):
        
        #Use first line as header
        dataframe = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)
        sample = dataframe.columns[-1]
        dataframe.rename(columns={sample:'sample'}, inplace=True)
        
        for index in dataframe.index:
            info_fields = re.findall(r';*([a-zA-Z]{1,20})=', dataframe.loc[index,'INFO'])
            info_values = re.findall(r'-?\d+\.?\d*e?[+-]?\d{0,2}', dataframe.loc[index,'INFO'])
            
            format_fields = dataframe.loc[index,'FORMAT'].split(":")
            format_values = dataframe.loc[index,'sample'].split(":")
                                    
            for ifield, ivalue in zip(info_fields,info_values):
                dataframe.loc[index,ifield] = ivalue
                
            for ffield, fvalue in zip(format_fields,format_values):
                dataframe.loc[index,ffield] = fvalue
            #if len(format_values[1].split(",")) != 2:
            #    print(format_values[1].split(","), index)
            #    print(dataframe.iloc[index])
        dataframe.rename(columns={'AF':'af'}, inplace=True)
        dataframe['REF_AD'] = dataframe['AD'].str.split(",").str[0]
        dataframe['ALT_AD'] = dataframe['AD'].str.split(",").str[1]
        # dataframe['REF_AD'] = dataframe['AD'].str.split(",").str[-2:].str[0] #When there is a minoritary third allele it places third in AD???
        #dataframe['ALT_AD'] = dataframe['AD'].str.split(",").str[-2:].str[1]
        
        to_float = ['QUAL', 'AC', 'af', 'AN', 'BaseQRankSum', 'DP', 'ExcessHet', 'FS',
       'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR','GQ','ALT_AD', 'REF_AD', 'InbreedingCoeff']
        
        to_int = ['POS', 'len_AD', 'gt0', 'gt1']
        
        to_str = ['#CHROM','REF','ALT', 'FILTER']
        
        for column in dataframe.columns:
            if column in to_float:
                dataframe[column] = dataframe[column].astype(float)
                
        for column in dataframe.columns:
            if column in to_int:
                dataframe[column] = dataframe[column].astype(int)
                
        for column in dataframe.columns:
            if column in to_str:
                dataframe[column] = dataframe[column].astype(str)

        dataframe['dp'] = (dataframe['REF_AD'] + dataframe['ALT_AD'])
        dataframe['aF'] = dataframe['REF_AD']/dataframe['dp']
        dataframe['AF'] = dataframe['ALT_AD']/dataframe['dp']

    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe

def add_snp_distance(vcf_df, max_length=False):
    """
    Calculate distance to the closest left and rigth SNP using a vcf imported as datafame
    Total reference length is inferred from vcf by default adding 100bp to the largest position
    in order to avoid reference parse, but it can me supplied by user
    """
    if max_length == False:
        max_length = max(vcf_df.POS.values.tolist()) + 100
        
    for index, _ in vcf_df[vcf_df.TYPE == "SNP"].iterrows():
        if index == 0:
            vcf_df.loc[index,'snp_left_distance'] = vcf_df.loc[index,'POS'] - 0
        elif index > 0:
            vcf_df.loc[index,'snp_left_distance'] = vcf_df.loc[index,'POS'] - vcf_df.loc[index - 1,'POS']
        if index == (len(vcf_df.index.values) - 1):
            vcf_df.loc[index,'snp_right_distance'] = max_length - vcf_df.loc[index,'POS']
        elif index < (len(vcf_df.index.values) - 1):
            vcf_df.loc[index,'snp_right_distance'] = vcf_df.loc[index + 1,'POS'] - vcf_df.loc[index,'POS']
            
    return vcf_df

def add_indel_distance(vcf_df, max_length=False):
    """
    Calculate distance to the closest left and rigth INDEL using a vcf imported as datafame
    Total reference length is inferred from vcf by default adding 100bp to the largest position
    in order to avoid reference parse, but it can me supplied by user
    """
    if max_length == False:
        max_length = max(vcf_df.POS.values.tolist()) + 100
        
    for index, _ in vcf_df[vcf_df.TYPE == "SNP"].iterrows():
        if index > 0 and index < max(vcf_df.index) and (vcf_df.loc[index - 1,'TYPE'] == 'INDEL'):
            if index == 0:
                vcf_df.loc[index,'indel_left_distance'] = vcf_df.loc[index,'POS'] - 0
            elif index > 0:
                vcf_df.loc[index,'indel_left_distance'] = vcf_df.loc[index,'POS'] - vcf_df.loc[index - 1,'POS']
        if index > 0 and index < max(vcf_df.index) and (vcf_df.loc[index + 1,'TYPE'] == 'INDEL'):
            if (index == (len(vcf_df.index.values) - 1)):
                vcf_df.loc[index,'indel_right_distance'] = max_length - vcf_df.loc[index,'POS']
            elif (index < (len(vcf_df.index.values) - 1)):
                vcf_df.loc[index,'indel_right_distance'] = vcf_df.loc[index + 1,'POS'] - vcf_df.loc[index,'POS']

    return vcf_df


def add_window_distance(vcf_df, window_size=10):
    """
    Add a column indicating the maximum number of SNPs in a windows of 10 or supplied distance
    """
    list_pos = vcf_df.POS.to_list() #all positions
    set_pos = set(list_pos) #to set for later comparing
    max_pos = max(vcf_df.POS.to_list()) #max to iter over positions (independent from reference)

    all_list = list(range(1, max_pos + 1)) #create a list to slide one by one
    
    df_header = "window_" + str(window_size)

    vcf_df[df_header] = 1 #Create all 1 by default

    #Slide over windows
    for i in range(0,max_pos,1):
        window_pos = all_list[i:i+window_size] #This splits the list in windows of determined length
        set_window_pos = set(window_pos)
        #How many known positions are in every window for later clasification
        num_conglomerate = set_pos & set_window_pos
        
        if len(num_conglomerate) > 1:
            for i in num_conglomerate:
                index = vcf_df.index[vcf_df["POS"] == i][0] #Retrieve index with the known position
                if vcf_df.loc[index,df_header] < len(num_conglomerate):
                    vcf_df.loc[index,df_header] = len(num_conglomerate)


def bed_to_dict(bed_file):
    dict_range_positions = {}
    with open(bed_file, 'r') as f:
        for line_number, line in enumerate(f):
            line_split = line.split(None) #This split by any blank character
            start = line_split[1]
            end = line_split[2]
            if len(line_split) == 3 and start.isdigit() and end.isdigit():
                start = int(start)
                end = int(end)
                dict_range_positions[start] = end
            else:
                if line_number != 0:
                    print("This file is not in bed format")
                    sys.exit(1)
    return dict_range_positions


def annotate_bed(dict_position, position):
    """
    Identify a position within a range
    credits: https://stackoverflow.com/questions/6053974/python-efficiently-check-if-integer-is-within-many-ranges
    """
    #dict_position = bed_to_dict(bed_file)
    if any(start <= position <= end for (start, end) in dict_position.items()):
        return True
    else:
        return False

def handle_polymorphism_freebayes(df):
    for index, _ in df[df.len_AD > 2].iterrows():
        split_AD_all = df.loc[index, 'AD'].split(",")
        split_AD = split_AD_all[1:]
        split_AD = [int(x) for x in split_AD]
        maxAD = max(split_AD)
        max_index = split_AD.index(maxAD)#Obtain index from highest value in list of positions
        df.loc[index, 'len_AD'] = 2 #reset number of alternatives as normal
        df.loc[index, 'AD'] = split_AD_all[0] + ',' + str(maxAD)
        repare_headers = ['ALT', 'AF', 'LEN', 'TYPE', 'ALT_QUAL', 'ALT_DP']
        for header in repare_headers:        
            df.loc[index, header] = df.loc[index, header].split(",")[max_index]  #split bases into list and retrieve the base using ps index

def import_VCF42_freebayes_to_df(vcf_file, sep='\t'):
    """
    Script to read vcf 4.2
    - now handle correct allele frequency calculated by summing REF reads + ALT reads instead from DP parameter
    - now retrieve the largest read number for ALT allele frequency in case is a heterozygous SNP (depends on calculate_ALT_AD())
    - now uses dataframe.iterrows() instead dataframe.index
    - remove snps with two alternate alleles, keeping the most abundant if this is more at least 3 times more frequent
    """

    header_lines = 0
    if vcf_file.endswith(".gz"):
        with gzip.open(vcf_file, 'rb') as f:
            first_line = f.readline().decode().strip()
            next_line = f.readline().decode().strip()
            while next_line.startswith("##"):
                header_lines = header_lines + 1
                next_line = f.readline().decode().strip()
    else:
        with open(vcf_file, 'r') as f:
            first_line = f.readline().strip()
            next_line = f.readline().strip()
            while next_line.startswith("##"):
                header_lines = header_lines + 1
                next_line = f.readline().strip()
    
    if first_line.endswith('VCFv4.2'):
        
        #Use first line as header
        if vcf_file.endswith(".gz"):
            df = pd.read_csv(vcf_file, compression='gzip', sep=sep, skiprows=[header_lines], header=header_lines)
        else:
            df = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)

        sample = df.columns[-1]
        df.rename(columns={sample:'sample'}, inplace=True)
        
        for index, data_row in df.iterrows():
            info_fields = [x.split('=')[0] for x in data_row.INFO.split(';')]
            info_values = [x.split('=')[1] for x in data_row.INFO.split(';')]
            
            format_fields = data_row['FORMAT'].split(":")
            format_values = data_row['sample'].split(":")
                                    
            for ifield, ivalue in zip(info_fields,info_values):
                df.loc[index,ifield] = ivalue
                
            for ffield, fvalue in zip(format_fields,format_values):
                df.loc[index,ffield] = fvalue

        #Rename columns to match iVar output to adapt covidma flow
        df.rename(columns={'#CHROM':'REGION','RO':'REF_DP', 'DP':'TOTAL_DP', 'AO':'ALT_DP', 'QR':'REF_QUAL', 'QA':'ALT_QUAL'}, inplace=True)
        
        df['len_AD'] = df['AD'].str.split(",").str.len()
        
        # this step remove false snps from cohort calling and reset index
        #df = df[df.ALT_AD > 0].reset_index(drop=True)

        handle_polymorphism_freebayes(df) #Leave the most common variation when teo alternatives are present        

        to_float = ['QUAL', 'AN', 'TOTAL_DP','GQ', 'REF_QUAL', 'ALT_QUAL']
        
        to_int = ['POS', 'REF_DP', 'ALT_DP', 'len_AD' ]
        
        to_str = ['REGION','REF','ALT', 'FILTER']
        
        for column in df.columns:
            if column in to_float:
                df[column] = df[column].astype(float)
                
        for column in df.columns:
            if column in to_int:
                df[column] = df[column].astype(int)
                
        for column in df.columns:
            if column in to_str:
                df[column] = df[column].astype(str)
                
        df['REF_FREQ'] = df['REF_DP']/df['TOTAL_DP']
        df['ALT_FREQ'] = df['ALT_DP']/df['TOTAL_DP']
        
        df = df[df.TYPE != 'complex']

        df = df.sort_values(by=['POS']).reset_index(drop=True)
        
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return df[['REGION', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'AF', 'AN', 'TOTAL_DP', 'DPB', 'LEN', 
       'NUMALT', 'ODDS', 'TYPE', 'GT', 'AD', 'len_AD', 'REF_DP', 'ALT_DP', 'REF_QUAL', 'ALT_QUAL', 'REF_FREQ', 'ALT_FREQ']]



