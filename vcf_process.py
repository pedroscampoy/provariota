#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np
import re
import gzip
import subprocess
from misc import check_file_exists, obtain_output_dir, check_create_dir, get_picard_path, execute_subprocess, check_remove_file, \
list_to_bed, count_lines


def calculate_ALT_AD(row):
    split_AD = row.AD.split(",")[1:]
    split_AD = [int(x) for x in split_AD]
    if row.len_AD > 2:
        max_AD = max(split_AD)
        #max_index = split_AD.index(max(split_AD))
        return max_AD
    else:
        #ALT_AD = row.AD.split(",")[1]
        return split_AD[0]

def calculate_true_ALT(row):
    split_AD = row.AD.split(",")[1:]
    split_AD = [int(x) for x in split_AD]
    split_ALT = row.ALT.split(",")
    if row.len_AD > 2:
        max_index = split_AD.index(max(split_AD))
        true_ALT = split_ALT[max_index]
        return true_ALT
    else:
        #ALT_AD = row.AD.split(",")[1]
        return split_ALT[0]

def handle_polymorphism(vcf_df):
    for index, _ in vcf_df[vcf_df.len_AD > 2].iterrows():
        split_AD = vcf_df.loc[index, 'AD'].split(",")[1:]
        split_AD = [int(x) for x in split_AD]
        if min(split_AD) == 0:
            max_index = split_AD.index(max(split_AD)) #Obtain index from highest value in list of positions
            split_ALT = vcf_df.loc[index, 'ALT'].split(",") #split bases into list
            vcf_df.loc[index, 'ALT'] = split_ALT[max_index] #Retrieve the base using ps index
            vcf_df.loc[index, 'len_AD'] = 2 #reset number of alternatives as normal
        elif max(split_AD)/min(split_AD) > 3: # the same process to avoid dividing by 0
            max_index = split_AD.index(max(split_AD))
            split_ALT = vcf_df.loc[index, 'ALT'].split(",")
            vcf_df.loc[index, 'ALT'] = split_ALT[max_index]
            vcf_df.loc[index, 'len_AD'] = 2

def define_var_type(row):
    len_ref = len(row.REF)
    len_alt = len(row.ALT)
    
    if len_ref == len_alt == 1:
        return "SNP"
    else:
        return "INDEL"

def import_VCF42_to_pandas(vcf_file, sep='\t'):
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
            dataframe = pd.read_csv(vcf_file, compression='gzip', sep=sep, skiprows=[header_lines], header=header_lines)
        else:
            dataframe = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)

        sample = dataframe.columns[-1]
        dataframe.rename(columns={sample:'sample'}, inplace=True)
        
        for index, data_row in dataframe.iterrows():
            info_fields = re.findall(r';*([a-zA-Z]{1,20})=', data_row.INFO)
            info_values = re.findall(r'-?\d+\.?\d*e?[+-]?\d{0,2}', data_row.INFO)
            
            format_fields = data_row['FORMAT'].split(":")
            format_values = data_row['sample'].split(":")
                                    
            for ifield, ivalue in zip(info_fields,info_values):
                dataframe.loc[index,ifield] = ivalue
                
            for ffield, fvalue in zip(format_fields,format_values):
                dataframe.loc[index,ffield] = fvalue
            
        dataframe.rename(columns={'AF':'af'}, inplace=True)
        
        dataframe['len_AD'] = dataframe['AD'].str.split(",").str.len()
        dataframe['REF_AD'] = dataframe['AD'].str.split(",").str[0]

        dataframe['ALT_AD'] = dataframe.apply(calculate_ALT_AD, axis=1)
        dataframe[['gt0','gt1']] = dataframe['GT'].str.split(r'[/|\|]', expand=True)
        
        # this step remove false snps from cohort calling and reset index
        #dataframe = dataframe[dataframe.ALT_AD > 0].reset_index(drop=True)

        handle_polymorphism(dataframe) #Leave the most common variation
        dataframe['TYPE'] = dataframe.apply(define_var_type, axis=1)
        

        to_float = ['QUAL', 'AC', 'af', 'AN', 'BaseQRankSum', 'DP', 'ExcessHet', 'FS',
       'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR','GQ','ALT_AD', 'REF_AD']
        
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
        
        dataframe = dataframe.sort_values(by=['POS']).reset_index(drop=True)
        
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe

def import_VCF41_to_pandas(vcf_file):
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
        #df[['Gene ID', 'Gene name', 'Gene start', 'Gene stop']] = df.GENE.str.split(":", expand=True)
        #df['GENE'] = df['INFO'].apply(lambda x: extract_gene_name(x))
        #df['Isreverse'] = df['GENE'].apply(lambda x: True if x.endswith("c") else False)

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

def add_window_distance_legacy(vcf_df, window_size=10):
    """
    DEPRECATED
    Add a column indicating the maximum number of SNPs in a windows of 10
    """
    list_pos = vcf_df.POS.to_list() #all positions
    set_pos = set(list_pos) #to set for later comparing
    max_pos = max(vcf_df.POS.to_list()) #max to iter over positions (independent from reference)

    all_list = list(range(1, max_pos + 1)) #create a list to slide one by one
    df_header = "window_" + str(window_size)

    #Create sets
    set_2 = set()
    set_3 = set()
    set_4 = set()
    set_5 = set()
    
    sets = [set_2, set_3, set_4, set_5]
    
    #Slide over windows
    for i in range(0,max_pos,1):
        window_pos = all_list[i:i+window_size] #This splits the list in windows of determined length
        set_window_pos = set(window_pos)
        #How many known positions are in every window for later clasification
        num_conglomerate = set_pos & set_window_pos
        
        if len(num_conglomerate) > 4:
            set_5.update(num_conglomerate)
        elif len(num_conglomerate) == 4:
            set_4.update(num_conglomerate)
        elif len(num_conglomerate) == 3:
            set_3.update(num_conglomerate)
        elif len(num_conglomerate) == 2:
            set_2.update(num_conglomerate)
    #Remove positions in a higher number of sets
    for set_num in range(0, len(sets)):
        if set_num < (len(sets) - 1):
            sets[set_num] = sets[set_num] - sets[set_num + 1]

    for index, _ in vcf_df.iterrows():
        if vcf_df.loc[index,'POS'] in sets[0]:
            vcf_df.loc[index, df_header] = 2
        elif vcf_df.loc[index,'POS'] in sets[1]:
            vcf_df.loc[index, df_header] = 3
        elif vcf_df.loc[index,'POS'] in sets[2]:
            vcf_df.loc[index, df_header] = 4
        elif vcf_df.loc[index,'POS'] in sets[3]:
            vcf_df.loc[index, df_header] = 5
        else:
            vcf_df.loc[index, df_header] = 1
            
    vcf_df[df_header] = vcf_df[df_header].astype(int)

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


def bed_to_df(bed_file):
    """
    Import bed file separated by tabs into a pandas dataframe
    -Handle header line
    -Handle with and without description (If there is no description adds true or false to annotated df)
    """
    header_lines = 0
    #Handle likely header by checking colums 2 and 3 as numbers
    with open(bed_file, 'r') as f:
        next_line = f.readline().strip()
        line_split = next_line.split(None) #This split by any blank character
        start = line_split[1]
        end = line_split[2]
        while not start.isdigit() and not end.isdigit():
            header_lines = header_lines + 1
            next_line = f.readline().strip()
            line_split = next_line.split(None) #This split by any blank character
            start = line_split[1]
            end = line_split[2]

    if header_lines == 0:
        dataframe = pd.read_csv(bed_file, sep="\t", header=None) #delim_whitespace=True
    else:
        dataframe = pd.read_csv(bed_file, sep="\t", skiprows=header_lines, header=None) #delim_whitespace=True
    if dataframe.shape[1] == 3:
        dataframe['description'] = True
        dataframe.columns = ["#CHROM", "start", "end", "description"]
    else:
        dataframe.columns = ["#CHROM", "start", "end", "description"]
        
    return dataframe

def add_bed_info(bed_df, position):
    """
    Identify a position within a range
    credits: https://stackoverflow.com/questions/6053974/python-efficiently-check-if-integer-is-within-many-ranges
    """
    #dict_position = bed_to_dict(bed_file)
    if any(start <= position <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())):
        description_out = bed_df.description[(bed_df.start <= position) & (bed_df.end >= position)].values[0]
        return description_out
    else:
        return False

def annotate_bed_s(vcf_annot, *bed_files):
    """
    More on list comprehension: https://stackoverflow.com/questions/3371269/call-int-function-on-every-list-element
    """
    print("ANNOTATING BED(S): ", bed_files)
    #bed_files = [ os.path.abspath(x) for x in bed_files ]
    #bed_files = list(map(os.path.abspath, bed_files)) #get full path for all files
    variable_list = [ x.split("/")[-1].split(".")[0] for x in bed_files ] #extract file name and use it as header
    
    for variable_name, bed_file in zip(variable_list,bed_files):
        bed_annot_df = bed_to_df(bed_file)
        vcf_annot[variable_name] = vcf_annot['POS'].apply(lambda x: add_bed_info(bed_annot_df,x))


def filter_vcf_list(raw_vcf, list_pos, name_out):
    """
    Remove positions supplies in a list and creates a different vcf with the name supplied
    """ 
    input_vcf = os.path.abspath(raw_vcf)
    input_vcf_dir_name = ("/").join(input_vcf.split("/")[:-1])
     
    vcf_output_file = os.path.join(input_vcf_dir_name, name_out)
        
    with open(input_vcf, "r") as f:
        with open(vcf_output_file, "w") as f1:
            for line in f:
                if line.startswith("#"):
                    f1.write(line)
                else:
                    position =  int(line.split("\t")[1])
                    if position not in list_pos:
                        f1.write(line)

def import_VCF42_cohort_pandas(vcf_file, sep='\t'):
    """
    Script to read vcf 4.2 cohort/join called vcf handling header lines
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
        dataframe = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)

    return dataframe

def identify_heterozygous(vcf_file, nocall_fr=0.5):
    
    df = import_VCF42_cohort_pandas(vcf_file)
    
    highly_hetz_positions = []

    #locate heterozygous positions in 0.2% or more samples
    for index, data_row in df.iloc[:,9:].iterrows():
        if any(bool(re.search(r'0[|\/][1-9]', x)) for x in data_row):
            #print(data_row.tolist())
            is_heterozygous = [bool(re.search(r'0[|\/][1-9]', x)) for x in data_row] #True False array
            #is_heterozygous = [x.startswith("0/1") for x in data_row] #True False array
            is_heterozygous_count = sum(is_heterozygous) #True = 1, False = 0
            #Drop positions
            if is_heterozygous_count / len(is_heterozygous) >= nocall_fr:
                highly_hetz_positions.append(df.loc[index, 'POS'])
                #print(df.loc[index, 'POS'], is_heterozygous_count, len(is_heterozygous))
        
    return highly_hetz_positions

def highly_hetz_to_bed(cohort_vcf_file, output_file_name, reference="CHROM", nocall_fr=0.5):
    """
    Determine positions with heterozygous positions from a cohort vcf and
    create a bed named highly_polimorphic.bed for further annotation
    """
    #Set output dir as input
    output_dir = ("/").join(cohort_vcf_file.split("/")[:-1])
    list_heterozygous = identify_heterozygous(cohort_vcf_file, nocall_fr)
    list_to_bed(list_heterozygous, output_dir, output_file_name, reference)

def identify_non_genotyped(vcf_file, nocall_fr=0.5):
    
    df = import_VCF42_cohort_pandas(vcf_file)
    
    non_genotyped_positions = []
    
    #sample_list = df.columns[9:].tolist()
    #remove positions which haven't been enotyped in 0.2% or more samples
    for index, data_row in df.iloc[:,9:].iterrows():
        if any(x.startswith("./.") for x in data_row):
            #print(data_row.tolist())
            non_genotyped = [x.startswith("./.") for x in data_row] #True False array
            non_genotyped_count = sum(non_genotyped) #True = 1, False = 0
            if non_genotyped_count / len(non_genotyped) >= nocall_fr:
                non_genotyped_positions.append(df.loc[index, 'POS'])
                #print(df.loc[index, 'POS'], is_heterozygous_count, len(is_heterozygous))
        
    return non_genotyped_positions

def non_genotyped_to_bed(cohort_vcf_file, output_file_name, reference="CHROM", nocall_fr=0.5):
    """
    Determine positions with non genotyped positions from a cohort vcf and
    create a bed named non_genotyped.bed for further annotation
    """
    #Set output dir as input
    output_dir = ("/").join(cohort_vcf_file.split("/")[:-1])
    list_non_genotyped = identify_non_genotyped(cohort_vcf_file, nocall_fr)
    list_to_bed(list_non_genotyped, output_dir, output_file_name, reference)


def coverage_to_list(input_file):
    sample_name = input_file.split("/")[-1].split(".")[0]
    coverage_list = []
    with open(input_file, 'r') as f:
            content = f.read()
            content_list = content.split('\n')
            while '' in content_list : content_list.remove('')
    coverage_list = [x.split("\t")[2] for x in content_list]

    return (sample_name, coverage_list)

def coverage_to_df(input_file,min_coverage=2, nocall_fr=0.5):
    sample_name = input_file.split("/")[-1].split(".")[0]
    min_cov_df = pd.DataFrame()
    coverage_list = []
    with open(input_file, 'r') as f:
            content = f.read()
            content_list = content.split('\n')
            while '' in content_list : content_list.remove('')
    coverage_list = [int(x.split("\t")[2]) for x in content_list]
    min_cov_df[sample_name] = coverage_list
    min_cov_df = min_cov_df[min_cov_df < min_coverage].dropna(how='all')
    
    return min_cov_df

def identify_uncovered(cov_folder, min_coverage=2, nocall_fr=0.5):
    cov_folder = os.path.abspath(cov_folder)
    len_files = set()
    #Create Position column and asign value
    cov_df = pd.DataFrame()
    
    for root, _, files in os.walk(cov_folder):
        for name in files:
            if name.endswith(".cov"):
                filename = os.path.join(root, name)
                #Add number of lines of file and compare to previous file
                len_files.add(count_lines(filename))
                #import to dataframe if they have the same positios(same reference)
                if len(len_files) == 1:
                    low_coverage_df = coverage_to_df(filename)
                    cov_df = cov_df.merge(low_coverage_df, how='outer', left_index=True, right_index=True)          
                else:
                    print("This file has different reference, please, check " + filename)
                    sys.exit(1)
                                   
    #Determine low covered positions in dataframe
    #Filter positions with values lower than min_cov, dro rows with all false and extract the indet to iterate
    df_any_uncovered = cov_df[cov_df < min_coverage].dropna(how='all')#.index.tolist()
    df_any_uncovered['N_uncovered'] = df_any_uncovered.count(axis=1)
    df_any_uncovered['Position'] = df_any_uncovered.index + 1
    
    n_samples = len(df_any_uncovered.columns) - 2
    
    df_half_uncovered_list = df_any_uncovered['Position'][df_any_uncovered.N_uncovered / n_samples >= nocall_fr].tolist()
    
    return df_half_uncovered_list


def poorly_covered_to_bed(coverage_folder, output_file_name, reference="CHROM", min_coverage=2, nocall_fr=0.5):
    """
    Determine shared low covered positions from their bedtools coverage files
    create a bed named poorly_covered.bed for further annotation
    """
    #Set output dir as input
    list_uncovered = identify_uncovered(coverage_folder, min_coverage=2, nocall_fr=nocall_fr)
    list_to_bed(list_uncovered, coverage_folder, output_file_name, reference)

def vcf_consensus_filter(vcf_file, distance=1, AF=0.75, QD=15, window_10=3, dp_limit=8, dp_AF=10, AF_dp=0.80, 
    highly_hetz=False, non_genotyped=False, poorly_covered=False, bed_to_filter=False, var_type="SNP"):
    """
    Apply custom filter to individual vcf based on:
    AF
    snp distance --> Replaced by window_10
    QD
    Window_10, 20 and 30
    gatk asigned genotype for diploid calls
    Highly heterozygous positions 
    Poorly covered positions
    """
    df_vcf = import_VCF42_to_pandas(vcf_file)

    vcf_path = os.path.abspath(vcf_file)
    output_dir = ("/").join(vcf_path.split("/")[:-2])
    vcf_name = vcf_path.split("/")[-1]

    tab_name = (".").join(vcf_name.split(".")[:-1])
    extend_raw = ".raw.tab"
    extend_final = "." + var_type + ".final.vcf"

    table_outputt_dir = os.path.join(output_dir, "Table")
    check_create_dir(table_outputt_dir)

    #Add polymorphic regions info (Phage, Transposon or PE/PPE regions for TB)
    if bed_to_filter == False:
        df_vcf['is_polymorphic'] = False
    else:
        annotate_bed_s(df_vcf, bed_to_filter)
    
    if highly_hetz != False:
        annotate_bed_s(df_vcf, highly_hetz)
    
    if non_genotyped != False:
        annotate_bed_s(df_vcf, non_genotyped)

    if poorly_covered != False:
        annotate_bed_s(df_vcf, poorly_covered)
    

    #Add info of nearby positions
    add_snp_distance(df_vcf)
    add_indel_distance(df_vcf)

    #Add info of clustered positions in sliding window
    add_window_distance(df_vcf, window_size=10)
    add_window_distance(df_vcf, window_size=20)
    add_window_distance(df_vcf, window_size=30)

    #Manage SNP INDEL filter
    if var_type == "SNP":
        var_to_filter = "INDEL"
    elif var_type == "INDEL":
        var_to_filter = "SNP"
    elif var_type == "ALL":
        var_to_filter = "*"
    else:
        print("Wrong variant type to filter, use SNP/INDEL/ALL")
        sys.exit(1)

    #output all raw info into a file in 'Table' folder
    new_out_file = tab_name + extend_raw
    output_raw_tab = os.path.join(table_outputt_dir, new_out_file)
    df_vcf.to_csv(output_raw_tab, sep='\t', index=False)
    
    #Apply all filters and extract positions as table to filer the final vcf
    list_positions_to_filter = df_vcf['POS'][((df_vcf.AF < AF) | 
                                (df_vcf.snp_left_distance <= distance)|
                                (df_vcf.snp_right_distance <= distance)|
                                (df_vcf.window_10 > window_10)|
                                (df_vcf.AF <= 0.0)|
                                (df_vcf.QD <= QD)|
                                (df_vcf.dp == 0)|
                                (df_vcf.len_AD > 2) |
                                (df_vcf.ALT_AD < 2) |
                                (df_vcf.ALT == '*') |
                                (df_vcf.TYPE == var_to_filter) |
                                (df_vcf.dp < dp_limit) |
                                (df_vcf.FILTER != "PASS") |
                                ((df_vcf.gt0 == 0) & (df_vcf.window_10 > 1)) |
                                ((df_vcf.gt0 == 0) & (df_vcf.window_20 >= 2)) |
                                ((df_vcf.gt0 == 0) & (df_vcf.window_30 >= 3)) |
                                ((df_vcf.dp < dp_AF) & (df_vcf.AF < AF_dp)) |
                                (df_vcf.highly_hetz == True) |
                                (df_vcf.poorly_covered == True) |
                                (df_vcf.non_genotyped == True) |
                                (df_vcf.is_polymorphic == True))].tolist()

    final_vcf_name = tab_name + extend_final
    filter_vcf_list(vcf_path, list_positions_to_filter, final_vcf_name)
