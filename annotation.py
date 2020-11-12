#!/usr/bin/env python

import os
import sys
import logging
import argparse
import pandas as pd
import numpy as np
import re
import subprocess
#from tabulate import tabulate
from misc import check_create_dir, execute_subprocess

logger = logging.getLogger()


##Import files containing annotation info and convert them to dictionary
#script_dir = os.path.dirname(os.path.realpath(__file__))


def tsv_to_vcf(tsv_file):
    df = pd.read_csv(tsv_file, sep="\t")
    is_empty = df.shape[0] == 0
    #df.insert(2, 'ID', '.')
    df.fillna(".", inplace=True)
    df["PASS"].replace({True: 'PASS'}, inplace=True)
    df.rename(columns={"REGION": "#CHROM", "GFF_FEATURE": "ID", "ALT_QUAL": "QUAL", "PASS": "FILTER"}, inplace=True)

    fial_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO']

    if not is_empty:
        df['INFO'] = df.apply(lambda x: "CODON={}-{};AA={}-{};DP={}".format(x.REF_CODON, x.ALT_CODON, x.REF_AA, x.ALT_AA, x.TOTAL_DP), axis=1)
    else:
        df = df.reindex(columns = fial_columns)
    df = df[fial_columns]
    
    return df

def snpeff_execution(vcf_file, annot_file, database=False):
    df_vcf = pd.read_csv(vcf_file, sep="\t")
    if df_vcf.shape[0] != 0:
        cmd = ["snpEff", "-noStats", database, vcf_file]
        with open(annot_file, "w+") as outfile:
            #calculate coverage and save it in th eoutput file
            subprocess.run(cmd,
            stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)
    else:
        with open(annot_file, "w+") as outfile:
            outfile.write('No annotation found')
            

def import_annot_to_pandas(vcf_file, sep='\t'):
    """
    Order several annoattion by:
    Putative impact: Effects having higher putative impact are first.
    Effect type: Effects assumed to be more deleterious effects first.
    Canonical transcript before non-canonical.
    Marker genomic coordinates (e.g. genes starting before first)
    https://pcingola.github.io/SnpEff/se_inputoutput/
    Parse vcf outputted by snpEFF which adds the ANN field
    Dependences: calculate_ALT_AD
                calculate_true_ALT
    """
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        if first_line == 'No annotation found':
            return pd.read_csv(vcf_file, sep=sep)
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            #print(next_line)
            next_line = f.readline()
        
    #Use first line as header
    df = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)

    ann_headers = ['Allele',
                    'Annotation',
                    'Annotation_Impact',
                    'Gene_Name',
                    'Gene_ID',
                    'Feature_Type',
                    'Feature_ID',
                    'Transcript_BioType',
                    'Rank',
                    'HGVS.c',
                    'HGVS.p',
                    'cDNA.pos / cDNA.length',
                    'CDS.pos / CDS.length',
                    'AA.pos / AA.length',
                    'ERRORS / WARNINGS / INFO']
    anlelle_headers = ['Codon_change', 'AA_change', 'DP', 'Allele']

    #Apply function to split and recover the first 15 fields = only first anotations, the most likely

    df['TMP_ANN_16'] = df['INFO'].apply(lambda x: ('|').join(x.split('|')[0:15]))
    df[ann_headers] = df['TMP_ANN_16'].str.split('|', expand=True)
    df['HGVS.c'] = df['HGVS.c'].str.split(".").str[-1]
    df['HGVS.p'] = df['HGVS.p'].str.split(".").str[-1]
    df[anlelle_headers] = df['Allele'].str.split(';', expand=True)
    
    for head in anlelle_headers:
        df[head] = df[head].str.split("=").str[-1]

    del df['TMP_ANN_16']

    #Send INFO column to last position
    df = df[ [ col for col in df.columns if col != 'INFO' ] + ['INFO'] ]

    return df

def annotate_snpeff(input_tsv_file, output_vcf_file, output_annot_file, database='NC_045512.2'):
    vcf_df = tsv_to_vcf(input_tsv_file)
    vcf_df.to_csv(output_vcf_file, sep="\t", index=False)
    #Execure snpEff
    snpeff_execution(output_vcf_file, output_annot_file, database=database)
    #Format annot vcf and remove vcf
    annot_df = import_annot_to_pandas(output_annot_file)
    annot_df.to_csv(output_annot_file, sep="\t", index=False)
    os.remove(output_vcf_file)


def annotate_pangolin(input_file, output_folder, output_filename, threads=8, max_ambig=0.6):
    cmd = ["pangolin", input_file, "--outdir", output_folder, "--outfile", output_filename, "--threads", str(threads), "--max-ambig", str(max_ambig)]
    execute_subprocess(cmd)

def get_reverse(nucleotyde):
    nucleotyde = str(nucleotyde)
    nucleotyde_rev = {'A' : 'T',
                     'T' : 'A',
                     'C' : 'G',
                     'G': 'C'}
    if len(nucleotyde) > 1:
        nucleotyde_str = nucleotyde[::-1] #Reverse nucleotide
        nucleotyde_str_fin = "".join([nucleotyde_rev[x] for x in nucleotyde_str]) #Complement nucleotide
        return nucleotyde_str_fin
    else:
        return nucleotyde_rev[nucleotyde]

def import_VCF_to_pandas(vcf_file):
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            #print(next_line)
            next_line = f.readline()

    if first_line.startswith('##'):
        df = pd.read_csv(vcf_file, sep='\t', skiprows=[header_lines], header=header_lines)
        
        df['ALT']=df['ALT'].str.upper()
        df['REF']=df['REF'].str.upper()
        #Check INFO
        if 'INFO' in df.columns:
            return df
        else:
            last_column = df.columns[-1]
            df = df.rename(columns={last_column: 'INFO'})
            return df
    else:
        print("This vcf file is not properly formatted")
        sys.exit(1)

def annotate_vcfs(tsv_df, vcfs):
    df = pd.read_csv(tsv_df, sep="\t")
    for vcf in vcfs:
        print("ANNOTATING VCF: ", vcf)
        header = (".").join(vcf.split("/")[-1].split(".")[0:-1])
        dfvcf = import_VCF_to_pandas(vcf)
        dfvcf = dfvcf[['POS', 'REF', 'ALT', 'INFO']]
        dfvcf = dfvcf.rename(columns={'INFO': header})
        df = df.merge(dfvcf, how='left')
    return df

def bed_to_df(bed_file):
    """
    Import bed file separated by tabs into a pandas df
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
        df = pd.read_csv(bed_file, sep="\t", header=None) #delim_whitespace=True
    else:
        df = pd.read_csv(bed_file, sep="\t", skiprows=header_lines, header=None) #delim_whitespace=True

    df = df.iloc[:,0:4]
    df.columns = ["#CHROM", "start", "end", "description"]
        
    return df

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
        return None

def annotate_bed_s(tsv_df, bed_files):
    df = pd.read_csv(tsv_df, sep="\t")

    variable_list = [ x.split("/")[-1].split(".")[0] for x in bed_files] #extract file name and use it as header
    
    for variable_name, bed_file in zip(variable_list,bed_files):
        print("ANNOTATING BED: ", bed_file)
        bed_annot_df = bed_to_df(bed_file)
        df[variable_name] = df['POS'].apply(lambda x: add_bed_info(bed_annot_df,x))
    return df

def user_annotation(tsv_file, output_file, vcf_files=[], bed_files=[]):
    bed_df = annotate_bed_s(tsv_file, bed_files)
    vcf_df = annotate_vcfs(tsv_file, vcf_files)

    df = bed_df.merge(vcf_df)

    df.to_csv(output_file, sep="\t", index=False)


if __name__ == '__main__':
    print("#################### ANNOTATION #########################")