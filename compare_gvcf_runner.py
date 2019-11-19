#!/usr/bin/env python

import os
import sys
import re
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, execute_subprocess, \
    extract_read_list, file_to_list, get_coverage, obtain_group_cov_stats, remove_low_covered, check_remove_file
from bbduk_trimmer import bbduk_trimming
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_recall import picard_dictionary, samtools_faidx, picard_markdup, haplotype_caller, call_variants, \
    select_variants, hard_filter, combine_gvcf, select_pass, select_pass_variants, recalibrate_bam, \
    split_vcf_saples
from vcf_process import vcf_consensus_filter

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 14 June 2019
REVISION: 

TODO:
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def get_arguments():

    parser = argparse.ArgumentParser(prog = 'snptb.py', description= 'Pipeline to call variants (SNVs) with any non model organism. Specialised in Mycobacterium Tuberculosis')
    
    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', metavar="input_cohort_vcf", type=str, required=True, help='REQUIRED.Input cohort vcf')
    input_group.add_argument('-s', '--sample_list', type=str, required=True, help='Sample names to analyse only in the file supplied')

    output_group = parser.add_argument_group('Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')

    vcf_group = parser.add_argument_group('VCF filters', 'parameters for variant filtering')

    vcf_group.add_argument('-m', '--maxnocallfr', type=str, required=False, default=0.2, help='maximun proportion of samples with non genotyped alleles')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=4, help='Threads to use')
    params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=8, help='MAx memory to use')



    #argcomplete.autocomplete(parser)
    arguments = parser.parse_args()

    return arguments

args = get_arguments()


sample_list_F = file_to_list(args.sample_list)
print("\n%d samples will be analysed: %s" % (len(sample_list_F), ",".join(sample_list_F)))


######################################################################
#####################START PIPELINE###################################
######################################################################
output = os.path.abspath(args.output)
group_name = args.input.split("/")[-1].split(".")[0]
out_vcf_dir = os.path.join(args.output, "VCF")
check_create_dir(out_vcf_dir)

output_vcf_file = os.path.abspath(args.input)

base_input = os.path.basename(args.input)

linked_file = os.path.join(out_vcf_dir, base_input)

check_remove_file(linked_file)

os.symlink(output_vcf_file, linked_file)


print("\n\n" + BLUE + BOLD + "STARTING COHORT GVCF TO SPLIT SAMPLE VCF IN GROUP: " + group_name + END_FORMATTING)


#SELECT VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#########################################################
out_vcfsnp_name = group_name + ".cohort.snp.vcf"
output_vcfsnp_file = os.path.join(out_vcf_dir, out_vcfsnp_name)

if os.path.isfile(output_vcfsnp_file):
    print(YELLOW + DIM + output_vcfsnp_file + " EXIST\nOmmiting Variant Selection (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Selecting Variants (Group) in group " + group_name + END_FORMATTING)
    select_variants(linked_file, select_type='SNP')

#HARD FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnp_name = group_name + ".cohort.snp.hf.vcf"
output_vcfhfsnp_file = os.path.join(out_vcf_dir, out_vcfhfsnp_name)


if os.path.isfile(output_vcfhfsnp_file):
    print(YELLOW + DIM + output_vcfhfsnp_file + " EXIST\nOmmiting Hard Filtering (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Hard Filtering Variants (Group) in group " + group_name + END_FORMATTING)
    hard_filter(output_vcfsnp_file, select_type='SNP')


#PASS FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnppass_name = group_name + ".cohort.snp.hf.pass.vcf"
out_vcfhfindelpass_name = group_name + ".cohort.indel.hf.pass.vcf"
output_vcfhfsnppass_file = os.path.join(out_vcf_dir, out_vcfhfsnppass_name)
output_vcfhfindelpass_file = os.path.join(out_vcf_dir, out_vcfhfindelpass_name)


if os.path.isfile(output_vcfhfindelpass_file) and os.path.isfile(output_vcfhfsnppass_file):
    print(YELLOW + DIM + output_vcfhfsnppass_file + " EXIST\nOmmiting PASS Filtering (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "PASS Filtering Variants (Group) in group " + group_name + END_FORMATTING)
    select_pass_variants(output_vcfhfsnp_file, nocall_fr=0.2)


split_vcf_saples(output_vcfhfsnppass_file, sample_list=sample_list_F)



for sample in sample_list_F:

    print("\n" + WHITE_BG + "FINAL FILTERING IN SAMPLE " + sample + END_FORMATTING)

    ################FINAL VCF FILTERING##################
    #####################################################
    out_final_name = sample + ".snp.hf.pass.final.vcf"
    in_final_name = sample + ".snp.hf.pass.vcf"
    output_final_vcf = os.path.join(out_vcf_dir, out_final_name)
    in_final_vcf = os.path.join(out_vcf_dir, in_final_name)

    if os.path.isfile(output_final_vcf):
        print(YELLOW + DIM + output_final_vcf + " EXIST\nOmmiting Final filter for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Final filter in sample " + sample + END_FORMATTING)
        vcf_consensus_filter(in_final_vcf,  distance=1, AF=0.75, QD=10, window_10=3)


print("\n\n" + MAGENTA + BOLD + "SAMPLE VARIANT CALL FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")

