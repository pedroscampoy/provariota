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
    samples_from_vcf, split_vcf_saples, combine_gvcf_folder, combine_vcf
from vcf_process import vcf_consensus_filter
from annotation import replace_reference, snpeff_annotation, final_annotation, create_report

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 05 July 2019
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
UNDERLINE = '\033[4m'
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

    input_group.add_argument('-i', '--input', metavar="input_cohort_vcf", type=str, required=True, help='REQUIRED.Input folder with all g.vcf')
    input_group.add_argument('-S', '--sample_list', type=str, required=False, help='Sample names to analyse only in the file supplied')
    input_group.add_argument('-r', '--reference', metavar="reference", type=str, required=True, help='REQUIRED. File to map against')
    input_group.add_argument('-s', '--sample', metavar="sample", type=str, required=False, default=False, help='Sample to identify further files')

    output_group = parser.add_argument_group('Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')

    vcf_group = parser.add_argument_group('VCF filters', 'parameters for variant filtering')
    
    vcf_group.add_argument('-b', '--bed_remove', type=str, required=False, default="TB", help='BED file with position ranges to filter from final vcf')
    vcf_group.add_argument('-m', '--maxnocallfr', type=str, required=False, default=0.1, help='maximun proportion of samples with non genotyped alleles')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=16, help='Threads to use')
    params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=32, help='MAx memory to use')



    #argcomplete.autocomplete(parser)
    arguments = parser.parse_args()

    return arguments

args = get_arguments()

print(args)

######################################################################
#####################START PIPELINE###################################
######################################################################
#Annotation related
script_dir = os.path.dirname(os.path.realpath(__file__))
annotation_dir = os.path.join(script_dir, "annotation/genes")
if args.bed_remove == "TB":
    bed_polymorphism = os.path.join(annotation_dir, "MTB_repeats_annot.bed")


output = os.path.abspath(args.output)
#input_dir = os.path.abspath(args.input)
group_name = output.split("/")[-1]
out_gvcf_dir = os.path.join(args.output, "GVCF")
out_vcf_dir = os.path.join(args.output, "VCF")
check_create_dir(out_vcf_dir)

gvcf_input_dir = os.path.abspath(args.input)



print("\n\n" + BLUE + BOLD + "STARTING COHORT GVCF TO SPLIT SAMPLE VCF IN GROUP: " + group_name + END_FORMATTING)

#CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_gvcf_name = group_name + ".cohort.g.vcf"
output_gvcf_file = os.path.join(out_gvcf_dir, out_gvcf_name)

if os.path.isfile(output_gvcf_file):
    print(YELLOW + DIM + output_gvcf_file + " EXIST\nOmmiting GVCF Combination for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "GVCF Combination in group " + group_name + END_FORMATTING)
    combine_gvcf_folder(args, gvcf_input_dir, sample_list=False)

#CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_vcf_name = group_name + ".cohort.raw.vcf"
output_vcf_file = os.path.join(out_vcf_dir, out_vcf_name)

if os.path.isfile(output_vcf_file):
    print(YELLOW + DIM + output_vcf_file + " EXIST\nOmmiting Variant Calling (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Variant Calling (Group) in group " + group_name + END_FORMATTING)
    call_variants(args, recalibrate=False, group=True)

#SELECT VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#########################################################
out_vcfsnp_name = group_name + ".cohort.snp.vcf"
out_vcfindel_name = group_name + ".cohort.indel.vcf"
output_vcfsnp_file = os.path.join(out_vcf_dir, out_vcfsnp_name)
output_vcfindel_file = os.path.join(out_vcf_dir, out_vcfindel_name)

if os.path.isfile(output_vcfsnp_file) and os.path.isfile(output_vcfindel_file):
    print(YELLOW + DIM + output_vcfsnp_file + " EXIST\nOmmiting Variant Selection (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Selecting Variants (Group) in group " + group_name + END_FORMATTING)
    select_variants(output_vcf_file, select_type='SNP')
    select_variants(output_vcf_file, select_type='INDEL')

#HARD FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnp_name = group_name + ".cohort.snp.hf.vcf"
out_vcfhfindel_name = group_name + ".cohort.indel.hf.vcf"
output_vcfhfsnp_file = os.path.join(out_vcf_dir, out_vcfhfsnp_name)
output_vcfhfindel_file = os.path.join(out_vcf_dir, out_vcfhfindel_name)


if os.path.isfile(output_vcfhfsnp_file) and os.path.isfile(output_vcfhfindel_file):
    print(YELLOW + DIM + output_vcfhfsnp_file + " EXIST\nOmmiting Hard Filtering (Group) for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Hard Filtering Variants (Group) in group " + group_name + END_FORMATTING)
    hard_filter(output_vcfsnp_file, select_type='SNP')
    hard_filter(output_vcfindel_file, select_type='INDEL')



#PASS FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfcombined_name = group_name + ".cohort.combined.hf.vcf"
output_vcfhfcombined_file = os.path.join(out_vcf_dir, out_vcfhfcombined_name)


if os.path.isfile(output_vcfhfcombined_file):
    print(YELLOW + DIM + output_vcfhfcombined_file + " EXIST\nOmmiting combination for group " + group_name + END_FORMATTING)
else:
    print(GREEN + "Combining both vcf SNP and INDEL in group " + group_name + END_FORMATTING)
    combine_vcf(output_vcfhfsnp_file, output_vcfhfindel_file, name_out=False)


split_vcf_saples(output_vcfhfcombined_file, sample_list=False, nocall_fr=args.maxnocallfr)



sample_list_F = samples_from_vcf(output_vcfhfcombined_file)

for sample in sample_list_F:

    print("\n" + WHITE_BG + "FINAL FILTERING IN SAMPLE " + sample + END_FORMATTING)

    ################FINAL VCF FILTERING##################
    #####################################################
    out_final_name = sample + ".combined.hf.SNP.final.vcf"
    in_final_name = sample + ".combined.hf.vcf"
    output_final_vcf = os.path.join(out_vcf_dir, out_final_name)
    in_final_vcf = os.path.join(out_vcf_dir, in_final_name)

    if os.path.isfile(output_final_vcf):
        print(YELLOW + DIM + output_final_vcf + " EXIST\nOmmiting Final filter for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Final filter in sample " + sample + END_FORMATTING)
        vcf_consensus_filter(in_final_vcf, distance=1, AF=0.75, QD=15, window_10=3, dp_limit=8, dp_AF=10, AF_dp=0.80, bed_to_filter=bed_polymorphism, var_type="SNP")


print("\n\n" + MAGENTA + BOLD + "SAMPLE VARIANT CALL FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")