#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging

# Third party imports
import argparse
import subprocess
import datetime


# Local application imports
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, execute_subprocess, \
    extract_read_list, file_to_list, get_coverage, obtain_group_cov_stats, remove_low_covered_mixed, clean_unwanted_files, \
    check_reanalysis, vcf_stats
from bbduk_trimmer import bbduk_trimming
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_recall import picard_dictionary, samtools_faidx, picard_markdup, haplotype_caller, call_variants, \
    select_variants, hard_filter, combine_gvcf, select_pass, select_pass_variants, recalibrate_bam, \
    samples_from_vcf,split_vcf_saples, combine_vcf
from vcf_process import vcf_consensus_filter, highly_hetz_to_bed, poorly_covered_to_bed, non_genotyped_to_bed
from annotation import replace_reference, snpeff_annotation, final_annotation, create_report, css_report
from species_determination import mash_screen, extract_species_from_screen
from compare_snp import ddtb_add, ddtb_compare

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM + CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 28 April 2019
REVISION: 
        20191119 - Adapt to other samples

TODO:
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

#COLORS AND AND FORMATTING

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



#ARGUMENTS

def get_arguments():

    parser = argparse.ArgumentParser(prog = 'snptb.py', description= 'Pipeline to call variants (SNVs) with any non model organism. Specialised in Mycobacterium Tuberculosis')
    
    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all fast[aq] files')
    input_group.add_argument('-r', '--reference', metavar="reference", type=str, required=True, help='REQUIRED. File to map against')
    input_group.add_argument('-s', '--sample', metavar="sample", type=str, required=False, help='Sample to identify further files')
    input_group.add_argument('-S', '--sample_list', type=str, required=False, help='Sample names to analyse only in the file supplied')
    input_group.add_argument('-B', '--annot_bed', type=str, required=False, action='append', help='bed file to annotate')
    input_group.add_argument('-V', '--annot_vcf', type=str, required=False, action='append', help='vcf file to annotate')
    
    output_group = parser.add_argument_group('Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')
    output_group.add_argument('-C', '--noclean', required=False, action='store_false', help='Clean unwanted files for standard execution')

    trimming_group = parser.add_argument_group('Trimming parameters', 'parameters for diferent triming conditions')

    trimming_group.add_argument('-H', '--hdist', type=str, required=False, help='Set hdist parameter, default 2')
    trimming_group.add_argument('-k', '--kmer', type=str, required=False, help='Set k parameter, default 21')

    gatk_group = parser.add_argument_group('GATK parameters', 'parameters for diferent variant calling')

    gatk_group.add_argument('-E', '--enrich_gvcf', required=False,  default=False, help='Point a directory with g.vcf files to enrich the analysis')
    gatk_group.add_argument('-A', '--all_cohort', required=False,  action='store_true', help='Output vcf of all samples instead of just the one inputted before cohort')

    vcf_group = parser.add_argument_group('VCF filters', 'parameters for variant filtering')

    vcf_group.add_argument('-b', '--bed_remove', type=str, required=False, default=False, help='BED file with position ranges to filter from final vcf')
    vcf_group.add_argument('-m', '--maxnocallfr', type=str, required=False, default=0.1, help='maximun proportion of samples with non genotyped alleles')

    annot_group = parser.add_argument_group('Annotation', 'parameters for variant annotation')

    annot_group.add_argument('--mash_database', type=str, required=False, default="/home/pjsola/REFERENCES/mash/refseq.genomes.k21s1000.msh", help='MASH ncbi annotation containing all species database')
    annot_group.add_argument('--snpeff_database', type=str, required=False, default=False, help='snpEFF annotation database')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    
    params_group.add_argument('-c', '--mincov', type=int, required=False, default=20, help='Minimun coverage to add samples into analysis')
    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=16, help='Threads to use')
    params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=32, help='Max memory to use')


    arguments = parser.parse_args()

    return arguments

args = get_arguments()


######################################################################
#####################START PIPELINE###################################
######################################################################
output = os.path.abspath(args.output)
group_name = output.split("/")[-1]

#LOGGING
#Create log file with date and time
right_now = str(datetime.datetime.now())
right_now_full = "_".join(right_now.split(" "))
log_filename = right_now_full + ".log"
log_full_path = os.path.join(output, log_filename)


logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(message)s')

file_handler = logging.FileHandler(log_full_path)
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
#stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)
logger.addHandler(file_handler)



logger.info("\n\n" + BLUE + BOLD + "STARTING PIPELINE IN GROUP: " + group_name + END_FORMATTING)

today = str(datetime.date.today())

logger.info("ARGUMENTS:")
logger.info(str(args))

check_reanalysis(args.output)

#Obtain all R1 and R2 from folder
r1, r2 = extract_read_list(args.input_dir)

#Check if there are samples to filter
sample_list_F = []
if args.sample_list == None:
    logger.info("\n" + "No samples to filter")
    for r1_file, r2_file in zip(r1, r2):
        sample = extract_sample(r1_file, r2_file)
        sample_list_F.append(sample)
else:
    logger.info("samples will be filtered")
    sample_list_F = file_to_list(args.sample_list)
logger.info("\n%d samples will be analysed: %s" % (len(sample_list_F), ",".join(sample_list_F)))




#PREPARE REFERENCE FOR MAPPING + FAI + DICT #########
#####################################################

picard_dictionary(args)
samtools_faidx(args)

#DECLARE FOLDERS CREATED IN PIPELINE ################
#AND KEY FILES ######################################
#####################################################
#Annotation related parameters
script_dir = os.path.dirname(os.path.realpath(__file__))

#Output related
out_trim_dir = os.path.join(output, "Trimmed")
out_map_dir = os.path.join(output, "Bam")
out_cov_dir = os.path.join(output, "Coverage")
out_gvcfr_dir = os.path.join(output, "GVCF_recal")
out_vcfr_dir = os.path.join(output, "VCF_recal")
out_gvcf_dir = os.path.join(output, "GVCF")
out_vcf_dir = os.path.join(output, "VCF")
out_annot_dir = os.path.join(output, "Annotation")
out_species_dir = os.path.join(output, "Species")
out_uncovered_dir = os.path.join(output, "Uncovered")
out_compare_dir = os.path.join(output, "Compare")
out_table_dir = os.path.join(output, "Table")

highly_hetz_bed = os.path.join(out_vcf_dir, "highly_hetz.bed")
non_genotyped_bed = os.path.join(out_vcf_dir, "non_genotyped.bed")
poorly_covered_bed = os.path.join(out_cov_dir, "poorly_covered.bed")

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    args.sample = sample
    if sample in sample_list_F:

        sample_number = str(sample_list_F.index(sample) + 1)
        sample_total = str(len(sample_list_F))

        out_bqsr_name = sample + ".bqsr.bam"
        output_bqsr_file = os.path.join(out_map_dir, out_bqsr_name)

        if not os.path.isfile(output_bqsr_file):
        
            args.r1_file = r1_file
            args.r2_file = r2_file

            logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

            ##############START PIPELINE#####################
            #################################################

            #INPUT ARGUMENTS
            ################
            check_file_exists(args.r1_file)
            check_file_exists(args.r2_file)

            args.output = os.path.abspath(args.output)
            check_create_dir(args.output)
            #QUALITY CHECK
            ##############
            """
            TODO: Quality check 
            TODO: Human filter
            """
                    
            #QUALITY TRIMMING AND ADAPTER REMOVAL WITH bbduk.sh
            ###################################################
            out_trim_name_r1 = sample + "_R1.clean.fastq.gz"
            out_trim_name_r2 = sample + "_R2.clean.fastq.gz"
            output_trimming_file_r1 = os.path.join(out_trim_dir, out_trim_name_r1)
            output_trimming_file_r2 = os.path.join(out_trim_dir, out_trim_name_r2)
            
            if os.path.isfile(output_trimming_file_r1) and os.path.isfile(output_trimming_file_r2):
                logger.info(YELLOW + DIM + output_trimming_file_r1 + " EXIST\nOmmiting Trimming for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Trimming sample " + sample + END_FORMATTING)
                bbduk_trimming(args)

            #MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
            #####################################################
            out_map_name = sample + ".rg.sorted.bam"
            output_map_file = os.path.join(out_map_dir, out_map_name)

            out_markdup_name = sample + ".rg.markdup.sorted.bam"
            output_markdup_file = os.path.join(out_map_dir, out_markdup_name)

            args.r1_file = output_trimming_file_r1
            args.r2_file = output_trimming_file_r2

            if os.path.isfile(output_map_file) or os.path.isfile(output_markdup_file):
                logger.info(YELLOW + DIM + output_map_file + " EXIST\nOmmiting Mapping for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Mapping sample " + sample + END_FORMATTING)
                logger.info("R1: " + output_trimming_file_r1 + "\nR2: " + output_trimming_file_r2 + "\nReference: " + args.reference)
                bwa_mapping(args)
                sam_to_index_bam(args)

            #MARK DUPLICATES WITH PICARDTOOLS ###################
            #####################################################
            #TO DO: remove output_map_file and include markdup in previous step checking for existence of .rg.markdup.sorted.bam
            out_markdup_name = sample + ".rg.markdup.sorted.bam"
            output_markdup_file = os.path.join(out_map_dir, out_markdup_name)

            args.input_bam = output_map_file

            if os.path.isfile(output_markdup_file):
                logger.info(YELLOW + DIM + output_markdup_file + " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Marking Dupes in sample " + sample + END_FORMATTING)
                logger.info("Input Bam: " + args.input_bam)
                picard_markdup(args)
            
            #CALCULATE COVERAGE FOR EACH POSITION##################
            #######################################################
            out_cov_name = sample + ".cov"
            output_cov_file = os.path.join(out_cov_dir, out_cov_name)

            if os.path.isfile(output_cov_file):
                logger.info(YELLOW + DIM + output_cov_file + " EXIST\nOmmiting coverage calculation for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Calculating coverage in sample " + sample + END_FORMATTING)
                get_coverage(args, output_markdup_file, output_fmt="-d")

            #SPECIES DETERMINATION USING mash #################
            ###################################################
            out_mash_name = sample + ".screen.tab"
            output_mash_file = os.path.join(out_species_dir, out_mash_name)
            
            if os.path.isfile(output_mash_file):
                logger.info(YELLOW + DIM + output_mash_file + " EXIST\nOmmiting Species calculation for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Determining species content in sample " + sample + END_FORMATTING)
                mash_screen(args, winner=True, r2=False, mash_database=args.mash_database)
            

            #HAPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
            #######################################################
            out_gvcfr_name = sample + ".g.vcf"
            output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

            if os.path.isfile(output_gvcfr_file):
                logger.info(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting Haplotype Call (Recall) for sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Haplotype Calling (Recall) in sample " + sample + END_FORMATTING)
                haplotype_caller(args, recalibrate=True, ploidy=2, bamout=False, forceactive=False)
        
        else:
            logger.info(YELLOW + DIM + "\nOMMITING BAM HANDLING FOR SAMPLE " + sample + END_FORMATTING)



#GROUP COVERAGE SUMMARY STATS##########################
#######################################################
group_name = output.split("/")[-1]
logger.info("\n\n" + BLUE + BOLD + "CHECKING LOW COVERED SAMPLES IN GROUP: " + group_name + END_FORMATTING + "\n")

out_cov_name = group_name + ".coverage.tab"
output_cov_file = os.path.join(out_cov_dir, out_cov_name)

if os.path.isfile(output_cov_file):
    logger.info(YELLOW + DIM + output_cov_file + " EXIST\nOmmiting group coverage calculation for group " + group_name + END_FORMATTING)
    saples_low_covered = []
else:
    logger.info(GREEN + "Group coverage stats in group " + group_name + END_FORMATTING)
    saples_low_covered = obtain_group_cov_stats(out_cov_dir, low_cov_threshold=args.mincov, unmmaped_threshold=20)


if os.path.isfile(poorly_covered_bed):
    logger.info(YELLOW + DIM + output_cov_file + " EXIST\nOmmiting poorly covered calculation for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Calculating low covered regions " + group_name + END_FORMATTING)
    poorly_covered_to_bed(out_cov_dir, "poorly_covered", reference="CHROM", min_coverage=2, nocall_fr=0.5)

if len(saples_low_covered) > 0:
    logger.info("\n" + YELLOW + BOLD + "There are sample(s) with low coverage that will be removed from the analysis: " + "\n"\
         + ",".join(saples_low_covered) + END_FORMATTING + "\n")
    remove_low_covered_mixed(args.output, saples_low_covered, "Uncovered")
    #Remove sample from the list of filtered samples
    ################################################
    for samples_to_remove in saples_low_covered:
        sample_list_F.remove(samples_to_remove)
else:
    logger.info("\n" + YELLOW + BOLD + "All samples have a decent depth of coverage according to threshold supplied" + "\n")




#ONCE ALL GVCF VARIANTS ARE CALLED, THEY ARE GATHERED AND FILTERED 
# TO RECALIBRATE ORIGINAL MARKDUPPED BAM
######################################################################
##############START GROUP CALLING FOR RECALIBRATION###################
######################################################################

group_name = output.split("/")[-1]
logger.info("\n\n" + BLUE + BOLD + "STARTING JOINT CALL FOR RECALIBATION IN GROUP: " + group_name + END_FORMATTING + "\n")

#CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_gvcfr_name = group_name + ".cohort.g.vcf"
output_gvcfr_file = os.path.join(out_gvcfr_dir, out_gvcfr_name)

if os.path.isfile(output_gvcfr_file):
    logger.info(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting GVCF Combination (Recall) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "GVCF Combination (Recall) in group " + group_name + END_FORMATTING)
    combine_gvcf(args, recalibrate=True, all_gvcf=False)

#CALL VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_vcfr_name = group_name + ".cohort.raw.vcf"
output_vcfr_file = os.path.join(out_vcfr_dir, out_vcfr_name)

if os.path.isfile(output_vcfr_file):
    logger.info(YELLOW + DIM + output_vcfr_file + " EXIST\nOmmiting Variant Calling (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Variant Calling (Recall-Group) in group " + group_name + END_FORMATTING)
    call_variants(args, recalibrate=True, group=True)

#SELECT VARIANTS 1/2 FOR HARD FILTERING AND RECALIBRATION
#########################################################
out_vcfsnpr_name = group_name + ".cohort.snp.vcf"
out_vcfindelr_name = group_name + ".cohort.indel.vcf"
output_vcfsnpr_file = os.path.join(out_vcfr_dir, out_vcfsnpr_name)
output_vcfindelr_file = os.path.join(out_vcfr_dir, out_vcfindelr_name)

if os.path.isfile(output_vcfsnpr_file) and os.path.isfile(output_vcfindelr_file):
    logger.info(YELLOW + DIM + output_vcfsnpr_file + " EXIST\nOmmiting Variant Selection (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Selecting Variants (Recall-Group) in group " + group_name + END_FORMATTING)
    select_variants(output_vcfr_file, select_type='SNP') #select_variants(output_vcfr_file, select_type='INDEL')
    select_variants(output_vcfr_file, select_type='INDEL')

#HARD FILTER VARIANTS 1/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnpr_name = group_name + ".cohort.snp.hf.vcf"
out_vcfhfindelr_name = group_name + ".cohort.indel.hf.vcf"
output_vcfhfsnpr_file = os.path.join(out_vcfr_dir, out_vcfhfsnpr_name)
output_vcfhfindelr_file = os.path.join(out_vcfr_dir, out_vcfhfindelr_name)


if os.path.isfile(output_vcfhfsnpr_file) and os.path.isfile(output_vcfhfindelr_file):
    logger.info(YELLOW + DIM + output_vcfhfsnpr_file + " EXIST\nOmmiting Hard Filtering (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Hard Filtering Variants (Recall-Group) in group " + group_name + END_FORMATTING)
    hard_filter(output_vcfsnpr_file, select_type='SNP')
    hard_filter(output_vcfindelr_file, select_type='INDEL')

#PASS FILTER VARIANTS 1/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnppass_name = group_name + ".cohort.snp.hf.pass.vcf"
out_vcfhfindelpass_name = group_name + ".cohort.indel.hf.pass.vcf"
output_vcfhfsnppass_file = os.path.join(out_vcfr_dir, out_vcfhfsnppass_name)
output_vcfhfindelpass_file = os.path.join(out_vcfr_dir, out_vcfhfindelpass_name)


if os.path.isfile(output_vcfhfsnpr_file) and os.path.isfile(output_vcfhfsnppass_file):
    logger.info(YELLOW + DIM + output_vcfhfsnppass_file + " EXIST\nOmmiting PASS Filtering (Recall-Group) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "PASS Filtering Variants (Recall-Group) in group " + group_name + END_FORMATTING)
    select_pass_variants(output_vcfhfsnpr_file, nocall_fr=0.2)
    select_pass_variants(output_vcfhfindelr_file, nocall_fr=0.2)


    ######################################################################
    ##############START RECALIBRATION AND FINAL CALL######################
    ######################################################################

logger.info("\n\n" + BLUE + BOLD + "STARTING RECALIBATION IN GROUP: " + group_name + END_FORMATTING)

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)

    args.sample = sample
    args.output = os.path.abspath(args.output)

    if sample in sample_list_F:

        sample_number = str(sample_list_F.index(sample) + 1)
        sample_total = str(len(sample_list_F))

        logger.info("\n" + WHITE_BG + "RECALIBRATION AND CALL ON SAMPLE: " + sample + " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

        ##############START BAM RECALIBRATION############
        #################################################

        ################BQSR AND APPLY BQSR##################
        #####################################################
        out_bqsr_name = sample + ".bqsr.bam"
        output_bqsr_file = os.path.join(out_map_dir, out_bqsr_name)

        if os.path.isfile(output_bqsr_file):
            logger.info(YELLOW + DIM + output_bqsr_file + " EXIST\nOmmiting Recalibration for sample " + sample + END_FORMATTING)
        else:
            logger.info(GREEN + "Recalibration in sample " + sample + END_FORMATTING)
            recalibrate_bam(args)

        #HAPLOTYPE CALL 1/2 FOR HARD FILTERING AND RECALIBRATION
        #######################################################
        out_gvcf_name = sample + ".g.vcf"
        output_gvcf_file = os.path.join(out_gvcf_dir, out_gvcf_name)

        #args.input_bam = output_bqsr_file

        if os.path.isfile(output_gvcf_file):
            logger.info(YELLOW + DIM + output_gvcf_file + " EXIST\nOmmiting Haplotype Call for sample " + sample + END_FORMATTING)
        else:
            logger.info(GREEN + "Haplotype Calling in sample " + sample + END_FORMATTING)
            haplotype_caller(args, recalibrate=False, ploidy=2, bamout=False, forceactive=False)

#ONCE ALL GVCF VARIANTS ARE CALLED, THEY ARE GATHERED AND FILTERED 
# FOR FINAL CALLING
######################################################################
##############START GROUP CALLING FOR FINAL CALL######################
######################################################################
group_name = args.output.split("/")[-1]
logger.info("\n\n" + BLUE + BOLD + "STARTING JOINT CALL FOR FINAL CALL IN GROUP: " + group_name + END_FORMATTING + "\n")

#CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_gvcf_name = group_name + ".cohort.g.vcf"
output_gvcf_file = os.path.join(out_gvcf_dir, out_gvcf_name)

if os.path.isfile(output_gvcf_file):
    logger.info(YELLOW + DIM + output_gvcfr_file + " EXIST\nOmmiting GVCF Combination for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "GVCF Combination in group " + group_name + END_FORMATTING)
    combine_gvcf(args, recalibrate=False, all_gvcf=args.enrich_gvcf)

#CALL VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#######################################################
out_vcf_name = group_name + ".cohort.raw.vcf"
output_vcf_file = os.path.join(out_vcf_dir, out_vcf_name)

if os.path.isfile(output_vcf_file):
    logger.info(YELLOW + DIM + output_vcf_file + " EXIST\nOmmiting Variant Calling (Group) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Variant Calling (Group) in group " + group_name + END_FORMATTING)
    call_variants(args, recalibrate=False, group=True)

#SELECT VARIANTS 2/2 FOR HARD FILTERING AND RECALIBRATION
#########################################################
out_vcfsnp_name = group_name + ".cohort.snp.vcf"
out_vcfindel_name = group_name + ".cohort.indel.vcf"
output_vcfsnp_file = os.path.join(out_vcf_dir, out_vcfsnp_name)
output_vcfindel_file = os.path.join(out_vcf_dir, out_vcfindel_name)

if os.path.isfile(output_vcfsnp_file) and os.path.isfile(output_vcfindel_file):
    logger.info(YELLOW + DIM + output_vcfsnp_file + " EXIST\nOmmiting Variant Selection (Group) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Selecting Variants (Group) in group " + group_name + END_FORMATTING)
    select_variants(output_vcf_file, select_type='SNP')
    select_variants(output_vcf_file, select_type='INDEL')

#HARD FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfsnp_name = group_name + ".cohort.snp.hf.vcf"
out_vcfhfindel_name = group_name + ".cohort.indel.hf.vcf"
output_vcfhfsnp_file = os.path.join(out_vcf_dir, out_vcfhfsnp_name)
output_vcfhfindel_file = os.path.join(out_vcf_dir, out_vcfhfindel_name)


if os.path.isfile(output_vcfhfsnp_file) and os.path.isfile(output_vcfhfindel_file):
    logger.info(YELLOW + DIM + output_vcfhfsnp_file + " EXIST\nOmmiting Hard Filtering (Group) for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Hard Filtering Variants (Group) in group " + group_name + END_FORMATTING)
    hard_filter(output_vcfsnp_file, select_type='SNP')
    hard_filter(output_vcfindel_file, select_type='INDEL')


#PASS FILTER VARIANTS 2/2 FOR RECALIBRATION #############
#########################################################
out_vcfhfcombined_name = group_name + ".cohort.combined.hf.vcf"
output_vcfhfcombined_file = os.path.join(out_vcf_dir, out_vcfhfcombined_name)


if os.path.isfile(output_vcfhfcombined_file):
    logger.info(YELLOW + DIM + output_vcfhfcombined_file + " EXIST\nOmmiting combination for group " + group_name + END_FORMATTING)
else:
    logger.info(GREEN + "Combining both vcf SNP and INDEL in group " + group_name + END_FORMATTING)
    combine_vcf(output_vcfhfsnp_file, output_vcfhfindel_file, name_out=False)

if args.all_cohort == True:
    split_vcf_saples(output_vcfhfcombined_file, sample_list=False, nocall_fr=args.maxnocallfr)
else:
    split_vcf_saples(output_vcfhfcombined_file, sample_list=sample_list_F, nocall_fr=args.maxnocallfr)
###########################################################################
###########################################################################
###########################################################################

logger.info(GREEN + "Determinind highly heterozygous and poorly genotyped regions in " + group_name + END_FORMATTING)
highly_hetz_to_bed(output_vcfhfcombined_file, "highly_hetz", reference="CHROM", nocall_fr=0.5)
non_genotyped_to_bed(output_vcfhfcombined_file, "non_genotyped", reference="CHROM", nocall_fr=0.5)

for r1_file, r2_file in zip(r1, r2):
    sample = extract_sample(r1_file, r2_file)
    args.sample = sample
    args.output = os.path.abspath(args.output)

    if sample in sample_list_F:

        logger.info("\n" + WHITE_BG + "FINAL FILTERING IN SAMPLE " + sample + END_FORMATTING)

        ################FINAL VCF FILTERING##################
        #####################################################
        out_final_name = sample + ".combined.hf.SNP.final.vcf"
        in_final_name = sample + ".combined.hf.vcf"
        output_final_vcf = os.path.join(out_vcf_dir, out_final_name)
        in_final_vcf = os.path.join(out_vcf_dir, in_final_name)

        if os.path.isfile(output_final_vcf):
            logger.info(YELLOW + DIM + output_final_vcf + " EXIST\nOmmiting Final filter for sample " + sample + END_FORMATTING)
        else:
            logger.info(GREEN + "Final filter in sample " + sample + END_FORMATTING)
            vcf_consensus_filter(in_final_vcf, distance=1, AF=0.75, QD=15, window_10=3, dp_limit=8, dp_AF=10, AF_dp=0.80,
             highly_hetz=highly_hetz_bed, 
             non_genotyped=non_genotyped_bed, 
             poorly_covered=poorly_covered_bed, 
             var_type="SNP")

#DETEMINING MIXED ORIGIN IN GROUP######################
#######################################################
output_vcfstat_file = os.path.join(out_table_dir, "vcf_stat.tab")
if os.path.isfile(output_vcfstat_file):
    logger.info("\n" + YELLOW + DIM + output_vcfstat_file + " EXIST\nOmmiting Mixed search in group " + group_name + END_FORMATTING)
    samples_mixed = []
else:
    logger.info(GREEN + "Finding Mixed samples in " + group_name + END_FORMATTING)
    samples_mixed = vcf_stats(out_table_dir, distance=15, quality=10)

if len(samples_mixed) > 0:
    logger.info("\n" + YELLOW + BOLD + "There are mixed sample(s): " + "\n"\
         + ",".join(samples_mixed) + END_FORMATTING + "\n")
    remove_low_covered_mixed(args.output, samples_mixed, "Mixed")
    #Remove sample from the list of filtered samples
    ################################################
    for samples_to_remove in samples_mixed:
        sample_list_F.remove(samples_to_remove)
else:
    logger.info("\n" + YELLOW + BOLD + "No mixed samples have been detected" + "\n")

logger.info("\n\n" + MAGENTA + BOLD + "VARIANT CALL FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")

#######################################################################################################################################
#################################END OF VARIANT CALLING################################################################################
#######################################################################################################################################
tuberculosis = False
if tuberculosis == True:
    logger.info("\n\n" + BLUE + BOLD + "STARTING ANNOTATION IN GROUP: " + group_name + END_FORMATTING + "\n")

    for root, _, files in os.walk(out_vcf_dir):
        for name in files:
            filename = os.path.join(root, name)
            output_path = os.path.join(out_annot_dir, name)
            if filename.endswith("combined.hf.vcf"):
                sample = name.split(".")[0]
                if sample in sample_list_F:
                    #ANNOTATION -AUTO AND SPECIFIC- ###################
                    ###################################################
                    out_annot_name = sample + ".combined.hf.annot.tsv"
                    output_annot_file = os.path.join(out_annot_dir, out_annot_name)

                    if os.path.isfile(output_annot_file):
                        logger.info(YELLOW + DIM + output_annot_file + " EXIST\nOmmiting Annotation for sample " + sample + END_FORMATTING)
                    else:
                        logger.info(GREEN + "Annotating snps in sample " + sample + END_FORMATTING)
                        replace_reference(filename, output_path)
                        snpeff_annotation(args, output_path, database=args.snpeff_database)
                        #Handle output vcf file from SnpEff annotation
                        vcf_path = (".").join(output_path.split(".")[:-1])
                        annot_vcf = vcf_path + ".annot"
                        #This function add SPECIFIC anotation
                        if args.annot_bed:
                            final_annotation(annot_vcf, *args.annot_bed)
                        else:
                            final_annotation(annot_vcf)


    logger.info("\n\n" + MAGENTA + BOLD + "ANNOTATION FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")
else:
    logger.info("NO TB Selected, snpEff won't be executed")




logger.info("\n\n" + BLUE + BOLD + "STARTING COMPARISON IN GROUP: " + group_name + END_FORMATTING + "\n")

check_create_dir(out_compare_dir)
folder_compare = today + "_" + group_name
path_compare = os.path.join(out_compare_dir, folder_compare)
check_create_dir(path_compare)
full_path_compare = os.path.join(path_compare, group_name)

compare_snp_matrix = full_path_compare + ".revised.tsv"

ddtb_add(out_vcf_dir, full_path_compare, recalibrate=args.output)
ddtb_compare(compare_snp_matrix)

logger.info("\n\n" + MAGENTA + BOLD + "COMPARING FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")


if args.noclean == True:
    logger.info("\n\n" + BLUE + BOLD + "STARTING CLEANING IN GROUP: " + group_name + END_FORMATTING + "\n")
    clean_unwanted_files(args)
else:
    logger.info("No cleaning was requested")

logger.info("\n\n" + MAGENTA + BOLD + "#####END OF PIPELINE SNPTB#####" + END_FORMATTING + "\n")

#./snptb_runner.py -i /home/laura/ANALYSIS/Lofreq/coinfection_designed/raw -r reference/MTB_ancestorII_reference.fasta -o /home/laura/ANALYSIS/Lofreq/coinfection_designed/TEST -s sample_list.txt
#/home/laura/DEVELOP/SNPTB/snptb_runner.py -i /home/laura/RAW/Mixtas_Italia/ -r /home/laura/DATABASES/REFERENCES/ancestorII/MTB_ancestorII_reference.fasta -o /home/laura/ANALYSIS/Lofreq/coinfection_italy/
