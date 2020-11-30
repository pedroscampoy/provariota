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
from misc import check_file_exists, extract_sample, check_create_dir, execute_subprocess, \
    extract_read_list, file_to_list, obtain_group_cov_stats, clean_unwanted_files, \
    check_reanalysis, vcf_stats, remove_low_quality, obtain_overal_stats
from preprocessing import fastqc_quality, fastp_trimming, format_html_image
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_variant import picard_dictionary, samtools_faidx, picard_markdup, ivar_consensus, \
    replace_consensus_header, create_bamstat, create_coverage, freebayes_variants
from vcf_process import filter_tsv_variants, vcf_to_ivar_tsv
from annotation import annotate_snpeff, annotate_pangolin, user_annotation
from compare_snp import ddtb_add, ddtb_compare, ddbb_create_intermediate, revised_df

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 22 Sep 2020
REVISION: 


TODO:
    Adapt check_reanalysis
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

logger = logging.getLogger()

def main():
    """
    Create main function to capture code errors: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python
    """

    #ARGUMENTS

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'provaryota.py', description= 'Pipeline to call variants (SNVs) with any non model organism. Specialised in prokaryotic organisms')
        
        input_group = parser.add_argument_group('Input', 'Input parameters')

        input_group.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all fast[aq] files')
        input_group.add_argument('-r', '--reference', metavar="reference", type=str, required=True, help='REQUIRED. File to map against')
        input_group.add_argument('-s', '--sample', metavar="sample", type=str, required=False, help='Sample to identify further files')
        input_group.add_argument('-S', '--sample_list', type=str, required=False, help='Sample names to analyse only in the file supplied')
        input_group.add_argument('-p', '--primers', type=str, default=False, required=False, help='Bed file including primers to trim')
        
        quality_group = parser.add_argument_group('Quality parameters', 'parameters for diferent triming conditions')

        quality_group.add_argument('-c', '--mincoverage', type=int, default=10, required=False, help='Minimum percentage of uncovered (Default 10)')
        quality_group.add_argument('-n', '--min_snp', type=int, required=False, default=3, help='SNP number to pass quality threshold')

        output_group = parser.add_argument_group('Output', 'Required parameter to output results')

        output_group.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')
        output_group.add_argument('-C', '--noclean', required=False, action='store_false', help='Clean unwanted files for standard execution')

        params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

        params_group.add_argument('-T', '--threads', type=str, dest="threads", required=False, default=16, help='Threads to use')
        params_group.add_argument('-M', '--memory', type=str, dest="memory", required=False, default=32, help='Max memory to use')

        annot_group = parser.add_argument_group('Annotation', 'parameters for variant annotation')

        annot_group.add_argument('-B', '--annot_bed', type=str, default=[], required=False, action='append', help='bed file to annotate')
        annot_group.add_argument('-V', '--annot_vcf', type=str, default=[], required=False, action='append', help='vcf file to annotate')
        
        annot_group = parser.add_argument_group('Annotation', 'parameters for variant annotation')

        annot_group.add_argument('--mash_database', type=str, required=False, default=False, help='MASH ncbi annotation containing all species database')
        annot_group.add_argument('--snpeff_database', type=str, required=False, default=False, help='snpEFF annotation database')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()


    ######################################################################
    #####################START PIPELINE###################################
    ######################################################################
    output = os.path.abspath(args.output)
    group_name = output.split("/")[-1]
    reference = os.path.abspath(args.reference)

    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = group_name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output, 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

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

    #Check if there are samples to filter out
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
    
    #picard_dictionary(args)
    samtools_faidx(args.reference)

    #DECLARE FOLDERS CREATED IN PIPELINE ################
    #AND KEY FILES ######################################
    #####################################################
    #Annotation related parameters
    #script_dir = os.path.dirname(os.path.realpath(__file__))

    #Output related
    #out_ref_dir = os.path.join(output, "Reference")
    out_qc_dir = os.path.join(output, "Quality")
    out_qc_pre_dir = os.path.join(out_qc_dir, "raw") #subfolder
    out_qc_post_dir = os.path.join(out_qc_dir, "processed") #subfolder
    out_trim_dir = os.path.join(output, "Trimmed")
    out_map_dir = os.path.join(output, "Bam")
    out_variant_dir = os.path.join(output, "Variants")
    out_variant_freebayes_dir = os.path.join(out_variant_dir, "freebayes_raw") #subfolder
    out_filtered_freebayes_dir = os.path.join(out_variant_dir, "freebayes_filtered") #subfolder
    out_consensus_dir = os.path.join(output, "Consensus")

    out_stats_dir = os.path.join(output, "Stats")
    out_stats_bamstats_dir = os.path.join(out_stats_dir, "Bamstats") #subfolder
    out_stats_coverage_dir = os.path.join(out_stats_dir, "Coverage") #subfolder
    out_compare_dir = os.path.join(output, "Compare")

    out_annot_dir = os.path.join(output, "Annotation")
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff") #subfolder
    out_annot_user_dir = os.path.join(out_annot_dir, "user") #subfolder


    for r1_file, r2_file in zip(r1, r2):
        #EXtract sample name
        sample = extract_sample(r1_file, r2_file)
        args.sample = sample
        if sample in sample_list_F:

            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            out_markdup_name = sample + ".rg.markdup.sorted.bam"
            output_markdup_file = os.path.join(out_map_dir, out_markdup_name)

            logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

            if not os.path.isfile(output_markdup_file):
            
                args.r1_file = r1_file
                args.r2_file = r2_file

                ##############START PIPELINE#####################
                #################################################

                #INPUT ARGUMENTS
                ################
                check_file_exists(r1_file)
                check_file_exists(r2_file)

                args.output = os.path.abspath(args.output)
                check_create_dir(args.output)

                ######################QUALITY CHECK in RAW with fastqc
                ######################################################
                check_create_dir(out_qc_dir)

                out_qc_raw_name_r1 = (".").join(r1_file.split('/')[-1].split('.')[0:-2]) + '_fastqc.html'
                out_qc_raw_name_r2 = (".").join(r2_file.split('/')[-1].split('.')[0:-2]) + '_fastqc.html'
                output_qc_raw_file_r1 = os.path.join(out_qc_pre_dir, out_qc_raw_name_r1)
                output_qc_raw_file_r2 = os.path.join(out_qc_pre_dir, out_qc_raw_name_r2)
                                
                if os.path.isfile(output_qc_raw_file_r1) and os.path.isfile(output_qc_raw_file_r2):
                    logger.info(YELLOW + DIM + output_qc_raw_file_r1 + " EXIST\nOmmiting QC for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Checking quality in sample " + sample + END_FORMATTING)
                    logger.info("R1: " + r1_file + "\nR2: " + r2_file)
                    fastqc_quality(r1_file, r2_file, out_qc_pre_dir, args.threads)

                """
                TODO: Human filter
                """
                
                ####QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
                ###################################################
                out_trim_name_r1 = sample + ".trimmed_R1.fastq.gz"
                out_trim_name_r2 = sample + ".trimmed_R2.fastq.gz"
                output_trimming_file_r1 = os.path.join(out_trim_dir, out_trim_name_r1)
                output_trimming_file_r2 = os.path.join(out_trim_dir, out_trim_name_r2)
                
                if os.path.isfile(output_trimming_file_r1) and os.path.isfile(output_trimming_file_r2):
                    logger.info(YELLOW + DIM + output_trimming_file_r1 + " EXIST\nOmmiting Trimming for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Trimming sample " + sample + END_FORMATTING)
                    fastp_trimming(r1_file, r2_file, sample, out_trim_dir, threads=args.threads, min_qual=20, window_size=10, min_len=35)


                ##################QUALITY CHECK in TRIMMED with fastqc
                ######################################################
                check_create_dir(out_qc_dir)
                
                out_qc_pos_r1 = sample + ".trimmed_R1_fastqc.html"
                out_qc_pos_r2 = sample + ".trimmed_R2_fastqc.html"
                output_qc_precessed_file_r1 = os.path.join(out_qc_post_dir, out_qc_pos_r1)
                output_qc_precessed_file_r2 = os.path.join(out_qc_post_dir, out_qc_pos_r2)
                                
                if os.path.isfile(output_qc_precessed_file_r1) and os.path.isfile(output_qc_precessed_file_r2):
                    logger.info(YELLOW + DIM + output_qc_raw_file_r1 + " EXIST\nOmmiting QC for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Checking quality in processed sample " + sample + END_FORMATTING)
                    logger.info("R1: " + output_trimming_file_r1 + "\nR2: " + output_trimming_file_r2)
                    fastqc_quality(output_trimming_file_r1, output_trimming_file_r2, out_qc_post_dir, args.threads)

                #MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
                #####################################################
                out_map_name = sample + ".rg.sorted.bam"
                output_map_file = os.path.join(out_map_dir, out_map_name)

                if os.path.isfile(output_map_file):
                    logger.info(YELLOW + DIM + output_map_file + " EXIST\nOmmiting Mapping for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Mapping sample " + sample + END_FORMATTING)
                    logger.info("R1: " + output_trimming_file_r1 + "\nR2: " + output_trimming_file_r2 + "\nReference: " + reference)
                    bwa_mapping(output_trimming_file_r1, output_trimming_file_r2, reference, sample, out_map_dir, threads=args.threads)
                    sam_to_index_bam(sample, out_map_dir, output_trimming_file_r1, threads=args.threads)


                #MARK DUPLICATES WITH PICARDTOOLS ###################
                #####################################################
                

                if os.path.isfile(output_markdup_file):
                    logger.info(YELLOW + DIM + output_markdup_file + " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Marking Dupes in sample " + sample + END_FORMATTING)
                    logger.info("Input Bam: " + output_map_file)
                    picard_markdup(output_map_file)
                
            else:
                logger.info(YELLOW + DIM + output_markdup_file + " EXIST\nOmmiting BAM mapping and BAM manipulation in sample " + sample + END_FORMATTING)
            
            ########################END OF MAPPING AND BAM MANIPULATION#####################################################################
            ################################################################################################################################
            
            #VARIANT CALLING WTIH ivar variants##################
            #####################################################
            check_create_dir(out_variant_dir)
            check_create_dir(out_variant_freebayes_dir)
            out_freebayes_variant_name = sample + ".vcf"
            out_freebayes_variant_file = os.path.join(out_variant_freebayes_dir, out_freebayes_variant_name)

            if os.path.isfile(out_freebayes_variant_file):
                logger.info(YELLOW + DIM + out_freebayes_variant_file + " EXIST\nOmmiting Variant call for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Calling variants with freebayes in sample " + sample + END_FORMATTING)
                freebayes_variants(reference, output_markdup_file, out_variant_freebayes_dir, sample)


            #VARIANT FORMAT ADAPTATION TO IVAR ##################
            #####################################################
            check_create_dir(out_filtered_freebayes_dir)
            out_ivar_variant_name = sample + ".tsv"
            out_freebayes_filtered_file = os.path.join(out_filtered_freebayes_dir, out_ivar_variant_name)
            

            if os.path.isfile(out_freebayes_filtered_file):
                logger.info(YELLOW + DIM + out_freebayes_filtered_file + " EXIST\nOmmiting format adaptation for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Adapting variants format in sample " + sample + END_FORMATTING)
                vcf_to_ivar_tsv(out_freebayes_variant_file, out_freebayes_filtered_file)
            
            #CREATE CONSENSUS with freebayes consensus##################
            #######################################################
            check_create_dir(out_consensus_dir)
            out_freebayes_consensus_name = sample + ".fa"
            out_freebayes_consensus_file = os.path.join(out_consensus_dir, out_freebayes_consensus_name)

            if os.path.isfile(out_freebayes_consensus_file):
                logger.info(YELLOW + DIM + out_freebayes_consensus_file + " EXIST\nOmmiting Consensus for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating consensus with freebayes in sample " + sample + END_FORMATTING)
                ivar_consensus(output_markdup_file, out_consensus_dir, sample, min_quality=20, min_frequency_threshold=0.8, min_depth=20, uncovered_character='N')
                logger.info(GREEN + "Replacing consensus header in " + sample + END_FORMATTING)
                replace_consensus_header(out_freebayes_consensus_file)


            ########################CREATE STATS AND QUALITY FILTERS########################################################################
            ################################################################################################################################
            #CREATE Bamstats#######################################
            #######################################################
            check_create_dir(out_stats_dir)
            check_create_dir(out_stats_bamstats_dir)
            out_bamstats_name = sample + ".bamstats"
            out_bamstats_file = os.path.join(out_stats_bamstats_dir, out_bamstats_name)

            if os.path.isfile(out_bamstats_file):
                logger.info(YELLOW + DIM + out_bamstats_file + " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating bamstats in sample " + sample + END_FORMATTING)
                create_bamstat(output_markdup_file, out_stats_bamstats_dir, sample, threads=args.threads)

            #CREATE Bamstats#######################################
            #######################################################
            check_create_dir(out_stats_coverage_dir)
            out_coverage_name = sample + ".cov"
            out_coverage_file = os.path.join(out_stats_coverage_dir, out_coverage_name)

            if os.path.isfile(out_coverage_file):
                logger.info(YELLOW + DIM + out_coverage_file + " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating coverage in sample " + sample + END_FORMATTING)
                create_coverage(output_markdup_file, out_stats_coverage_dir, sample)

            
    
    ###################fastqc OUTPUT FORMAT FOR COMPARISON
    ######################################################
    logger.info(GREEN + "Creating summary report for quality result " + END_FORMATTING)
    format_html_image(out_qc_dir)

    ###############################coverage OUTPUT SUMMARY
    ######################################################
    logger.info(GREEN + "Creating summary report for coverage result " + END_FORMATTING)
    #obtain_group_cov_stats(out_stats_coverage_dir, group_name)

    #####################READS and VARIANTS OUTPUT SUMMARY
    ######################################################
    logger.info(GREEN + "Creating overal summary report " + END_FORMATTING)
    #obtain_overal_stats(output, group_name)

    ######################################REMOVE UNCOVERED
    ##############################################################################################################################
    logger.info(GREEN + "Removing low quality samples" + END_FORMATTING)
    #remove_low_quality(output, min_coverage=args.mincoverage, min_hq_snp=args.min_snp, type_remove='Uncovered')

    #ANNOTATION WITH SNPEFF AND USER INPUT ##############
    #####################################################
    logger.info("\n\n" + BLUE + BOLD + "STARTING ANNOTATION IN GROUP: " + group_name + END_FORMATTING + "\n")
    check_create_dir(out_annot_dir)
    check_create_dir(out_annot_snpeff_dir)
    ####SNPEFF
    if args.snpeff_database != False:
        for root, _, files in os.walk(out_filtered_freebayes_dir):
            if root == out_filtered_freebayes_dir: 
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(out_annot_snpeff_dir, sample + ".annot")
                        if os.path.isfile(out_annot_file):
                            logger.info(YELLOW + DIM + out_annot_file + " EXIST\nOmmiting snpEff Annotation for sample " + sample + END_FORMATTING)
                        else:
                            logger.info(GREEN + "Annotating sample with snpEff: " + sample + END_FORMATTING)
                            output_vcf = os.path.join(out_annot_snpeff_dir, sample + '.vcf')
                            annotate_snpeff(filename, output_vcf, out_annot_file, database=args.snpeff_database)
    ####USER DEFINED
    if not args.annot_bed and not args.annot_vcf:
        logger.info(YELLOW + BOLD + "Ommiting User Annotation, no BED or VCF files supplied" + END_FORMATTING)
    else:
        check_create_dir(out_annot_user_dir)
        for root, _, files in os.walk(out_filtered_freebayes_dir):
            if root == out_filtered_freebayes_dir: 
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(out_annot_user_dir, sample + ".tsv")
                        user_annotation(filename, out_annot_file, vcf_files=args.annot_vcf, bed_files=args.annot_bed)


    ################SNP COMPARISON using tsv variant files
    ######################################################
    logger.info("\n\n" + BLUE + BOLD + "STARTING COMPARISON IN GROUP: " + group_name + END_FORMATTING + "\n")

    check_create_dir(out_compare_dir)
    folder_compare = today + "_" + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)

    #ddtb_add(out_filtered_freebayes_dir, full_path_compare)
    compare_snp_matrix_recal = full_path_compare + ".revised.final.tsv"
    compare_snp_matrix_recal_intermediate = full_path_compare + ".revised_intermediate.tsv"
    recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(out_filtered_freebayes_dir, out_stats_coverage_dir, min_freq_discard=0.1, min_alt_dp=4)
    recalibrated_snp_matrix_intermediate.to_csv(compare_snp_matrix_recal_intermediate, sep="\t", index=False)
    recalibrated_revised_df = revised_df(recalibrated_snp_matrix_intermediate, path_compare, min_freq_include=0.7, min_threshold_discard_uncov_sample=0.4, min_threshold_discard_uncov_pos=0.4, min_threshold_discard_htz_sample=0.7, min_threshold_discard_htz_pos=0.4, remove_faulty=True, drop_samples=True, drop_positions=True)
    recalibrated_revised_df.to_csv(compare_snp_matrix_recal, sep="\t", index=False)
    ddtb_compare(compare_snp_matrix_recal, distance=5)


    logger.info("\n\n" + MAGENTA + BOLD + "COMPARING FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")

    """
 
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

    

    #ddtb_add(out_vcf_dir, full_path_compare)
    ddtb_add(out_vcf_dir, full_path_compare, recalibrate=args.output)

    compare_snp_matrix = full_path_compare + ".revised.tsv"
    
    ddtb_compare(compare_snp_matrix)

    logger.info("\n\n" + MAGENTA + BOLD + "COMPARING FINISHED IN GROUP: " + group_name + END_FORMATTING + "\n")


    if args.noclean == True:
        logger.info("\n\n" + BLUE + BOLD + "STARTING CLEANING IN GROUP: " + group_name + END_FORMATTING + "\n")
        clean_unwanted_files(args)
    else:
        logger.info("No cleaning was requested")

    
    """
    logger.info("\n\n" + MAGENTA + BOLD + "#####END OF PIPELINE SNPTB#####" + END_FORMATTING + "\n")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise