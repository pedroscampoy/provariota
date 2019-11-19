#!/usr/bin/env python

import os
import argparse
import subprocess
import shutil
from misc import check_file_exists, obtain_output_dir, check_create_dir, get_picard_path, execute_subprocess, check_remove_file, \
    longest_common_suffix



"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 18 April 2019
REVISION:

TODO
    -Handle parameters
    -Explore samtools markdup

================================================================
END_OF_HEADER
================================================================
"""


#COLORS AND AND FORMATTING
"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
"""
END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
DIM = '\033[2m'

def get_arguments():

    parser = argparse.ArgumentParser(prog = 'bam_recalibrate.py', description= 'Call variants, get HQ SNPs and use them to realibrate')
    
    input_group = parser.add_argument_group('Input', 'Required files to map')

    input_group.add_argument('-b', '--bam', dest="input_bam", metavar="input.bam", type=str, required=True, help='Input bam to recalibrate and call')
    input_group.add_argument('-r', '--reference', metavar="ref.fasta", type=str, required=True, help='File to map against')

    output_group = parser.add_argument_group('Output')

    output_group.add_argument('-o', '--output', type=str, required=False, help='Output directory')
    output_group.add_argument('-s', '--sample', type=str, required=False, help='Sample name to handle output files ')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    
    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=1, help='Threads to use')


    #argcomplete.autocomplete(parser)
    arguments = parser.parse_args()

    return arguments

#args = get_arguments()

def samtools_markdup(args):
    #http://www.htslib.org/doc/samtools.html
    # Add ms and MC tags for markdup to use later
    #samtools fixmate -m namesort.bam fixmate.bam
    # Markdup needs position order
    #samtools sort -o positionsort.bam fixmate.bam
    # Finally mark duplicates
    #samtools markdup positionsort.bam markdup.bam
    pass

def picard_markdup(args):
    #java -jar picard.jar MarkDuplicates \
    #  I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
    picard_jar = get_picard_path()
    
    input_bam = os.path.abspath(args.input_bam)
    in_param = "I=" + input_bam
    
    path_file_name = input_bam.split(".")[0]
    file_name = path_file_name.split("/")[-1]
    output_markdup = path_file_name + ".rg.markdup.bam"
    output_markdup_sorted = path_file_name + ".rg.markdup.sorted.bam"
    out_param = "O=" + output_markdup

    stat_output_dir = obtain_output_dir(args, "Stats")
    stat_output_file = file_name + ".markdup.metrics.txt"
    stat_output_full = os.path.join(stat_output_dir, stat_output_file)
    stats_param = "M=" + stat_output_full

    check_create_dir(stat_output_dir)

    cmd_markdup = ["java", "-jar", picard_jar, "MarkDuplicates", 
    in_param, out_param, stats_param]
    execute_subprocess(cmd_markdup)
    

    #samtools sort: samtools sort $output_dir/$sample".sorted.bam" -o $output_dir/$sample".sorted.bam"
    cmd_sort = ["samtools", "sort", output_markdup, "-o", output_markdup_sorted]
    execute_subprocess(cmd_sort)

    #Handled in Haplotype Caller function
    #samtools index: samtools index $output_dir/$sample".sorted.bam"
    subprocess.run(["samtools", "index", output_markdup_sorted], 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    check_remove_file(input_bam)
    check_remove_file(output_markdup)

def picard_dictionary(args):
    #java -jar picard.jar CreateSequenceDictionary\
    # R=reference.fasta O=reference.dict
    picard_jar = get_picard_path()

    input_reference = os.path.abspath(args.reference)
    ref_param = "R=" + input_reference

    path_file_list = input_reference.split(".")[:-1]
    path_file_name = ".".join(path_file_list)
    dict_file_name = path_file_name + ".dict"
    out_param = "O=" + dict_file_name

    if os.path.exists(dict_file_name):
        print(dict_file_name + " already EXIST")
    else:
        cmd = ["java", "-jar", picard_jar, "CreateSequenceDictionary", 
        ref_param, out_param]
        execute_subprocess(cmd)


def samtools_faidx(args):
    #samtools faidx reference.fa

    input_reference = os.path.abspath(args.reference)
    fai_file_name = input_reference + ".fai"

    if os.path.exists(fai_file_name):
        print(fai_file_name + " already EXIST")
    else:
        cmd = ["samtools", "faidx", input_reference]
        execute_subprocess(cmd)


def haplotype_caller(args, recalibrate=False, ploidy=2, bamout=False, forceactive=False, intervals=False):
    #base_quality=13, 
    """
    #No excuses
    https://software.broadinstitute.org/gatk/documentation/article?id=11081
    """
    #input_bam = os.path.abspath(args.input_bam)
    input_reference = os.path.abspath(args.reference)
    
    bam_output_dir = obtain_output_dir(args, "Bam")
    #file_name = path_file_name.split("/")[-1] #sample_name
    file_name = args.sample
    #path_file_name = os.path.join(output_dir, gvcf_output_file)

    if recalibrate:
        input_bam_to_call_name = file_name + ".rg.markdup.sorted.bam" 

        gvcf_output_dir = obtain_output_dir(args, "GVCF_recal")
        gvcf_output_file = file_name + ".g.vcf"
    else:
        input_bam_to_call_name = file_name + ".bqsr.bam"

        gvcf_output_dir = obtain_output_dir(args, "GVCF")
        gvcf_output_file = file_name + ".g.vcf"

    check_create_dir(gvcf_output_dir)

    input_bam_to_call = os.path.join(bam_output_dir, input_bam_to_call_name)
    gvcf_output_full = os.path.join(gvcf_output_dir, gvcf_output_file)

    memory_param = "-Xmx" + str(args.memory) + "g"

    hc_args = ["gatk", "HaplotypeCaller",
    "--java-options", memory_param,
    "--reference", input_reference,
    "--input", input_bam_to_call,
    "--output", gvcf_output_full,
    "--emit-ref-confidence", "GVCF",
    "--annotation-group", "AS_StandardAnnotation",
    "--sample-ploidy", str(ploidy)
    ]

#"--min-base-quality-score", str(base_quality),

    #Create bam index
    #cmd_index = ["samtools", "index", input_bam_to_call]
    #execute_subprocess(cmd_index)

    if bamout:
        bamout_output_dir = obtain_output_dir(args, "Bamout")
        bamout_output_file = file_name + ".p" + str(ploidy) + ".out.bam"
        bamout_output_full = os.path.join(bamout_output_dir, bamout_output_file)
        check_create_dir(bamout_output_dir)
        bamout_params = ["--bam-output", bamout_output_full]
        hc_args.extend(bamout_params)
    
    if forceactive:
        force_params = ["--force-active", "--disable-optimizations"]
        hc_args.extend(force_params)

    execute_subprocess(hc_args)

    """
    gatk --java-options "-Xmx4g" HaplotypeCaller --reference ref.fasta --input markdup.bam --output sample.g.vcf" \
   --intervals MTB_anc_PAIR1_50_50:1*-100000 --emit-ref-confidence GVCF --bam-output "Bamout/"$i".out.bam" \
   --force-active --disable-optimizations --dont-trim-active-regions --max-num-haplotypes-in-population 2 \
   --assembly-region-padding 5000 --max-assembly-region-size 10000 --assembly-region-out "/home/laura/ANALYSIS/GATK/coinfection_PAIR15050/Assembly/"$i"_assembly.tsv" \
   --annotation-group AS_StandardAnnotation --annotation-group StandardHCAnnotation
    """

def call_variants(args, recalibrate=False, group=True):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
    #Call variants:
    gatk --java-options "-Xmx4g" GenotypeGVCFs -R Homo_sapiens_assembly38.fasta -V input.g.vcf.gz -O output.vcf.gz
    """
    output = os.path.abspath(args.output)

    input_reference = os.path.abspath(args.reference)

    if not args.sample:
        args.sample = "nosample"

    file_name = args.sample #sample_name
    group_name = output.split("/")[-1] #group_name

    if recalibrate:
        
        gvcf_input_dir = obtain_output_dir(args, "GVCF_recal")
        vcf_output_dir = obtain_output_dir(args, "VCF_recal")
    else:
        gvcf_input_dir = obtain_output_dir(args, "GVCF")
        vcf_output_dir = obtain_output_dir(args, "VCF")

    if group:
        gvcf_input_file = group_name + ".cohort.g.vcf"
        vcf_output_file = group_name + ".cohort.raw.vcf"
    else:
        gvcf_input_file = file_name + ".g.vcf"
        vcf_output_file = file_name + ".raw.vcf"

    gvcf_input_full = os.path.join(gvcf_input_dir, gvcf_input_file)
    vcf_output_full = os.path.join(vcf_output_dir, vcf_output_file)

    check_create_dir(gvcf_input_dir)
    check_create_dir(vcf_output_dir)

    memory_param = "-Xmx" + str(args.memory) + "g"

    cmd = ["gatk", "GenotypeGVCFs",
    "--java-options", memory_param,
    "--reference", input_reference,
    "--variant", gvcf_input_full,
    "--output", vcf_output_full]
    
    execute_subprocess(cmd)
    

def select_variants(raw_vcf, select_type='SNP'):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php
    gatk SelectVariants -V cohort.vcf.gz -select-type SNP -O snps.vcf.gz
    """
    if select_type == "SNP":
        extension = ".snp.vcf"
    elif select_type == "INDEL":
        extension = ".indel.vcf"
    else:
        print(RED + BOLD + "Choose a correct type to filter" + END_FORMATTING)

    input_vcf = os.path.abspath(raw_vcf)
    check_file_exists(input_vcf)
    
    raw_vcf_file_name = (".").join(input_vcf.split(".")[:-2])
    #file_name = raw_vcf_file_name.split("/")[-1] #sample_name

    vcf_selected_output_file = raw_vcf_file_name + extension

    #memory_param = "-Xmx" + str(args.memory) + "g"
    #"--java-options", memory_param,

    cmd = ["gatk", "SelectVariants", 
    "--variant", input_vcf,
    "--select-type-to-include", select_type,
    "--select-type-to-include", "MIXED",
    "--output", vcf_selected_output_file]

#    "--remove-unused-alternates",



    execute_subprocess(cmd)

def hard_filter(selected_vcf, select_type='SNP'):
    """
    https://software.broadinstitute.org/gatk/documentation/article.php?id=6925
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.php
    https://software.broadinstitute.org/gatk/documentation/article?id=23216
    SNP:
    gatk VariantFiltration -V snps.vcf.gz "--filter-expression", "QD < 2.0", "--filter-name", "QD2" \
    "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30" "--filter-expression", "SOR > 3.0", "--filter-name", "SOR3" "--filter-expression", "FS > 60.0", "--filter-name", "FS60" \
    "--filter-expression", "MQ < 40.0", "--filter-name", "MQ40" "--filter-expression", "MQRankSum < -12.5", "--filter-name", "MQRankSum-12.5" "--filter-expression", "ReadPosRankSum < -8.0" \
   , "--filter-name", "ReadPosRankSum-8" -O snps_filtered.vcf.gz
    INDEL:
    gatk VariantFiltration -V indels.vcf.gz "--filter-expression", "QD < 2.0", "--filter-name", "QD2" "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30" \
    -"--filter-expression", "FS > 200.0", "--filter-name", "FS200" -"--filter-expression", "ReadPosRankSum < -20.0", "--filter-name", "ReadPosRankSum-20" -O indels_filtered.vcf.gz
    #--filterExpression "QD<2.0||FS>60.0||MQ<40.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filterName "my_snp_filter" 
    """

    input_vcf = os.path.abspath(selected_vcf)
    check_file_exists(input_vcf)
    
    selected_vcf_file_name = (".").join(input_vcf.split(".")[:-2])

    if select_type == "SNP":
        extension = ".snp.hf.vcf"
        vcf_hard_filtered_output_file = selected_vcf_file_name + extension
        cmd = ["gatk", "VariantFiltration", 
            "--variant", input_vcf,
            "--filter-expression", "QD < 2.0", "--filter-name", "QD2",
            "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30",
            "--filter-expression", "SOR > 3.5", "--filter-name", "SOR3",
            "--filter-expression", "FS > 60.0", "--filter-name", "FS60",
            "--filter-expression", "MQ < 40.0", "--filter-name", "MQ40",
            "--filter-expression", "DP < 10", "--filter-name", "DP10",
            "--filter-expression", "MQRankSum < -12.5", "--filter-name", "MQRankSum-12.5",
            "--filter-expression", "ReadPosRankSum < -8.0", "--filter-name", "ReadPosRankSum-8",
            "--output", vcf_hard_filtered_output_file]

    elif select_type == "INDEL":
        extension = ".indel.hf.vcf"
        vcf_hard_filtered_output_file = selected_vcf_file_name + extension
        cmd = ["gatk", "VariantFiltration", 
            "--variant", input_vcf,
            "--filter-expression", "QD < 2.0", "--filter-name", "QD2",
            "--filter-expression", "QUAL < 30.0", "--filter-name", "QUAL30",
            "--filter-expression", "SOR > 10.0", "--filter-name", "SOR10",
            "--filter-expression", "FS > 200.0", "--filter-name", "FS200",
            "--filter-expression", "ReadPosRankSum < -20.0", "--filter-name", "ReadPosRankSum-20",
            "--output", vcf_hard_filtered_output_file]
    else:
        print(RED + BOLD + "Choose a correct type to filter" + END_FORMATTING)

    execute_subprocess(cmd)


def combine_gvcf(args, recalibrate=False, all_gvcf=False):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_CombineGVCFs.php
    #combined multi-sample gVCF:
    gatk CombineGVCFs -R reference.fasta --variant sample1.g.vcf.gz --variant sample2.g.vcf.gz -O cohort.g.vcf.gz
    """
    output = os.path.abspath(args.output)
    input_reference = os.path.abspath(args.reference)
    
    group_name = output.split("/")[-1] #group_name

    if recalibrate:
        gvcf_input_dir = obtain_output_dir(args, "GVCF_recal")
    else:
        gvcf_input_dir = obtain_output_dir(args, "GVCF")

    gvcf_output_file = group_name + ".cohort.g.vcf"
    gvcf_output_full = os.path.join(gvcf_input_dir, gvcf_output_file)

    check_create_dir(gvcf_input_dir)

    memory_param = "-Xmx" + str(args.memory) + "g"

    cmd = ["gatk", "CombineGVCFs",
    "--java-options", memory_param,
    "--reference", input_reference,
    "--output", gvcf_output_full]

    for root, _, files in os.walk(gvcf_input_dir):
        for name in files:
            filename = os.path.join(root, name)
            if filename.endswith(".g.vcf"):
                cmd.append("--variant")
                cmd.append(filename)
    if all_gvcf != False:
        if os.path.isdir(all_gvcf):
            all_gvcf = os.path.abspath(all_gvcf)
            print("Using gvcf from enricment folder:" + all_gvcf)
            for root, _, files in os.walk(all_gvcf):
                for name in files:
                    filename = os.path.join(root, name)
                    if filename.endswith(".g.vcf"):
                        cmd.append("--variant")
                        cmd.append(filename)
        else:
            print("GVCF enrichment folder does not exist")

    execute_subprocess(cmd)

def select_pass_variants(raw_vcf, nocall_fr=0.1):
    """
    Filter a vcf file. Output a vcf file with PASS positions adding a .pass to the output file
    Used since it creates the neccesasary vcf index
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php
    https://gatkforums.broadinstitute.org/gatk/discussion/13127/do-gatk4-tools-ignore-vcf-sites-marked-as-filtered-or-must-they-be-removed-from-the-file
    """
    #max_nocall=2, 

    input_vcf = os.path.abspath(raw_vcf)
    check_file_exists(input_vcf)

    raw_vcf_file_name = (".").join(input_vcf.split(".")[:-1])
     
    extension = ".pass.vcf"
    vcf_selected_output_file = raw_vcf_file_name + extension

    cmd = ["gatk", "SelectVariants", 
    "--variant", input_vcf,
    "--max-nocall-fraction", str(nocall_fr),
    "--exclude-filtered",
    "--remove-unused-alternates",
    "--output", vcf_selected_output_file]

    #"--max-nocall-number", str(max_nocall),
    execute_subprocess(cmd)

def select_pass(raw_vcf):
    """
    Homemade script
    Filter a vcf file. Output a vcf file with PASS positions adding a .pass to the output file
    """
    input_vcf = os.path.abspath(raw_vcf)
    raw_vcf_file_name = (".").join(input_vcf.split(".")[:-1])
     
    extension = ".pass.vcf"
    vcf_selected_output_file = raw_vcf_file_name + extension
    
    check_file_exists(input_vcf)
    
    with open(input_vcf, "r") as f:
        with open(vcf_selected_output_file, "w") as f1:
            for line in f:
                if line.startswith("#"):
                    f1.write(line)
                else:
                    if line.split("\t")[6] == "PASS":
                        f1.write(line)


def recalibrate_bam(args, tb=False):
    """
    BaseRecalibrator
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php
    #Recalibrate bam:
    gatk BaseRecalibrator --input my_reads.bam --reference reference.fasta --known-sites sites_of_variation.vcf \
    --known-sites another/optional/setOfSitesToMask.vcf --output recal_data.table
    ApplyBQSR
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php
    gatk ApplyBQSR --reference reference.fasta --input input.bam --bqsr-recal-file recalibration.table --output output.bam
    """
    #output = os.path.abspath(args.output)
    input_reference = os.path.abspath(args.reference)
    
    #Automate M. tuberculosis reference for aditional recalibraion positions
    if ("NC_000962.3" in input_reference) or ("h37rv" in input_reference.lower()) or ("ancestor" in input_reference):
        tb = True
        script_dir = os.path.dirname(os.path.realpath(__file__))
        reference_dir = os.path.join(script_dir, "reference")
        if ("NC_000962.3" in input_reference) or ("h37rv" in input_reference.lower()):
            reference_file = os.path.join(reference_dir, "190508_ddtb.NC_000962.3.BQSR.table")
        elif ("ancestor" in input_reference):
            reference_file = os.path.join(reference_dir, "190508_ddtb.BQSR.table")

    #group_name = output.split("/")[-1] #group_name
    sample_name = args.sample
    bam_input_dir = obtain_output_dir(args, "Bam")
    vcf_input_dir = obtain_output_dir(args, "VCF_recal")
    
    bam_input_file_name = sample_name + ".rg.markdup.sorted.bam" 
    bam_input_file = os.path.join(bam_input_dir, bam_input_file_name)

    table_output_file_name = sample_name + ".recall.table"
    table_output_file = os.path.join(vcf_input_dir, table_output_file_name)

    memory_param = "-Xmx" + str(args.memory) + "g"

    #BaseRecalibrator
    
    cmd_bqsr = ["gatk", "BaseRecalibrator",
    "--java-options", memory_param, 
    "--reference", input_reference,
    "--input", bam_input_file,
    "--output", table_output_file]

    if tb == True:
        cmd_bqsr.append("--known-sites")
        cmd_bqsr.append(reference_file)

    for root, _, files in os.walk(vcf_input_dir):
        for name in files:
            filename = os.path.join(root, name)
            if filename.endswith(".hf.pass.vcf"):
                cmd_bqsr.append("--known-sites")
                cmd_bqsr.append(filename)
    
    execute_subprocess(cmd_bqsr)

    #ApplyBQSR

    bam_output_file_name = sample_name + ".bqsr.bam"
    bam_output_file = os.path.join(bam_input_dir, bam_output_file_name)

    cmd_apply = ["gatk", "ApplyBQSR", 
    "--reference", input_reference,
    "--input", bam_input_file,
    "--bqsr-recal-file", table_output_file,
    "--output", bam_output_file]

    execute_subprocess(cmd_apply)


def samples_from_vcf(vcf_file):
    samples = subprocess.run(["bcftools", "query", "-l", vcf_file],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True)
    sample_list = samples.stdout.split("\n")[:-1]
    return sample_list


def combine_vcf(vcf_file_1, vcf_file_2, name_out=False):
    """
    Merge vcf files by position (POS) and ALT variants
    """ 
    input_vcf_1 = os.path.abspath(vcf_file_1)
    input_vcf_2 = os.path.abspath(vcf_file_2)
    #input_vcf_dir_name = os.path.dirname(vcf_file_1)
    
    if name_out == False:
        prefix = os.path.commonprefix([input_vcf_1, input_vcf_2])
        suffix = longest_common_suffix([input_vcf_1, input_vcf_2])
        output_file = prefix + "combined" + suffix
        
    else:
        output_file = os.path.abspath(name_out)
    
    header_lines_list = []
    header_lines_list_f1 = []
    #Extract filter info from header since is the only difference between headers
    with open(input_vcf_2, "r") as f2:
            for line in f2:
                if line.startswith("#"):
                    header_lines_list.append(line)
                    
    #Extract filter info from file1
    with open(input_vcf_1, "r") as f1:
        for line in f1:
            if line.startswith("##FILTER") and line not in header_lines_list:
                header_lines_list_f1.append(line)
                
    #Combine header info, addiing filter info together
    #Extract all lines starting with ##FILTER
    filter_fields = [i for i in header_lines_list if i.startswith('##FILTER')]
    #Obtain the index of the first line with index
    filter_index = header_lines_list.index(filter_fields[0])
    #Include the list within the header
    header_lines_list[filter_index:filter_index] = header_lines_list_f1

    variant_lines = []
    
    with open(input_vcf_1, "r") as f1:
        with open(input_vcf_2, "r") as f2:
            with open(output_file, "w+") as fout:
                fout.write("".join(header_lines_list))
                for line in f1:
                    if not line.startswith("#"):
                        variant_lines.append(line)
                for line in f2:
                    if not line.startswith("#") and line not in variant_lines:
                        variant_lines.append(line)
                fout.write("".join(variant_lines))


def split_vcf_saples(vcf_file, sample_list=False, nocall_fr=0.1):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php
    https://www.biostars.org/p/224702/
    #TODO: check if argument --exclude-filtered is suitable here. It would save select_pass_variants() step
    """
    
    if sample_list == False:
        #samples = subprocess.run(["bcftools", "query", "-l", vcf_file],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True)
        #sample_list = samples.stdout.split("\n")[:-1]
        sample_list = samples_from_vcf(vcf_file)
    else:
        sample_list = sample_list
    
    vcf_file_path = os.path.abspath(vcf_file)
    vcf_dir_name = os.path.dirname(vcf_file)
    vcf_file_name = vcf_file_path.split("/")[-1]
    vcf_file_extension = (".").join(vcf_file_name.split(".")[2:])

    for sample_name in sample_list:
        output_vcf_name = sample_name + "." + vcf_file_extension
        output_vcf_file = os.path.join(vcf_dir_name, output_vcf_name)
        cmd = ["gatk", "SelectVariants",
        "--max-nocall-fraction", str(nocall_fr),
        "--variant", vcf_file,
        "--sample-name", sample_name,
        "--exclude-non-variants",
        "--output", output_vcf_file]
        #"--exclude-non-variants", #remove non genotyped variants
        #"--remove-unused-alternates", #avoid poblational polymorphism
        #--preserve-alleles
        #"--keep-original-dp",
        #"--keep-original-ac",
        #"--select-type-to-include", "SNP",
        #"--select-type-to-include", "MIXED",

        if not os.path.isfile(output_vcf_file):
            execute_subprocess(cmd)

def combine_gvcf_folder(args, gvcf_input_dir, sample_list=False):
    """
    https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_CombineGVCFs.php
    #combined multi-sample gVCF:
    gatk CombineGVCFs -R reference.fasta --variant sample1.g.vcf.gz --variant sample2.g.vcf.gz -O cohort.g.vcf.gz
    """
    output = os.path.abspath(args.output)
    gvcf_input_dir = os.path.abspath(gvcf_input_dir)
    input_reference = os.path.abspath(args.reference)
    
    group_name = output.split("/")[-1] #group_name

    gvcf_output_dir = obtain_output_dir(args, "GVCF")
    gvcf_output_file = group_name + ".cohort.g.vcf"
    gvcf_output_full = os.path.join(gvcf_output_dir, gvcf_output_file)

    check_create_dir(gvcf_output_dir)

    memory_param = "-Xmx" + str(args.memory) + "g"

    cmd = ["gatk", "CombineGVCFs", 
    "--java-options", memory_param,
    "--reference", input_reference,
    "--output", gvcf_output_full]

    if os.path.isdir(gvcf_input_dir):
        for root, _, files in os.walk(gvcf_input_dir):
            for name in files:
                filename = os.path.join(root, name)
                if filename.endswith(".g.vcf"):
                    cmd.append("--variant")
                    cmd.append(filename)
    else:
        print("GVCF enrichment folder does not exist")

    execute_subprocess(cmd)

if __name__ == '__main__':
    print("#################### BAM RECALL #########################")