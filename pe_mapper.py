#!/usr/bin/env python

import os
import gzip
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, get_picard_path, execute_subprocess, check_remove_file


"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 27 March 2019
REVISION:

TODO

================================================================
END_OF_HEADER
================================================================
"""


def get_arguments():

    parser = argparse.ArgumentParser(prog = 'pe_mapper.py', description= 'Creates an index, map reads and convert sam to ordered bam')
    
    input_group = parser.add_argument_group('Input', 'Required files to map')

    input_group.add_argument('-1', '--read_1', dest="r1_file", metavar="R1_file[.fastq]", type=str, required=True, help='Input file for R1')
    input_group.add_argument('-2', '--read_2', dest="r2_file", metavar="R1_file[.fastq]", type=str, required=True, help='Input file for R2')
    input_group.add_argument('-r', '--reference', type=str, required=True, help='File to map against')

    output_group = parser.add_argument_group('Output')

    output_group.add_argument('-o', '--output', type=str, required=False, help='Output file name, default same as input directory')
    output_group.add_argument('-s', '--sample', type=str, required=False, help='Sample name to handle output files ')

    params_group = parser.add_argument_group('Parameters', 'map and arguments for each mapper')

    mapper_exclusive = params_group.add_mutually_exclusive_group()

    mapper_exclusive.add_argument("--bowtie2",  dest = "bowtie2_mapper", action="store_false", required= False, help="Uses bowtie2 to map reads (Default)")
    mapper_exclusive.add_argument("--bwa",  dest = "bwa_mapper", action="store_true", required= False, help="Uses bwa to map reads")

    params_group.add_argument('-a', '--extensive_mapping', action="store_true", required=False, help='Use extensive mapping (Default False)')
    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=4, help='Threads to use')
    params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=8, help='MAx memory to use')


    arguments = parser.parse_args()

    return arguments



    # Determine mapper used in pipeline and get the command
    # bwa index ref.fa && bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
    # bowtie2-build $database $database &&  bowtie2 -1 $R1 -2 $R2 -S $output_dir/$sample.sam 
    # -q --very-sensitive-local $a_mapping -p $threads -x $database
    
def bowtie2_mapping(args):
    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)
    reference = os.path.abspath(args.reference)

    sample = extract_sample(r1,r2)
    output_dir = obtain_output_dir(args, "Bam")
    sample_name = sample + ".sam"
    output_file = os.path.join(output_dir, sample_name)

    check_create_dir(output_dir)

    if args.extensive_mapping:
        extensive_command = "-a"
    else:
        extensive_command = ""
    #bowtie2 index
    cmd_index = ["bowtie2-build", reference, reference]
    execute_subprocess(cmd_index)
    
    #bowtie map
    cmd_map = ["bowtie2", "-1", r1, "-2", r2, "-S", output_file, "-q", "--very-sensitive-local", "-p", str(args.threads), "-x", reference, extensive_command]
    execute_subprocess(cmd_map)



def bwa_mapping(args):
    """
    #Store output in a file when it is outputted in stdout
    https://stackoverflow.com/questions/4965159/how-to-redirect-output-with-subprocess-in-python
    """
    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)
    reference = os.path.abspath(args.reference)

    sample = extract_sample(r1,r2)
    output_dir = obtain_output_dir(args, "Bam")
    sample_name = sample + ".sam"
    output_file = os.path.join(output_dir, sample_name)

    check_create_dir(output_dir)
    
    cmd_index = ["bwa", "index", reference]
    execute_subprocess(cmd_index)
    
    cmd_map = ["bwa", "mem", "-t", str(args.threads), "-o", output_file, reference, r1, r2]
    execute_subprocess(cmd_map)
    """
    Create file whew it outputs thoug stdout --> Easier with -o param
    with open(output_file, "w") as outfile:
        #map reads and save it in th eoutput file
        subprocess.run(["bwa", "mem", "-t", str(args.threads), reference, r1, r2], 
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)
    """

def add_SG(args, input_bam, output_bg_sorted):
    """
    @MN00227:45:000H255J3:1:11102:21214:1110 1:N:0:18
    @NS500454:48:HKG57BGXX:1:11101:17089:1032 2:N:0:TCCTGAGC+TCTTACGC
    @NS500454:27:HJJ32BGXX:1:11101:12392:1099 1:N:0:2

    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:
    <is filtered>:<control number>:<sample number | barcode1'+barcode2'>
    ID = Read group identifier {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE} 
    PU = Platform Unit #optional
    SM = Sample
    PL = Platform/technology used to produce the read (ILLUMINA, SOLID, LS454, HELICOS and PACBIO)
    LB = DNA preparation library identifier
    """
    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)

    sample = extract_sample(r1,r2)

    with gzip.open(r1) as f:
        first_line = f.readline().strip().decode()
    #print(first_line)
    first_line_list = first_line.split(":")
 
    rg_id = ".".join([first_line_list[2],first_line_list[3],first_line_list[-1]])
    rg_pu = ".".join([first_line_list[2],first_line_list[3],first_line_list[-1]])
    rg_sm = sample
    rg_pl = "ILLUMINA"
    rg_lb = "lib_" + sample

    rg_id_param = "RGID=" + rg_id
    rg_pu_param = "RGPU=" + rg_pu
    rg_sm_param = "RGSM=" + rg_sm
    rg_pl_param = "RGPL=" + rg_pl
    rg_lb_param = "RGLB=" + rg_lb

    picard_jar = get_picard_path()

    input_param = "INPUT=" + input_bam
    output_param = "OUTPUT=" + output_bg_sorted


    # java -jar picard.jar AddOrReplaceReadGroups \ 
    # INPUT=reads.bam \ OUTPUT=reads_addRG.bam \ RGID=H0164.2 \ #be sure to change from default of 1
    # RGLB= library1 \ RGPL=illumina \ RGPU=H0164ALXX140820.2 \ RGSM=sample1 \ 
    # SORT_ORDER=coordinate \ CREATE_INDEX=true

    cmd = ["java", "-jar", picard_jar, "AddOrReplaceReadGroups", 
    input_param, output_param, rg_id_param, rg_lb_param, rg_pl_param, rg_pu_param, rg_sm_param,
    "SORT_ORDER=coordinate"]
    execute_subprocess(cmd)

def sam_to_index_bam(args):
    # input_sam_path = os.path.abspath(input_sam)
    # if output_bam == "inputdir":
    #     output_bam = os.path.dirname(input_sam_path)
    # else:
    #     output_bam = output_bam

    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)

    sample = extract_sample(r1,r2)
    output_dir = obtain_output_dir(args, "Bam")
    sample_name = sample + ".sam"
    input_sam_path = os.path.join(output_dir, sample_name)

    input_name = (".").join(os.path.basename(input_sam_path).split(".")[:-1])

    output_bam_name = input_name + ".bam"
    output_bam_path = os.path.join(output_dir, output_bam_name)

    output_bg_sorted_name = input_name + ".rg.sorted.bam"
    output_bg_sorted_path = os.path.join(output_dir, output_bg_sorted_name)

    check_create_dir(output_dir)
    """
    #sam to bam: samtools view -Sb $input_file -o $output_dir/$sample.bam
    with open(output_bam_path, "w") as outfile:
        #map reads and save it in th eoutput file
        subprocess.run(["samtools", "view", "-Sb", input_sam_path], 
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)
    """
    cmd = ["samtools", "view", "-Sb", input_sam_path, "-o", output_bam_path, "--threads", str(args.threads)]
    execute_subprocess(cmd)

    check_remove_file(input_sam_path)

    add_SG(args, output_bam_path, output_bg_sorted_path)

    check_remove_file(output_bam_path)

    """
    output_sorted_name = input_name + ".sorted.bam"
    output_sorted_path = os.path.join(output_dir, output_sorted_name)
    if os.path.exists(output_sorted_path):
        os.remove(output_sorted_path)
    
    REPLACED BY PICARD TOOLS
    #samtools sort: samtools sort $output_dir/$sample".sorted.bam" -o $output_dir/$sample".sorted.bam"
    subprocess.run(["samtools", "sort", output_bam_path, "-o", output_sorted_path], 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    os.remove(output_bam_path)

    #samtools index: samtools index $output_dir/$sample".sorted.bam"
    subprocess.run(["samtools", "index", output_sorted_path], 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    """
