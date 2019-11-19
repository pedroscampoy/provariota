# PYTHON_ARGCOMPLETE_OK

import os
import sys
import re
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, obtain_output_dir, check_create_dir, execute_subprocess


"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 12 April 2019
REVISION: 
16 April 2019 : RAdd function get_bbduk_adapters() to handle adapter location
23 April 2019 : Add function execute_subprocess() to handle errors properly

TODO
    -Handle parameters

================================================================
END_OF_HEADER
================================================================
"""

"""
def get_arguments():

    parser = argparse.ArgumentParser(prog = 'bbduk_trimmer.py', description= 'Filter reads by quality and removes adapters')
    
    input_group = parser.add_argument_group('Input', 'Required files to map')

    input_group.add_argument('-1', '--read_1', dest="r1_file", metavar="R1_file[.fastq]", type=str, required=True, help='Input file for R1')
    input_group.add_argument('-2', '--read_2', dest="r2_file", metavar="R1_file[.fastq]", type=str, required=True, help='Input file for R2')

    output_group = parser.add_argument_group('Output')

    output_group.add_argument('-o', '--output', type=str, required=False, help='Output directory')
    output_group.add_argument('-s', '--sample', type=str, required=False, help='Sample name to handle output files ')

    params_group = parser.add_argument_group('Parameters', 'parameters for diferent stringent conditions')

    params_group.add_argument('-H', '--hdist', type=str, required=False, help='Set hdist parameter, default 2')
    params_group.add_argument('-k', '--kmer', type=str, required=False, help='Set k parameter, default 21')
    params_group.add_argument('-T', '--threads', type=str, dest = "threads", required=False, default=4, help='Threads to use')
    params_group.add_argument('-M', '--memory', type=str, dest = "memory", required=False, default=8, help='MAx memory to use')



    #argcomplete.autocomplete(parser)
    arguments = parser.parse_args()

    return arguments

args = get_arguments()

"""
    # Filter reads by quality and remove adapters
    # http://seqanswers.com/forums/showthread.php?t=42776
    # bbduk.sh in1=read1.fq in2=read2.fq out1=clean1.fq out2=clean2.fq ref=adapters.fa trimq=10 qtrim=w minlen=40 ktrim=r k=21 mink=11 hammingdistance=2 threads=auto tpe=f tbo=f stats=stats.txt

def get_bbduk_adapters():
    type_route = subprocess.run(["whereis", "bbduk.sh"],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
    regex = re.compile(r'\/.*\.sh')
    adapter_route = re.search(regex, type_route.stdout)

    partial_path_list = adapter_route.group().split("/")[0:-1]
    partial_path = "/".join(partial_path_list)

    final_path_adapters = os.path.join(partial_path, "resources/adapters.fa") #/home/$user/programs/bbmap/bbmap/resources/adapters.fa
        
    return final_path_adapters
    
def bbduk_trimming(args):
    """
    TODO : handle params
    """
    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)
    
    output_dir = obtain_output_dir(args, "Trimmed")

    in1_param = "in1=" + r1
    in2_param = "in2=" + r2

    sample = extract_sample(r1,r2)

    out1_param = "out1=" + output_dir + "/" + sample + "_R1.clean.fastq.gz"
    out2_param = "out2=" + output_dir + "/" + sample + "_R2.clean.fastq.gz"

    stats_param = "stats=" + output_dir + "/" + sample + "_trim.stats"

    adapter_path = "ref=" + get_bbduk_adapters()

    memory_param = "-Xmx" + str(args.memory) + "g"
    threads_param = "threads=" + str(args.threads)

    check_create_dir(output_dir)

    #bbduk.sh
    cmd = ["bbduk.sh", memory_param, in1_param, in2_param, out1_param, out2_param, adapter_path,
        "trimq=15" , "qtrim=rl", "minlen=40",
        "ktrim=r", "k=21", "mink=11", "hammingdistance=2",
        threads_param, "tpe", "tbo", stats_param]
    
    execute_subprocess(cmd)



if __name__ == '__main__':
    print("#################### TRIMMING #########################")
    #print(args)

    #get_bbduk_adapters()
    #bbduk_trimming(args)

#python bbduk_trimmer.py -1 ../RAW/AL14621_R1.fastq.gz -2 ../RAW/AL14621_R2.fastq.gz -o .
#python bbduk_trimmer.py -1 /home/pedro/analysis/Mixed/RAW/1mixed_S58_R1_001.fastq.gz -2 /home/pedro/analysis/Mixed/RAW/1mixed_S58_R2_001.fastq.gz -o /home/pedro/analysis/Mixed