
import os
import sys
import re
import subprocess
import pandas as pd
import numpy as np


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

def check_file_exists(file_name):
    """
        Check file exist and is not 0 Kb, if not program exit.
    """
    file_info = os.stat(file_name) #Retrieve the file info to check if has size > 0

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        print(RED + BOLD + "File: %s not found or empty\n" % file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)


def check_remove_file(file_name):
    """
    Check file exist and remove it.
    """
    if os.path.exists(file_name):
        os.remove(file_name)
    

def import_to_pandas(file_table, header=False, sep='\t'):
    if header == False:
        #exclude first line, exclusive for vcf outputted by PipelineTB
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        #Use first line as header
        dataframe = pd.read_csv(file_table, sep=sep, header=0)
    
    return dataframe


def extract_sample(R1_file, R2_file):
    """
    Extract sample from R1, R2 files.
    """
    basename_R1 = os.path.basename(R1_file)
    basename_R2 = os.path.basename(R2_file)

    sample_name_R = os.path.commonprefix([basename_R1, basename_R2])
  
    long_suffix = re.search('_S.*', sample_name_R)
    short_suffix = re.search('_R.*', sample_name_R)
    bar_suffix = re.search('_$', sample_name_R)
    
    if long_suffix:
        match = long_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif short_suffix:
        match = short_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif bar_suffix:
        match = bar_suffix.group()
        sample_name = sample_name_R.rstrip("_")
    else:
        sample_name = sample_name_R

    return sample_name


def obtain_output_dir(args, subfolder=None):
    """
    Returns output folder and output file depending on the output supplied.
    """
    if args.output != None:
        output_dir_arg = os.path.abspath(args.output)
        output_dir = os.path.join(output_dir_arg, subfolder)
    elif args.r1_file:
        r1 = os.path.abspath(args.r1_file)
        output_dir_arg = os.path.dirname(r1)
        output_dir = os.path.join(output_dir_arg, subfolder)
    elif args.input_bam:
        bam = os.path.abspath(args.input_bam)
        output_dir_arg = os.path.dirname(bam)
        output_dir = os.path.join(output_dir_arg, subfolder)
    return output_dir
    
def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def get_picard_path():
    type_route = subprocess.run(["whereis", "picard.jar"],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
    regex = re.compile(r'\/.*\.jar')
    picard_route = re.search(regex, type_route.stdout)

    return picard_route.group()


def get_snpeff_path():
    type_route = subprocess.run(["whereis", "snpEff.jar"],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
    snpeff_route = type_route.stdout.split(" ")
    snpeff_jar = False
    snpeff_config = False
    for executable in snpeff_route:
        executable = executable.strip()
        if executable.endswith(".jar"):
            snpeff_jar = executable
        elif executable.endswith(".config"):
            snpeff_config = executable
    
    if snpeff_jar and snpeff_config:
        return snpeff_config, snpeff_jar
    else:
        print("some snpEff executables are missing")
        sys.exit(1)


def execute_subprocess(cmd):
    """
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """
    if cmd[0] == "java":
        prog = cmd[2].split("/")[-1] + " " + cmd[3]
        param = cmd[4:]
    elif cmd[0] == "samtools" or cmd[0] == "bwa" or cmd[0] == "gatk":
        prog = " ".join(cmd[0:2])
        param = cmd[3:]
    else:
        prog = cmd[0]
        param = cmd[1:]
    
    try:
        command = subprocess.run(cmd , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if command.returncode == 0:
            print(GREEN + DIM + "Program %s successfully executed" % prog + END_FORMATTING)
        else:
            print (RED + BOLD + "Command %s FAILED\n" % prog + END_FORMATTING
                + BOLD + "WITH PARAMETERS: " + END_FORMATTING + " ".join(param) + "\n"
                + BOLD + "EXIT-CODE: %d\n" % command.returncode +
                "ERROR:\n" + END_FORMATTING + command.stderr.decode().strip())
    except OSError as e:
        sys.exit(RED + BOLD + "failed to execute program '%s': %s" % (prog, str(e)) + END_FORMATTING)


def extract_read_list_legacy(input_dir):
    """
    Search files in a directory sort by name and extract comon name of R1 and R2
    with extract_sample() function
    190615 - Limit only parent folder, not subdirectories
    """
    input_dir = os.path.abspath(input_dir)
    r1_list = []
    r2_list = []
    for root, _, files in os.walk(input_dir):
        if root == input_dir: # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                is_fasta = re.match(r'.*\.f(ast)*[aq](\.gz)*',name)
                r1 = re.match(r'.*(_R1_|_1|_1_|_R1).*\.f(ast)*[aq](\.gz)*$',name)
                r2 = re.match(r'.*(_R2_|_2|_2_|_R2).*\.f(ast)*[aq](\.gz)*$',name)
                if is_fasta:
                    if r1:
                        r1_list.append(filename)
                    elif r2:
                        r2_list.append(filename)
                    else:
                        print(RED + "ERROR, file is not R1 nor R2" + END_FORMATTING)
                        sys.exit(1)
    r1_list = sorted(r1_list)
    r2_list = sorted(r2_list)
    return r1_list, r2_list

def extract_read_list(input_dir):
    """
    Search files in a directory sort by name and extract comon name of R1 and R2
    with extract_sample() function
    190615 - Limit only parent folder, not subdirectories
    """
    input_dir = os.path.abspath(input_dir)
    all_fasta = []
    r1_list = []
    r2_list = []
    for root, _, files in os.walk(input_dir):
        if root == input_dir: # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                is_fasta = re.match(r'.*\.f(ast)*[aq](\.gz)*',filename)
                if is_fasta:
                    all_fasta.append(filename)
    all_fasta = sorted(all_fasta)
    if len(all_fasta) % 2 == 0:
        for index, fasta_file in enumerate(all_fasta):
            if index % 2 == 0:
                r1_list.append(fasta_file)
            elif index % 1 == 0:
                r2_list.append(fasta_file)          
    else:
        print('ERROR: The number of fastq sequence are not paired')
        
    r1_list = sorted(r1_list)
    r2_list = sorted(r2_list)
    
    return r1_list, r2_list


def extract_sample_list():
    #sample_list = []
    # for r1, r2 in zip(r1_list, r2_list):
    #     sample = extract_sample(r1, r2)
    #     sample_list.append(sample)
    pass

def return_codon_position(number):
    position = number % 3
    if position == 0:
        position = 3
    print("number= %s, pos= %s" % (number,position))

def file_to_list(file_name):
    list_F = []
    file_name_abs = os.path.abspath(file_name)
    with open(file_name_abs, "r") as f:
        for line in f:
            list_F.append(line.strip())
    return list_F


def get_coverage(args, input_bam, output_fmt="-d"):
    """
    #Calculate genome coverage at each position using bedtools and an input bam
    https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
    """
    #reference = os.path.abspath(args.reference)

    input_bam = os.path.abspath(input_bam)
    input_bam_base = os.path.basename(input_bam)

    sample = input_bam_base.split(".")[0]
    output_dir = obtain_output_dir(args, "Coverage")
    sample_name = sample + ".cov"
    output_file = os.path.join(output_dir, sample_name)

    check_create_dir(output_dir)

    #execute_subprocess(cmd)
    with open(output_file, "w") as outfile:
        #calculate coverage and save it in th eoutput file
        subprocess.run(["genomeCoverageBed", "-ibam", input_bam, output_fmt], 
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)

def calculate_cov_stats(file_cov):
    df = pd.read_csv(file_cov, sep="\t", names=["#CHROM", "POS", "COV" ])
    unmmaped_pos = len(df.POS[df.COV == 0].tolist())
    pos_0_10 = len(df.POS[(df.COV > 0) & (df.COV <= 10)].tolist())
    pos_10_20 = len(df.POS[(df.COV > 10) & (df.COV <= 20)].tolist())
    pos_high20 = len(df.POS[(df.COV > 20)].tolist())
    total_pos = df.shape[0]
    unmmaped_prop = "%.2f" % ((unmmaped_pos/total_pos)*100)
    prop_0_10 = "%.2f" % ((pos_0_10/total_pos)*100)
    prop_10_20 = "%.2f" % ((pos_10_20/total_pos)*100)
    prop_high20 = "%.2f" % ((pos_high20/total_pos)*100)
    
    mean_cov = "%.2f" % (df.COV.mean())
    
    return mean_cov, unmmaped_prop, prop_0_10, prop_10_20, prop_high20

def obtain_group_cov_stats(directory, low_cov_threshold=20, unmmaped_threshold=20):
    directory_path = os.path.abspath(directory)
    
    if directory_path.endswith("Coverage"):
        file_name = directory_path.split("/")[-2]
    else:
        file_name = "samples"

    output_file_name = file_name + ".coverage.tab"
    output_file = os.path.join(directory_path,output_file_name)

    saples_low_covered = []

    with open(output_file, "w") as outfile:
            outfile.write("#SAMPLE" + "\t" + "MEAN_COV" + "\t" + "UNMMAPED_PROP" + "\t" + "COV1-10X" + "\t" + "COV10-20X" + "\t" + "COV>20X" + "\t" + "\n")
            #print("#SAMPLE" + "\t" + "MEAN_COV" + "\t" + "UNMMAPED_PROP" + "\t" + "COV1-10X" + "\t" + "COV10-20X" + "\t" + "COV>20X" + "\t" + "\n")
            for root, _, files in os.walk(directory_path):
                for name in files:
                    filename = os.path.join(root, name)
                    file_name_cov = os.path.basename(filename)
                    sample = file_name_cov.split(".")[0]
                    if filename.endswith("cov") and (os.path.getsize(filename) > 0):
                        coverage_stats = calculate_cov_stats(filename)
                        mean_cov = coverage_stats[0]
                        unmmaped_prop = coverage_stats[1]
                        if float(mean_cov) < low_cov_threshold or float(unmmaped_prop) > unmmaped_threshold:
                            saples_low_covered.append(sample)
                        outfile.write(sample + "\t" + ("\t").join(coverage_stats) + "\n")
                        #print((sample + "\t" + ("\t").join(coverage_stats)) + "\n")

    return saples_low_covered

def edit_sample_list(file_list, sample_list):
    with open(file_list, 'r') as f:
        content = f.read()
        content_list = content.split('\n')
        while '' in content_list : content_list.remove('')
        
    with open (file_list, 'w+') as fout:
            for line in content_list:
                if line not in sample_list:
                    fout.write(line + "\n")

def remove_low_covered_mixed(output_dir, sample_list, type_remove):
    output_dir = os.path.abspath(output_dir)
    group = output_dir.split("/")[-1]
    uncovered_dir = os.path.join(output_dir, type_remove) #Uncovered or Mixed
    check_create_dir(uncovered_dir)

    sample_list_file = os.path.join(output_dir, "sample_list.txt")
    
    for root, _, files in os.walk(output_dir):
        #Any previous file created except for Table for mixed samples
        # and Species for both uncovered and mixed
        if root.endswith('GVCF_recal') or root.endswith('Coverage') \
        or root.endswith('VCF') or root.endswith('VCF_recal') \
        or root.endswith('Bam') or root.endswith('GVCF') \
        or root.endswith('Table'):
            for name in files:
                filename = os.path.join(root, name)
                for sample_low in sample_list:
                    sample_dot = sample_low + "." #Adapt name to the posibility that two samples starts with the same name
                    if name.startswith(sample_dot):
                        if os.path.isfile(sample_list_file):
                            os.remove(filename)
 
        #Place low covered samples in a specific folder to analize them with different parameters
        if root.endswith(group):
            for name in files:
                filename = os.path.join(root, name)
                for sample_low in sample_list:
                    sample_lowbar = sample_low + "_"
                    if name.startswith(sample_lowbar) and name.endswith("fastq.gz"):
                        dest_uncovered_path = os.path.join(uncovered_dir, name)
                        if os.path.isfile(sample_list_file):
                            os.rename(filename, dest_uncovered_path)
    if os.path.isfile(sample_list_file):
        edit_sample_list(sample_list_file, sample_list)


def clean_unwanted_files(args):
    Trimmed_dir = ""
    for root, _, files in os.walk(args.output):
        if 'Trimmed' in root:
            Trimmed_dir = root
        for name in files:
            filename = os.path.join(root, name)
            if root.endswith("Bam") and not "bqsr" in filename:
                print("Removed: " + filename)
                os.remove(filename)
            #elif filename.endswith("cohort.g.vcf") or filename.endswith("cohort.g.vcf.idx"):
            #    print("Removed: " + filename)
            #    os.remove(filename)
            elif root.endswith("Annotation") and (filename.endswith("annot.genes.txt") or filename.endswith(".vcf") or filename.endswith(".annot.html")):
                print("Removed: " + filename)
                os.remove(filename)
            elif root.endswith("Trimmed"):
                print("Removed: " + filename)
                os.remove(filename)
                
    if Trimmed_dir:
        print("Removed folder: " + Trimmed_dir)
        os.rmdir(Trimmed_dir)
                
def longest_common_suffix(list_of_strings):
    """
    Return the longest common suffix in a list of strings
    Adapted from https://gist.github.com/willwest/ca5d050fdf15232a9e67
    """
    reversed_strings = [s[::-1] for s in list_of_strings]
    reversed_lcs = os.path.commonprefix(reversed_strings)
    lcs = reversed_lcs[::-1]
    return lcs

def list_to_bed(input_list, output_dir, output_file_name, reference="CHROM"):
    """
    Turn a list into a bed file with start and end position having the same value
    """
    output_dir = os.path.abspath(output_dir)
    
    output_bed_file = output_file_name + ".bed"
    
    final_output_path = os.path.join(output_dir, output_bed_file)
    
    with open (final_output_path, 'w+') as f:
        for position in input_list:
            line = ("\t").join([reference, str(position), str(position)]) + "\n"
            f.write(line)

def count_lines(input_file):
    with open(input_file, 'r') as f:
        content = f.read()
        content_list = content.split('\n')
        while '' in content_list : content_list.remove('')
    return len(content_list)

def check_reanalysis(output_dir):
    output_dir = os.path.abspath(output_dir)
    #group = output_dir.split("/")[-1]
    
    bam_dir = os.path.join(output_dir, "Bam")
    vcf_dir = os.path.join(output_dir, "VCF")
    gvcf_dir = os.path.join(output_dir, "GVCF")
    gvcfr_dir = os.path.join(output_dir, "GVCF_recal")
    vcfr_dir = os.path.join(output_dir, "VCF_recal")
    cov_dir = os.path.join(output_dir, "Coverage")
    table_dir = os.path.join(output_dir, "Table")
    
    previous_files = [bam_dir, vcf_dir, gvcf_dir, gvcfr_dir]
    
    #check how many folders exist
    file_exist = sum([os.path.exists(x) for x in previous_files]) #True = 1, False = 0
    
    #Handle reanalysis: First time; reanalysis o reanalysis with aditional samples
    if file_exist > 0: #Already analysed
        
        samples_analyzed = os.listdir(bam_dir)
        samples_analyzed = len([ x for x in samples_analyzed if ".bai" not in x and "bqsr" in x])

        samples_fastq = os.listdir(output_dir)
        samples_fastq = len([ x for x in samples_fastq if x.endswith('fastq.gz')]) / 2
        
        if samples_analyzed >= samples_fastq:
            print(MAGENTA + "\nPREVIOUS ANALYSIS DETECTED, NO NEW SEQUENCES ADDED\n" + END_FORMATTING)
        
        else:
            print(MAGENTA + "\nPREVIOUS ANALYSIS DETECTED, NEW SEQUENCES ADDED\n" + END_FORMATTING)
            for root, _, files in os.walk(output_dir):
                    if root ==  gvcf_dir or root == gvcfr_dir or root == vcfr_dir:
                        for name in files:
                            filename = os.path.join(root, name)
                            if (("GVCF_recal" in filename) or ("/VCF_recal" in filename)) and "cohort" in filename and samples_analyzed < 100:
                                os.remove(filename)
                            elif "cohort" in filename and "/GVCF/" in filename:
                                os.remove(filename)
                    elif root == vcf_dir or root == table_dir:
                        for name in files:
                            filename = os.path.join(root, name)
                            if "cohort" in filename or filename.endswith(".bed") or filename.endswith(".tab"):
                                os.remove(filename)
                    elif root == cov_dir:
                        for name in files:
                            filename = os.path.join(root, name)
                            if "coverage.tab" in filename:
                                os.remove(filename)
                            if "poorly_covered.bed" in filename and samples_analyzed < 100:
                                os.remove(filename)
            #print(file_exist, samples_analyzed, samples_fastq)

def extrach_variants_summary(vcf_table, distance=15, quality=10 ):
    sample = vcf_table.split("/")[-1].split(".")[0]
    
    df = pd.read_csv(vcf_table, sep="\t", header=0)
    
    total_snp = len(df[df.TYPE == "SNP"].index)
    total_indels = len(df[df.TYPE == "INDEL"].index)
    total_homozygous = len(df[(df.TYPE == "SNP") & (df.gt0 == 1)].index)
    total_heterozygous = len(df[(df.TYPE == "SNP") & (df.gt0 == 0)].index)
    median_allele_freq = "%.2f" % (df.AF[df.TYPE == "SNP"].median())
    mean_allele_freq = "%.2f" % (df.AF[df.TYPE == "SNP"].mean())
    
    distance = distance
    QD = quality
    position_to_filter = df['POS'][((df.snp_left_distance <= distance)|
                                (df.snp_right_distance <= distance)|
                                (df.window_10 >= 2)|
                                (df.AF <= 0.0) |
                                (df.len_AD > 2) |
                                (df.TYPE != "SNP") |
                                (df.QD <= QD) |
                                (df.highly_hetz == True) |
                                (df.poorly_covered == True) |
                                (df.non_genotyped == True) |
                                (df.is_polymorphic == True))].tolist()
    
    filtered_df = df[~df.POS.isin(position_to_filter)]
    
    filtered_df_htz = filtered_df[filtered_df.gt0 == 0]
    
    ftotal_snp = len(filtered_df[filtered_df.TYPE == "SNP"].index)
    ftotal_homozygous = len(filtered_df[(filtered_df.TYPE == "SNP") & (filtered_df.gt0 == 1)].index)
    ftotal_heterozygous = len(filtered_df[(filtered_df.TYPE == "SNP") & (filtered_df.gt0 == 0)].index)
    fmedian_allele_freq = "%.2f" % (filtered_df.AF[filtered_df.TYPE == "SNP"].median())
    fmean_allele_freq = "%.2f" % (filtered_df.AF[filtered_df.TYPE == "SNP"].mean())
    fmean_allele_freq_htz = "%.2f" % (filtered_df_htz.AF[filtered_df_htz.TYPE == "SNP"].mean())
    
    output = [sample,
              total_snp,
              total_indels,
              total_homozygous,
              total_heterozygous,
              median_allele_freq,
              mean_allele_freq,
              ftotal_snp,
              ftotal_homozygous,
              ftotal_heterozygous,
              fmedian_allele_freq,
              fmean_allele_freq,
              fmean_allele_freq_htz]
    output = [str(x) for x in output]
    
    return "\t".join(output)

def vcf_stats(folder_table, distance=15, quality=10):
    
    out_file = os.path.join(folder_table, "vcf_stat.tab")
    mixed_samples = []
    
    with open(out_file, 'w+') as fout:
        fout.write("\t".join(["SAMPLE", 
                              "#SNP", 
                              "#INDELS", 
                              "#HOMOZ_SNP", 
                              "#HETZ_SNP", 
                              "MEDIAN_AF_SNP", 
                              "MEAN_AF_SNP", 
                              "#FSNP", 
                              "#FHOMOZ_SNP", 
                              "#FHETZ_SNP", 
                              "FMEDIAN_AF_SNP",
                              "FMEAN_AF_SNP",
                              "FMEAN_AF_SNP_HTZ"]))
        fout.write("\n")
        for root, _, files in os.walk(folder_table):
            for name in files:
                filename = os.path.join(root, name)
                if filename.endswith("raw.tab"):
                    line = extrach_variants_summary(filename)
                    line_split = line.split("\t")
                    sample = line_split[0]
                    htz_filtered = line_split[9]
                    if int(htz_filtered) > 100:
                        mixed_samples.append(sample)
                    fout.write(line)
                    fout.write("\n")
    return mixed_samples