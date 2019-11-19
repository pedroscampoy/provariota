#!/usr/bin/env python

import os
import pandas as pd
import argparse
import sys
import subprocess
from sklearn.metrics import jaccard_similarity_score, pairwise_distances, accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as ssd #pdist

from misc import check_file_exists, import_to_pandas, check_create_dir
from vcf_process import import_VCF42_cohort_pandas

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
    
    parser.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all vcf files')
    parser.add_argument('-s', '--sample_list', type=str, required=False, help='File with sample names to analyse instead of all samples')
    parser.add_argument("-r", "--recalibrate", required= False, type=str, default=False, help='Pipeline folder where Bam and VCF subfolders are present')
    parser.add_argument("-R", "--reference", required= False, type=str, default=False, help='Reference fasta file used in original variant calling')

    parser.add_argument('-o', '--output', type=str, required=True, help='Name of all the output files, might include path')

    arguments = parser.parse_args()

    return arguments

def blank_database():
    new_pandas_ddtb = pd.DataFrame(columns=['Position','N', 'Samples'])
    return new_pandas_ddtb

def import_VCF4_to_pandas(vcf_file, sep='\t'):
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
        dataframe['POS'] = dataframe['POS'].astype(int)
        
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe

def recheck_variant(format_sample):
    #GT:AD:DP:GQ:PGT:PID:PL:PS
    list_format = format_sample.split(":")
    gt = list_format[0]
    #gt0 = gt[0]
    #gt1 = gt[1]
    ad = list_format[1]
    ref = int(ad.split(',')[0])
    alt = max(int(x) for x in ad.split(',')[0:])
    
    if gt == "0/0":
        value = 0
    elif gt == "1/1":
        value = 1
    else:
        if gt == "./.":
            value = "!"
        elif "2" in gt:
            value = "!"
        elif (ref > alt):
            value = 0
        elif (alt > ref):
            value = 1
        else:
            value = "!"
            
    return value

def recheck_variant_mpileup(reference_file, position, sample, bam_folder):
    #Find reference name
    with open(reference_file) as f:
        reference = f.readline().split(" ")[0].strip(">").strip()
    #Identify correct bam
    for root, _, files in os.walk(bam_folder):
        for name in files:
            filename = os.path.join(root, name)
            if name.startswith(sample) and name.endswith(".bqsr.bam"):
                bam_file = filename
    #format position for mpileuo execution (NC_000962.3:632455-632455)
    position = reference + ":" + str(position) + "-" + str(position)
    
    #Execute command and retrieve output
    cmd = ["samtools", "mpileup", "-f", reference_file, "-aa", "-r", position, bam_file]
    text_mpileup = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
    
    #Extract 5th column to find variants
    variant = text_mpileup.stdout.split()[4]
    var_list = list(variant)
    
    most_freq_var = max(set(var_list), key = var_list.count).upper()
        
    if most_freq_var == "." or most_freq_var == "," or most_freq_var == "*":
        return 0
    else:
        return 1

def identify_nongenotyped_mpileup(reference_file, row_position, sample_list_matrix, list_presence, bam_folder):
    """
    Replace nongenotyped ("!") with the most abundant genotype
    """
    #mode = max(set(list_presence), key = list_presence.count)
    
    count_ng = list_presence.count("!")
    sample_number = len(list_presence)
    
    if "!" not in list_presence:
        return list_presence
    elif count_ng/sample_number > 0.2:
        return 'delete'
    else:
        indices_ng = [i for i, x in enumerate(list_presence) if x == "!"]
        for index in indices_ng:
            #print(reference_file, row_position, sample_list_matrix[index], bam_folder)
            list_presence[index] = recheck_variant_mpileup(reference_file, row_position, sample_list_matrix[index], bam_folder)
        #new_list_presence = [mode if x == "!" else x for x in list_presence]
        return list_presence

def extract_recalibrate_params(pipeline_folder, reference=False):
    cohort_file = ""
    pipeline_folder = os.path.abspath(pipeline_folder)
    for root, dirs, _ in os.walk(pipeline_folder):
        #print(pipeline_folder, root)
        if root == pipeline_folder:
            for directory in dirs:
                subfolder = os.path.join(root, directory)
                if subfolder.endswith("/VCF"):
                    for file_vcf in os.listdir(subfolder):
                        if file_vcf.endswith("cohort.combined.hf.vcf"):
                            cohort_file = os.path.join(subfolder, file_vcf)
                            if reference == False:
                                with open(cohort_file, 'r') as f:
                                    for line in f:
                                        if line.startswith("#"):
                                            if "--reference " in line:
                                                reference_file = line.split("--reference ")[1].strip().split(" ")[0].strip()
                            else:
                                reference_file = reference
                elif subfolder.endswith("/Bam"):
                    bam_folder = subfolder
    if cohort_file:
        return (cohort_file, bam_folder, reference_file)
    else:
        print(RED + "cohort.combined.hf.vcf not found, wait for pipeline to finish" + END_FORMATTING)
        sys.exit(1)
                    

def recalibrate_ddbb_vcf(snp_matrix_ddbb, vcf_cohort, bam_folder, reference_file):
    
    vcf_cohort = os.path.abspath(vcf_cohort)
    #snp_matrix_ddbb = os.path.abspath(snp_matrix_ddbb)
    
    df_matrix = snp_matrix_ddbb
    df_cohort = import_VCF42_cohort_pandas(vcf_cohort)
    
    sample_list_matrix = df_matrix.columns[3:]
    n_samples = len(sample_list_matrix)
    #sample_list_cohort = df_cohort.columns.tolist()[9:]
    
    list_index_dropped = []
    #Iterate over non unanimous positions 
    for index, data_row in df_matrix[df_matrix.N < n_samples].iloc[:,3:].iterrows():
        #Extract its position
        row_position = int(df_matrix.loc[index,"Position"])
        #print(data_row.values)
        #Use enumerate to retrieve column index (column ondex + 3)
        presence_row = [recheck_variant(df_cohort.loc[df_cohort.POS == row_position, df_matrix.columns[n + 3]].item()) \
                           for n,x in enumerate(data_row)]
        #print(presence_row, row_position)
        #Resolve non genotyped using gvcf files
        new_presence_row = identify_nongenotyped_mpileup(reference_file, row_position, sample_list_matrix, presence_row, bam_folder)
        
        #find positions with 20% of nongenotyped and delete them OR
        #reasign positions without nongenotyped positions 
        if new_presence_row == 'delete':
            list_index_dropped.append(index)
        else:
            df_matrix.iloc[index, 3:] = new_presence_row
            df_matrix.loc[index, 'N'] = sum(new_presence_row)
        #print(new_presence_row)
        #print("\n")
    #Remove all rows at once to avoid interfering with index during for loop
    df_matrix.drop(index=list_index_dropped, axis=0, inplace=True)
    
    return df_matrix


def ddtb_add(input_folder, output_filename, recalibrate=False, sample_filter=False, vcf_suffix=".combined.hf.SNP.final.vcf" ):
    directory = os.path.abspath(input_folder)
    output_filename = os.path.abspath(output_filename)

    #Make sure output exist to force change name
    if os.path.isfile(output_filename):
        print(YELLOW + "ERROR: " + BOLD + "output database EXIST, choose a different name or manually delete" + END_FORMATTING)
        sys.exit(1)

    final_ddbb = blank_database()
    sample_filter_list = []


    #Handle sample filter
    if sample_filter == False:
        sample_filter_list = [x.split(".")[0] for x in os.listdir(input_folder) if x.endswith(vcf_suffix)]
    else:
        if os.path.isfile(sample_filter):
            with open(sample_filter, 'r') as f:
                for line in f:
                    sample_filter_list.append(line.strip())
        else:
            "Sample file don't exist"
            sys.exit(1)
    
    print(sample_filter_list)

    if len(sample_filter_list) < 1:
        print("prease provide 2 or more samples")
        sys.exit(1)

    #print("Previous final database contains %s rows and %s columns\n" % final_ddbb.shape)
    print("The directory selected is: %s" % directory)
    

    all_samples = 0
    new_samples = 0
    for filename in os.listdir(directory):
        if not filename.startswith('.') and filename.endswith(vcf_suffix):
            
            all_samples = all_samples + 1
            positions_shared = []
            positions_added = []
            
            sample = filename.split(".")[0] #Manage sample name

            if sample in sample_filter_list:
                print("\nThe file is: %s" % filename)

            
                file = os.path.join(directory, filename) #Whole file path
                check_file_exists(file) #Manage file[s]. Check if file exist and is greater than 0

                new_sample = import_VCF4_to_pandas(file) #Import files in annotated vcf format

                #Check if sample exist
                ######################
                if sample not in final_ddbb.columns.tolist():
                    print("Adding new sample %s to %s" % (sample, os.path.basename(output_filename)))
                    new_samples = new_samples + 1
                    new_colum_index = len(final_ddbb.columns) #extract the number of columns to insert a new one
                    #final_ddbb[sample] = sample #adds a new column but fills all blanks with the value sample
                    final_ddbb.insert(new_colum_index, sample, 0) #add a new column with defauls values = 0
                    
                    #Check if position exist
                    ########################
                    for position in new_sample['POS'].unique(): #extract first column in file
                        
                        if position not in final_ddbb["Position"].values:
                            positions_added.append(int(position)) #Count new positions for stats
                            
                            new_row = len(final_ddbb.index)
                            final_ddbb.loc[new_row,'Position'] = int(position)
                            final_ddbb.loc[new_row,'Samples'] = sample
                            final_ddbb.loc[new_row,'N'] = int(1)
                            final_ddbb.loc[new_row,sample] = str(1)
                        else:
                            positions_shared.append(int(position)) #Count shared positions for stats
                            
                            #Check whether the column matches the value and retrieve the first position [0]
                            #of the object index generated
                            index_position = final_ddbb.index[final_ddbb["Position"] == int(position)][0]
                            #Add sample to corresponding cell [position, samples]
                            number_samples_with_position = final_ddbb.loc[index_position,'N']
                            names_samples_with_position = final_ddbb.loc[index_position,'Samples']
                            new_names_samples = names_samples_with_position + "," + sample
                            #Sum 1 to the numbes of samples containing the position
                            final_ddbb.loc[index_position,'N'] = number_samples_with_position + 1
                            final_ddbb.loc[index_position,'Samples'] = new_names_samples
                            final_ddbb.loc[index_position,sample] = str(1) #Add "1" in cell with correct position vs sample (indicate present)

                    print("\nSAMPLE:\t%s\nTOTAL Variants:\t%s\nShared Variants:\t%s\nNew Variants:\t%s\n"
                    % (sample, len(new_sample.index), len(positions_shared), len(positions_added)))
                else:
                    print(YELLOW + "The sample " + sample + " ALREADY exist" + END_FORMATTING)

    #final_ddbb = final_ddbb.fillna(0).sort_values("Position") 
    final_ddbb["Position"] = final_ddbb["Position"].astype(int) #TO REMOVE when nucleotides are added
    final_ddbb['N'] = final_ddbb['N'].astype(int)
    #final_ddbb = final_ddbb.reset_index(drop=True)

    print("Final database now contains %s rows and %s columns" % final_ddbb.shape)
    if recalibrate == False:
        output_filename = output_filename + ".tsv"
        final_ddbb.to_csv(output_filename, sep='\t', index=False)
    else:
        recalibrate = os.path.abspath(recalibrate)
        if os.path.exists(recalibrate):
            recalibrate_params = extract_recalibrate_params(recalibrate)
            print("\n" + MAGENTA + "Recalibration selected" + END_FORMATTING)
            print(output_filename)
            output_filename = output_filename + ".revised.tsv"

            final_ddbb_revised = recalibrate_ddbb_vcf(final_ddbb, recalibrate_params[0], recalibrate_params[1], recalibrate_params[2])
            """
            if args.reference and args.reference != False:
                final_ddbb_revised = recalibrate_ddbb_vcf(final_ddbb, recalibrate_params[0], recalibrate_params[1], args.reference)
                
            else:
                final_ddbb_revised = recalibrate_ddbb_vcf(final_ddbb, recalibrate_params[0], recalibrate_params[1], recalibrate_params[2])
            """
            final_ddbb_revised.to_csv(output_filename, sep='\t', index=False)
        else:
            print("The directory supplied for recalculation does not exixt")
            sys.exit(1)
    print(output_filename)

    #Create small report with basic count
    #####################################
            
    print("\n" + GREEN + "Position check Finished" + END_FORMATTING)
    print(GREEN + "Added " + str(new_samples) + " samples out of " + str(all_samples) + END_FORMATTING + "\n")


    ###########################COMPARE FUNCTIONS#####################################################################

def compare_snp_columns(sample1, sample2, df):
    jaccard_similarity = accuracy_score(df[sample1], df[sample2]) #similarities between colums
    hamming_similarity = 1 - jaccard_similarity #disagreements between colums
    snp_distance = int(hamming_similarity * (len(df.index)+1))
    return snp_distance

def snp_distance_pairwise(dataframe, output_file):
    if os.path.exists(output_file):
        os.remove(output_file)
    with open(output_file, "a") as f:
        for sample1 in dataframe.iloc[:,3:].columns: #remove first 3 colums
            for sample2 in dataframe.iloc[:,3:].columns:
                if sample1 != sample2:
                    snp_distance = compare_snp_columns(sample1, sample2, dataframe)
                    line_distance = "%s\t%s\t%s\n" % (sample1, sample2, snp_distance)
                    f.write(line_distance)

def snp_distance_matrix(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    snp_distance_df = pd.DataFrame(hamming_distance * len(dataframe_only_samples.index), index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    snp_distance_df = snp_distance_df.astype(int)
    snp_distance_df.to_csv(output_file, sep='\t', index=True)

def hamming_distance_matrix(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    hamming_distance_df = pd.DataFrame(hamming_distance, index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    hamming_distance_df.to_csv(output_file, sep='\t', index=True)

def clustermap_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    sns.clustermap(dataframe_only_samples, annot=False, cmap="YlGnBu", figsize=(13, 13))
    plt.savefig(output_file, format="png")

def dendogram_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average') #method='single'

    plt.rcParams['lines.linewidth'] = 8 #Dendrogram line with
    plt.rcParams['xtick.major.size'] = 10 #Only affect to tick (line) size
    plt.rcParams.update({'font.size': 30}) #Increase x tick label size
    #plt.tick_params(labelsize=30)
    plt.figure(figsize=(30, 50))
    plt.ylabel('samples', fontsize=30)
    plt.xlabel('snp distance', fontsize=30)

    shc.dendrogram(Z, labels=labelList, orientation='left', distance_sort='descending', show_leaf_counts=True, color_threshold=10, leaf_font_size=20)

    
    plt.savefig(output_file, format="png")

# Convert dendrogram to Newick
def linkage_to_newick(dataframe, output_file):
    """
    Thanks to https://github.com/biocore/scikit-bio/issues/1579
    Input :  Z = linkage matrix, labels = leaf labels
    Output:  Newick formatted tree string
    """
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average')

    tree = shc.to_tree(Z, False)
    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            #print("%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick))
            return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = f"):{(parentdist - node.dist)/2}{newick}"
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = buildNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            #print(newick)
            return newick

    with open(output_file, 'w') as f:
        f.write(buildNewick(tree, "", tree.dist, labelList))
    return buildNewick(tree, "", tree.dist, labelList)

def matrix_to_rdf(snp_matrix, output_name):
    #snp_matrix = import_to_pandas(tsv_matrix, header=True)
    #tsv_matrix = os.path.abspath(tsv_matrix)
    #output_name = ".".join(tsv_matrix.split(".")[:-1]) + ".rdf"
    #output_name = output_name + ".rdf"

    max_samples = max(snp_matrix.N.tolist())
    snp_matrix = snp_matrix[snp_matrix.N < max_samples]
    
    with open(output_name, 'w+') as fout:
        snp_number = snp_matrix.shape[0]
        first_line = "  ;1.0\n"
        #print(first_line)
        fout.write(first_line)
        snp_list = snp_matrix.Position.astype(int).tolist()
        snp_list = " ;".join([str(x) for x in snp_list]) + " ;\n"
        #print(snp_list)
        fout.write(snp_list)
        third_line = ("10;" * snp_number) + "\n"
        #print(third_line)
        fout.write(third_line)
        transposed_snp_matrix = snp_matrix.T
        for index, row in transposed_snp_matrix.iloc[3:,:].iterrows():
            sample_header = ">"+ index+";1;;;;;;;\n"
            #print(sample_header)
            fout.write(sample_header)
            snp_row = "".join([str(x) for x in row.tolist()]) + "\n"
            #print(snp_row)
            fout.write(snp_row)
        ref_header = ">REF;1;;;;;;;\n"
        #print(ref_header)
        fout.write(ref_header)
        ref_snp = "0" * snp_number
        #print(ref_snp)
        fout.write(ref_snp)

def matrix_to_common(snp_matrix, output_name):
    #snp_matrix = import_to_pandas(tsv_matrix, header=True)
    #tsv_matrix = os.path.abspath(tsv_matrix)
    #output_name = ".".join(tsv_matrix.split(".")[:-1]) + ".rdf"
    #output_name = output_name + ".rdf"

    max_samples = max(snp_matrix.N.tolist())
    total_samples = len(snp_matrix.columns[3:])
    if max_samples == total_samples:
        with open(output_name, 'w+') as fout: 
            common_snps = snp_matrix['Position'][snp_matrix.N == max_samples].astype(int).astype(str).tolist()
            line = "\n".join(common_snps)
            fout.write("Position\n")
            fout.write(line)
    else:
        print("No common SNPs were found")

def ddtb_compare(final_database):

    database_file = os.path.abspath(final_database)
    check_file_exists(database_file)
    presence_ddbb = import_to_pandas(database_file, header=True)

    output_path = database_file.split(".")[0]

    print("Output path is: " + output_path)


    print(BLUE + BOLD + "Comparing all samples in " + database_file + END_FORMATTING)
    prior_pairwise = datetime.datetime.now()

    #Calculate pairwise snp distance for all and save file
    print(CYAN + "Pairwise distance" + END_FORMATTING)
    pairwise_file = output_path + ".snp.pairwise.tsv"
    snp_distance_pairwise(presence_ddbb, pairwise_file)
    after_pairwise = datetime.datetime.now()
    print("Done with pairwise in: %s" % (after_pairwise - prior_pairwise))

    #Calculate snp distance for all and save file
    print(CYAN + "SNP distance" + END_FORMATTING)
    snp_dist_file = output_path + ".snp.tsv"
    snp_distance_matrix(presence_ddbb, snp_dist_file)

    #Calculate hamming distance for all and save file
    print(CYAN + "Hamming distance" + END_FORMATTING)
    hmm_dist_file = output_path + ".hamming.tsv"
    hamming_distance_matrix(presence_ddbb, hmm_dist_file)
    """
    #Represent pairwise snp distance for all and save file
    print(CYAN + "Drawing distance" + END_FORMATTING)
    prior_represent = datetime.datetime.now()
    png_dist_file = output_path + ".snp.distance.png"
    #clustermap_dataframe(presence_ddbb, png_dist_file)
    after_represent = datetime.datetime.now()
    print("Done with distance drawing in: %s" % (after_represent - prior_represent))
    """
    #Represent dendrogram snp distance for all and save file
    print(CYAN + "Drawing dendrogram" + END_FORMATTING)
    png_dend_file = output_path + ".snp.dendrogram.png"
    dendogram_dataframe(presence_ddbb, png_dend_file)

    #Output a Newick file distance for all and save file
    print(CYAN + "Newick dendrogram" + END_FORMATTING)
    newick_file = output_path + ".nwk"
    linkage_to_newick(presence_ddbb, newick_file)

    #Output a binary snp matrix distance in rdf format
    print(CYAN + "rdf format" + END_FORMATTING)
    rdf_file = output_path + ".rdf"
    matrix_to_rdf(presence_ddbb, rdf_file)

    #Output a list of all common snps in group compared
    print(CYAN + "Common SNPs" + END_FORMATTING)
    common_file = output_path + ".common.txt"
    matrix_to_common(presence_ddbb, common_file)

    


if __name__ == '__main__':
    print("#################### COMPARE SNPS #########################")

    args = get_arguments()
    print(args)

    if args.recalibrate == False:
        compare_snp_matrix = args.output + ".tsv"
    else:
        compare_snp_matrix = args.output + ".revised.tsv"

    ddtb_add(args.input_dir, args.output, sample_filter=args.sample_list, recalibrate=args.recalibrate)
    ddtb_compare(compare_snp_matrix)