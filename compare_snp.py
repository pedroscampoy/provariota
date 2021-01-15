#!/usr/bin/env python

import os
import re
import logging
import pandas as pd
import numpy as np
import argparse
import sys
import subprocess
from sklearn.metrics import pairwise_distances, accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as ssd #pdist

logger = logging.getLogger()


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
    
    parser.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=False, help='REQUIRED.Input directory containing all vcf files')
    parser.add_argument('-s', '--sample_list', default=False, required=False, help='File with sample names to analyse instead of all samples')
    parser.add_argument('-d', '--distance', default=0, required=False, help='Minimun distance to cluster groups after comparison')
    parser.add_argument('-c', '--only-compare', dest="only_compare", required=False, default=False, help='Add already calculated snp binary matrix')
    parser.add_argument('-r', '--recalibrate', required= False, type=str, default=False, help='Coverage folder')
    parser.add_argument('-b', '--bam_folder', required= False, type=str, default=False, help='Bam folder')
    parser.add_argument('-R', '--reference', required= False, type=str, default=False, help='Reference fasta file used in original variant calling')

    parser.add_argument('-o', '--output', type=str, required=True, help='Name of all the output files, might include path')

    arguments = parser.parse_args()

    return arguments

def check_file_exists(file_name):
    """
        Check file exist and is not 0 Kb, if not program exit.
    """
    file_info = os.stat(file_name) #Retrieve the file info to check if has size > 0

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        logger.info(RED + BOLD + "File: %s not found or empty\n" % file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)

def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def blank_database():
    new_pandas_ddtb = pd.DataFrame(columns=['Position','N', 'Samples'])
    return new_pandas_ddtb

####COMPARE SNP 2.0
def import_tsv_variants(tsv_file,  min_total_depth=4, min_alt_dp=4, only_snp=True):
    base_file = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    sample = base_file.split(".")[0]

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")

    df = df[((df.TOTAL_DP >= min_total_depth) &
                    (df.ALT_DP >= min_alt_dp))]

    df = df[['REGION','POS', 'REF', 'ALT', 'ALT_FREQ', 'TYPE']]
    df = df.rename(columns={'ALT_FREQ' : sample})
    if only_snp == True:
        df = df[df.TYPE == 'snp']
        df = df.drop(['TYPE'], axis=1)
        return df
    else:
        return df

def extract_lowfreq(tsv_file,  min_total_depth=4, min_alt_dp=4, min_freq_include=0.7, only_snp=True):
    base_file = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    sample = base_file.split(".")[0]

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")

    df = df[(df.ALT_DP < min_alt_dp) &
            (df.ALT_FREQ >= min_freq_include)]

    df = df[['REGION','POS', 'REF', 'ALT', 'ALT_FREQ', 'TYPE']]
    df['ALT_FREQ'] = '?'
    df = df.rename(columns={'ALT_FREQ' : sample})
    if only_snp == True:
        df = df[df.TYPE == 'snp']
        df = df.drop(['TYPE'], axis=1)
        return df
    else:
        return df

def extract_uncovered(cov_file, min_total_depth=4):
    base_file = os.path.basename(cov_file)
    input_file = os.path.abspath(cov_file)
    sample = base_file.split(".")[0]

    df = pd.read_csv(input_file, sep="\t", header=None)
    df.columns = ['REGION', 'POS', sample]
    df = df[df[sample] == 0]
    df = df.replace(0,'!')
    return df

def ddbb_create_intermediate(variant_dir, coverage_dir, min_freq_discard=0.1, min_alt_dp=10):
    df = pd.DataFrame(columns=['REGION','POS', 'REF', 'ALT'])
    #Merge all raw
    for root, _, files in os.walk(variant_dir):
        if root == variant_dir:
            for name in files:
                if name.endswith('.tsv'):
                    logger.debug("Adding: " + name)
                    filename = os.path.join(root, name)
                    dfv = import_tsv_variants(filename)
                    df = df.merge(dfv, how='outer')
    #Rounf frequencies
    df = df.round(2)
    #Remove <= 0.1 (parameter in function)
    handle_lowfreq = lambda x: None if x <= min_freq_discard else x # IF HANDLE HETEROZYGOUS CHANGE THIS 0 for X or 0.5
    df.iloc[:,4:] = df.iloc[:,4:].applymap(handle_lowfreq)
    #Drop all NaN rows
    df['AllNaN'] = df.apply(lambda x: x[4:].isnull().values.all(), axis=1)
    df = df[df.AllNaN == False]
    df = df.drop(['AllNaN'], axis=1).reset_index(drop=True)

    #Include poorly covered
    for root, _, files in os.walk(variant_dir):
        if root == variant_dir:
            for name in files:
                if name.endswith('.tsv'):
                    filename = os.path.join(root, name)
                    sample = name.split('.')[0]
                    logger.debug("Adding lowfreqs: " + sample)
                    dfl = extract_lowfreq(filename, min_total_depth=4, min_alt_dp=min_alt_dp, only_snp=True)
                    df[sample].update(df[['REGION', 'POS', 'REF', 'ALT']].merge(dfl, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])

    #Include uncovered
    samples_coverage = df.columns.tolist()[4:]
    for root, _, files in os.walk(coverage_dir):
        for name in files:
            if name.endswith('.cov'):
                filename = os.path.join(root, name)
                sample = name.split('.')[0]
                if sample in df.columns[4:]:
                    samples_coverage.remove(sample)
                    logger.debug("Adding uncovered: " + sample)
                    dfc = extract_uncovered(filename)
                    #df.update(df[['REGION', 'POS']].merge(dfc, on=['REGION', 'POS'], how='left'))
                    df[sample].update(df[['REGION', 'POS']].merge(dfc, on=['REGION', 'POS'], how='left')[sample])
                    #df.combine_first(df[['REGION', 'POS']].merge(dfc, how='left'))
    if len(samples_coverage) > 0:
        logger.info("WARNING: " + (',').join(samples_coverage) + " coverage file not found" )
    #Asign 0 to rest (Absent)
    df = df.fillna(0)

    #Determine N (will help in poorly covered determination)
    def estract_sample_count(row):
        count_list = [i not in ['!',0,'0'] for i in row[4:]]
        samples = np.array(df.columns[4:])
        #samples[np.array(count_list)] filter array with True False array
        return (sum(count_list), (',').join(samples[np.array(count_list)]))

    if 'N' in df.columns:
        df = df.drop(['N','Samples'], axis=1)
    if 'Position' in df.columns:
        df = df.drop('Position', axis=1)

    df[['N', 'Samples']] = df.apply(estract_sample_count, axis=1, result_type='expand')

    df['Position'] = df.apply(lambda x: ('|').join([x['REGION'],x['REF'],str(x['POS']),x['ALT']]), axis=1)

    df = df.drop(['REGION','REF','POS','ALT'], axis=1)

    df = df[['Position', 'N', 'Samples'] + [ col for col in df.columns if col not in ['Position', 'N', 'Samples']]]

    return df



################################ END COMPARE SNP 2.0

def import_VCF4_to_pandas(vcf_file, sep='\t'):
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            #logger.info(next_line)
            next_line = f.readline()
    
    if first_line.endswith('VCFv4.2'):
        
        #Use first line as header
        dataframe = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)
        sample = dataframe.columns[-1]
        dataframe.rename(columns={sample:'sample'}, inplace=True)
        dataframe['POS'] = dataframe['POS'].astype(int)
        
    else:
        logger.info("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe

def import_to_pandas(file_table, header=False, sep='\t'):
    if header == False:
        #exclude first line, exclusive for vcf outputted by PipelineTB
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        #Use first line as header
        dataframe = pd.read_csv(file_table, sep=sep, header=0)
    
    return dataframe

def import_tsv_pandas(vcf_file, sep='\t'):
    if check_file_exists(vcf_file):
        dataframe = pd.read_csv(vcf_file, sep=sep, header=0)
        dataframe['POS'] = dataframe['POS'].astype(int)
    else:
        logger.info("This vcf file is empty or does not exist")
        sys.exit(1)
           
    return dataframe

def ddtb_add(input_folder, output_filename, sample_filter=False, vcf_suffix=".tsv" ):
    directory = os.path.abspath(input_folder)
    output_filename = os.path.abspath(output_filename)

    #Make sure output exist to force change name
    if os.path.isfile(output_filename):
        logger.info(YELLOW + "ERROR: " + BOLD + "output database EXIST, choose a different name or manually delete" + END_FORMATTING)
        sys.exit(1)

    final_ddbb = blank_database()
    sample_filter_list = []

    #Handle sample filter
    if sample_filter == False:
        sample_filter_list = [x.split(".")[0] for x in os.listdir(directory) if x.endswith(vcf_suffix)]
    else:
        if os.path.isfile(sample_filter):
            with open(sample_filter, 'r') as f:
                for line in f:
                    sample_filter_list.append(line.strip())
        else:
            "Sample file don't exist"
            sys.exit(1)
    
    if len(sample_filter_list) < 1:
        logger.info("prease provide 2 or more samples")
        sys.exit(1)

    #logger.info("Previous final database contains %s rows and %s columns\n" % final_ddbb.shape)
    logger.info("The directory selected is: %s" % directory)
    

    all_samples = 0
    new_samples = 0
    for filename in os.listdir(directory):
        if not filename.startswith('.') and filename.endswith(vcf_suffix):
            
            all_samples = all_samples + 1
            positions_shared = []
            positions_added = []
            
            sample = filename.split(".")[0] #Manage sample name

            if sample in sample_filter_list:
                logger.info("\nThe file is: %s" % filename)

                file = os.path.join(directory, filename) #Whole file path
                check_file_exists(file) #Manage file[s]. Check if file exist and is greater than 0

                new_sample = import_tsv_pandas(file) #Import files in annotated tsv format

                #Check if sample exist
                ######################
                if sample not in final_ddbb.columns.tolist():
                    logger.info("Adding new sample %s to %s" % (sample, os.path.basename(output_filename)))
                    new_samples = new_samples + 1
                    new_colum_index = len(final_ddbb.columns) #extract the number of columns to insert a new one
                    #final_ddbb[sample] = sample #adds a new column but fills all blanks with the value sample
                    final_ddbb.insert(new_colum_index, sample, 0) #add a new column with defauls values = 0
                    
                    #Check if position exist
                    ########################
                    for _, row in new_sample.iterrows():

                        position = ('|').join([row['REGION'],row['REF'],str(row['POS']),row['ALT']])
                        
                        if position not in final_ddbb["Position"].values:
                            positions_added.append(position) #Count new positions for stats
                            
                            new_row = len(final_ddbb.index)
                            final_ddbb.loc[new_row,'Position'] = position
                            final_ddbb.loc[new_row,'Samples'] = sample
                            final_ddbb.loc[new_row,'N'] = int(1)
                            final_ddbb.loc[new_row,sample] = str(1)
                        else:
                            positions_shared.append(position) #Count shared positions for stats
                            
                            #Check whether the column matches the value and retrieve the first position [0]
                            #of the object index generated
                            index_position = final_ddbb.index[final_ddbb["Position"] == position][0]
                            #Add sample to corresponding cell [position, samples]
                            number_samples_with_position = final_ddbb.loc[index_position,'N']
                            names_samples_with_position = final_ddbb.loc[index_position,'Samples']
                            new_names_samples = names_samples_with_position + "," + sample
                            #Sum 1 to the numbes of samples containing the position
                            final_ddbb.loc[index_position,'N'] = number_samples_with_position + 1
                            final_ddbb.loc[index_position,'Samples'] = new_names_samples
                            final_ddbb.loc[index_position,sample] = str(1) #Add "1" in cell with correct position vs sample (indicate present)

                    logger.info("\nSAMPLE:\t%s\nTOTAL Variants:\t%s\nShared Variants:\t%s\nNew Variants:\t%s\n"
                    % (sample, len(new_sample.index), len(positions_shared), len(positions_added)))
                else:
                    logger.info(YELLOW + "The sample " + sample + " ALREADY exist" + END_FORMATTING)

    final_ddbb = final_ddbb.fillna(0).sort_values("Position") 
    #final_ddbb["Position"] = final_ddbb["Position"].astype(int) #TO REMOVE when nucleotides are added
    final_ddbb['N'] = final_ddbb['N'].astype(int)
    #final_ddbb = final_ddbb.reset_index(drop=True)

    logger.info("Final database now contains %s rows and %s columns" % final_ddbb.shape)
    
    output_filename = output_filename + ".tsv"
    final_ddbb.to_csv(output_filename, sep='\t', index=False)
    logger.info(output_filename)

    #Create small report with basic count
    #####################################
            
    logger.info("\n" + GREEN + "Position check Finished" + END_FORMATTING)
    logger.info(GREEN + "Added " + str(new_samples) + " samples out of " + str(all_samples) + END_FORMATTING + "\n")

###########################MATRIX TO CLUSTER FUNCTIONS###########################################################

def pairwise_to_cluster(pw,threshold = 0):
    groups = {}
    columns = pw.columns.tolist()
    sorted_df = pw[(pw[columns[0]] != pw[columns[1]]) & (pw[columns[2]] <= threshold)].sort_values(by=[columns[2]])
    
    def rename_dict_clusters(cluster_dict):
        reordered_dict = {}
        for i, k in enumerate(list(cluster_dict)):
            reordered_dict[i] = cluster_dict[k]
        return reordered_dict
    
    def regroup_clusters(list_keys, groups_dict, both_samples_list):
        #sum previous clusters
        list_keys.sort()
        new_cluster = sum([groups_dict[key] for key in list_keys], [])
        #add new cluster
        cluster_asign = list(set(new_cluster + both_samples_list))
        #Remove duped cluster
        first_cluster = list_keys[0]
        groups_dict[first_cluster] = cluster_asign
        rest_cluster = list_keys[1:]
        for key in rest_cluster:
            del groups_dict[key]
        groups_dict = rename_dict_clusters(groups_dict)
        return groups_dict
        
    for _, row in sorted_df.iterrows():
        group_number = len(groups)
        sample_1 = str(row[0])
        sample_2 = str(row[1])
        both_samples_list = row[0:2].tolist()
                
        if group_number == 0:
            groups[group_number] = both_samples_list
        
        all_samples_dict = sum(groups.values(), [])
                
        if sample_1 in all_samples_dict or sample_2 in all_samples_dict:
            #extract cluster which have the new samples
            key_with_sample = {key for (key,value) in groups.items() if (sample_1 in value or sample_2 in value)}
            
            cluster_with_sample = list(key_with_sample)
            cluster_with_sample_name = cluster_with_sample[0]
            number_of_shared_clusters = len(key_with_sample)
            if number_of_shared_clusters > 1:
                groups = regroup_clusters(cluster_with_sample, groups, both_samples_list)
            else:
                groups[cluster_with_sample_name] = list(set(groups[cluster_with_sample_name] + both_samples_list))
        else:
            groups[group_number] = both_samples_list
            
    for _, row in pw[(pw[pw.columns[0]] != pw[pw.columns[1]]) & (pw[pw.columns[2]] > threshold)].iterrows():
        sample_1 = str(row[0])
        sample_2 = str(row[1])
        all_samples_dict = sum(groups.values(), [])
        if sample_1 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_1]
        
        if sample_2 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_2]
            
    cluster_df = pd.DataFrame(groups.values(),index=list(groups))
    
    cluster_df_return = cluster_df.stack().droplevel(1).reset_index().rename(columns={'index': 'group', 0: 'id'})
            
    return cluster_df_return

def calculate_N(row):
    return len(row.samples)

def calculate_mean_distance(row, df):
    if row.N > 1:
        list_sample = row.samples
        list_sample = [str(x) for x in list_sample]
        list_sample = [x.split(".")[0] for x in list_sample]
        dataframe = df.loc[list_sample,list_sample]
        stacked_df = dataframe.stack()
        mean_distance = stacked_df.mean(skipna = True)
        min_distance = stacked_df.min(skipna = True)
        max_distance = stacked_df.max(skipna = True)
        return round(mean_distance, 2), min_distance, max_distance
    else:
        return 'NaN', 'NaN', 'NaN'

def matrix_to_cluster(pairwise_file, matrix_file, distance=0):
    output_dir = ('/').join(pairwise_file.split('/')[0:-1])

    logger.info('Reading Matrix')
    dfdist = pd.read_csv(matrix_file, index_col=0, sep='\t', )
    dfdist.columns = dfdist.columns.astype(str)
    dfdist.index = dfdist.index.astype(str)
    logger.info('Reading Pairwise')
    pairwise = pd.read_csv(pairwise_file, sep="\t", names=['sample_1', 'sample_2', 'dist'])
    logger.info('Creating Clusters')
    clusters = pairwise_to_cluster(pairwise,threshold=distance)

    cluster_summary = clusters.groupby('group')['id'].apply(list).reset_index(name='samples')
    cluster_summary['N'] = cluster_summary.apply(calculate_N, axis=1)
    cluster_summary = cluster_summary.sort_values(by=['N'], ascending=False)

    logger.info('Reseting group number by length')
    sorted_index = cluster_summary.index.to_list()
    sorted_index.sort()
    sorted_index = [x + 1 for x in sorted_index]
    cluster_summary['group'] = sorted_index
    cluster_summary = cluster_summary.sort_values(by=['N'], ascending=False)

    cluster_summary[['mean', 'min', 'max']] = cluster_summary.apply(lambda x: calculate_mean_distance(x, dfdist), axis=1, result_type="expand")

    final_cluster = cluster_summary[["group", "samples"]].explode("samples").reset_index(drop=True)
    final_cluster = final_cluster.sort_values(by=['group'], ascending=True)

    final_cluster_file = os.path.join(output_dir, "group_table_" + str(distance) + ".tsv")
    cluster_summary_file = os.path.join(output_dir, "group_summary_" + str(distance) + ".tsv")

    cluster_summary.to_csv(cluster_summary_file, sep='\t', index=False)
    final_cluster.to_csv(final_cluster_file, sep='\t', index=False)


def revised_df(df, out_dir=False, min_freq_include=0.8, min_threshold_discard_uncov_sample=0.4, min_threshold_discard_uncov_pos=0.4, min_threshold_discard_htz_sample=0.4, min_threshold_discard_htz_pos=0.4, min_threshold_discard_all_pos=0.6, min_threshold_discard_all_sample=0.6, remove_faulty=True, drop_samples=True, drop_positions=True):
    if remove_faulty == True:

        uncovered_positions = df.iloc[:,3:].apply(lambda x:  sum([i in ['!','?'] for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=1)
        heterozygous_positions = df.iloc[:,3:].apply(lambda x: sum([(i not in ['!','?',0,1, '0', '1']) and (float(i) < min_freq_include) and (float(i) > 0.1) for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=1)
        report_position = pd.DataFrame({'Position': df.Position, 'uncov_fract': uncovered_positions, 'htz_frac': heterozygous_positions, 'faulty_frac': uncovered_positions + heterozygous_positions})
        faulty_positions = report_position['Position'][(report_position.uncov_fract >= min_threshold_discard_uncov_pos) | (report_position.htz_frac >= min_threshold_discard_htz_pos) | (report_position.faulty_frac >= min_threshold_discard_all_pos)].tolist()


        uncovered_samples = df.iloc[:,3:].apply(lambda x: sum([i in ['!','?'] for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=0)
        heterozygous_samples = df.iloc[:,3:].apply(lambda x: sum([(i not in ['!','?',0,1, '0', '1']) and (float(i) < min_freq_include) and (float(i) > 0.1) for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=0)
        report_samples = pd.DataFrame({'sample': df.iloc[:,3:].columns, 'uncov_fract': uncovered_samples, 'htz_frac': heterozygous_samples, 'faulty_frac': uncovered_samples + heterozygous_samples})
        faulty_samples = report_samples['sample'][(report_samples.uncov_fract >= min_threshold_discard_uncov_sample) | (report_samples.htz_frac >= min_threshold_discard_htz_sample) | (report_samples.faulty_frac >= min_threshold_discard_all_sample)].tolist()

        if out_dir != False:
            out_dir = os.path.abspath(out_dir)
            report_samples_file = os.path.join(out_dir, 'report_samples.tsv')
            report_faulty_samples_file = os.path.join(out_dir, 'faulty_samples.tsv')
            report_positions_file = os.path.join(out_dir, 'report_positions.tsv')
            report_faulty_positions_file = os.path.join(out_dir, 'faulty_positions.tsv')
            intermediate_cleaned_file = os.path.join(out_dir, 'intermediate.highfreq.tsv')
            report_position.to_csv(report_positions_file, sep="\t", index=False)
            report_samples.to_csv(report_samples_file, sep="\t", index=False)
            with open(report_faulty_samples_file, 'w+') as f:
                f.write(('\n').join(faulty_samples))
            with open(report_faulty_positions_file, 'w+') as f2:
                f2.write(('\n').join(faulty_positions))

        if drop_positions == True:
            df = df[~df.Position.isin(faulty_positions)]
        if drop_samples == True:
            df = df.drop(faulty_samples, axis=1)

        logger.info('FAULTY POSITIONS:\n{}\n\nFAULTY SAMPLES:\n{}'.format(("\n").join(faulty_positions), ("\n").join(faulty_samples)))

    #Uncovered to 0
    #Number of valid to remove o valid and replace lowfreq
    df['valid'] = df.apply(lambda x: sum([i != '?' and i != '!' and float(i) >= min_freq_include for i in x[3:]]), axis=1)
    df = df[df.valid >= 1]
    df = df.drop('valid', axis=1)

    if out_dir != False:
        df.to_csv(intermediate_cleaned_file, sep="\t", index=False)

    df = df.replace('!', 0)

    df = df.replace('?', 1)

    df.iloc[:,3:] = df.iloc[:,3:].astype(float)

    #Replace Htz to 0
    #f = lambda x: 1 if x >= min_freq_include else 0 # IF HANDLE HETEROZYGOUS CHANGE THIS 0 for X or 0.5
    #df.iloc[:,3:] = df.iloc[:,3:].applymap(f)

    #Replace Htz with 1
    f = lambda x: 1 if x > 0.5 else 0 # IF HANDLE HETEROZYGOUS CHANGE THIS 0 for X or 0.5
    df.iloc[:,3:] = df.iloc[:,3:].applymap(f)

    df.N = df.apply(lambda x: sum(x[3:]), axis=1)

    def extract_sample_name(row):
        count_list = [i not in ['!',0,'0'] for i in row[3:]]
        samples = np.array(df.columns[3:])
        #samples[np.array(count_list)] filter array with True False array
        return ((',').join(samples[np.array(count_list)]))


    df['Samples'] = df.apply(extract_sample_name, axis=1)

    #Remove positions with 0 samples after htz
    df = df[df.N > 0]
    
    return df

def recheck_variant_mpileup_intermediate(reference_id, position, alt_snp, sample, previous_binary, bam_folder, min_cov_low_freq = 10):

    """
    http://www.htslib.org/doc/samtools-mpileup.html
    www.biostars.org/p/254287
    In the pileup format (without -u or -g), each line represents a genomic position, consisting of chromosome name, 1-based coordinate,
    reference base, the number of reads covering the site, read bases, base qualities and alignment mapping qualities. Information on match,
    mismatch, indel, strand, mapping quality and start and end of a read are all encoded at the read base column. At this column, a dot stands
    for a match to the reference base on the forward strand, a comma for a match on the reverse strand, a '>' or '<' for a reference skip,
    ACGTN for a mismatch on the forward strand and acgtn for a mismatch on the reverse strand. A pattern +[0-9]+[ACGTNacgtn]+ indicates there
    is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the
    pattern, followed by the inserted sequence. Similarly, a pattern -[0-9]+[ACGTNacgtn]+ represents a deletion from the reference. The deleted
    bases will be presented as * in the following lines. Also at the read base column, a symbol ^ marks the start of a read. The ASCII of the
    character following ^ minus 33 gives the mapping quality. A symbol $ marks the end of a read segment
    """
    if previous_binary != 0 and previous_binary != '0':
        #logger.info('NON0: {} in sample {} is not 0: {}'.format(position, sample, previous_binary))
        return previous_binary
    else:
        #previous_binary = int(previous_binary)
        position = int(position)

        #Identify correct bam
        for root, _, files in os.walk(bam_folder):
            for name in files:
                filename = os.path.join(root, name)
                sample_file = name.split('.')[0]
                if name.startswith(sample) and sample_file == sample and name.endswith(".bam"):
                    bam_file = filename
        #format position for mpileup execution (NC_000962.3:632455-632455)
        position_format = reference_id + ":" + str(position) + "-" + str(position)
        
        #Execute command and retrieve output
        cmd = ["samtools", "mpileup", "-aa", "-r", position_format, bam_file]
        text_mpileup = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
        split_mpileup = text_mpileup.stdout.split()
        #Extract 5th column to find variants
        mpileup_reference = split_mpileup[0]
        mpileup_position = int(split_mpileup[1])
        mpileup_depth = int(split_mpileup[3])
        mpileup_variants = split_mpileup[4]
        variant_list = list(mpileup_variants)
        variants_to_account = ['A', 'T', 'C', 'G', '.', ',']
        variant_upper_list = [x.upper() for x in variant_list]
        variant_upper_list = [x for x in variant_upper_list if x in variants_to_account]

        if len(variant_upper_list) == 0:
            most_counted_variant = "*"
            mpileup_depth = 0

        if reference_id != mpileup_reference:
            logger.info('ERROR: References are different')
            sys.exit(1)
        else:
            if  mpileup_depth == 0:
                logger.debug('WARNING: SAMPLE: {} has 0 depth in position {}. BEFORE: {}'.format(sample, position, previous_binary))
                return '!'
            elif mpileup_depth > 0:
                most_counted_variant = max(set(variant_upper_list), key = variant_upper_list.count)
                count_all_variants = {x:variant_upper_list.count(x) for x in variant_upper_list}
                freq_most_frequent = count_all_variants[most_counted_variant]/len(variant_upper_list)
                freq_most_frequent = round(freq_most_frequent,2)

                if (most_counted_variant == alt_snp) and (freq_most_frequent >= 0.1) and (mpileup_depth > min_cov_low_freq):
                    logger.debug('WARNING: SAMPLE: {} has heterozygous position at {} with frequency {}. BEFORE: {}: {}'.format(sample, position, freq_most_frequent, previous_binary, text_mpileup.stdout))
                    return freq_most_frequent
                elif (most_counted_variant == alt_snp) and (freq_most_frequent >= 0.6) and (mpileup_depth <= min_cov_low_freq):
                    logger.debug('WARNING: SAMPLE: {} has lowcov position at {} with frequency {}. BEFORE: {}: {}'.format(sample, position, freq_most_frequent, previous_binary, text_mpileup.stdout))
                    return '?'
                elif (most_counted_variant == alt_snp):
                    return freq_most_frequent
                elif (alt_snp in variant_upper_list):
                    very_low_freq = round(count_all_variants[alt_snp]/len(variant_upper_list), 2)
                    logger.debug('SAMPLE: {} has VERY lowcov position at {} with frequency {}. BEFORE: {}: {}'.format(sample, position, very_low_freq, previous_binary, text_mpileup.stdout))
                    return very_low_freq
                else:
                    logger.debug('ELSE {}, most_counted_variant {}, freq_most_frequent {}'.format(text_mpileup.stdout, most_counted_variant, freq_most_frequent))
                    return 0

def recalibrate_ddbb_vcf_intermediate(snp_matrix_ddbb_file, bam_folder, min_cov_low_freq = 10):
    
    df_matrix = pd.read_csv(snp_matrix_ddbb_file, sep="\t")
    
    sample_list_matrix = df_matrix.columns[3:]
    n_samples = len(sample_list_matrix)
    
    #Iterate over non unanimous positions 
    for index, data_row in df_matrix[df_matrix.N < n_samples].iloc[:,3:].iterrows():
        #Extract its position
        whole_position = df_matrix.loc[index,"Position"]
        row_reference = whole_position.split('|')[0]
        row_position = int(whole_position.split('|')[2])
        row_alt_snp = whole_position.split('|')[3]

        #Use enumerate to retrieve column index (column ondex + 3)
        #find positions with frequency >80% in mpileup execution
        #Returns ! for coverage 0
        new_presence_row = [recheck_variant_mpileup_intermediate(row_reference, row_position, row_alt_snp, df_matrix.columns[n + 3], x, bam_folder, min_cov_low_freq = 10) for n,x in enumerate(data_row)]
        
        df_matrix.iloc[index, 3:] = new_presence_row
        #df_matrix.loc[index, 'N'] = sum([x == 1 for x in new_presence_row])

    def estract_sample_count(row):
        count_list = [i not in ['!',0,'0'] for i in row[3:]]
        samples = np.array(df_matrix.columns[3:])
        #samples[np.array(count_list)] filter array with True False array
        return (sum(count_list), (',').join(samples[np.array(count_list)]))


    df_matrix[['N', 'Samples']] = df_matrix.apply(estract_sample_count, axis=1, result_type='expand')

    return df_matrix


###########################COMPARE FUNCTIONS#####################################################################
#################################################################################################################

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

def snp_distance_matrix(dataframe, output_matrix, output_pairwise):
    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    snp_distance_df = pd.DataFrame(hamming_distance * len(dataframe_only_samples.index), index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    snp_distance_df = snp_distance_df.astype(int)
    pairwise = snp_distance_df.stack().reset_index(name='distance').rename(columns={'level_0': 'sample_1', 'level_1': 'sample_2'})
    snp_distance_df.to_csv(output_matrix, sep='\t', index=True)
    pairwise.to_csv(output_pairwise, sep='\t', header=False, index=False)

def hamming_distance_matrix(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    hamming_distance_df = pd.DataFrame(hamming_distance, index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    hamming_distance_df.to_csv(output_file, sep='\t', index=True)

def clustermap_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    sns.clustermap(dataframe_only_samples, annot=False, cmap="YlGnBu", figsize=(13, 13))
    plt.savefig(output_file, format="png")

def dendogram_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
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
    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average')

    tree = shc.to_tree(Z, False)
    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            #logger.info("%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick))
            return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = f"):{(parentdist - node.dist)/2}{newick}"
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = buildNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            #logger.info(newick)
            return newick

    with open(output_file, 'w') as f:
        f.write(buildNewick(tree, "", tree.dist, labelList))
    return buildNewick(tree, "", tree.dist, labelList)

def matrix_to_rdf(snp_matrix, output_name):
    #snp_matrix = import_to_pandas(tsv_matrix, header=True)
    #tsv_matrix = os.path.abspath(tsv_matrix)
    #output_name = ".".join(tsv_matrix.split(".")[:-1]) + ".rdf"
    #output_name = output_name + ".rdf"

    #max_samples = max(snp_matrix.N.tolist())
    #snp_matrix = snp_matrix[snp_matrix.N < max_samples]
    
    with open(output_name, 'w+') as fout:
        snp_number = snp_matrix.shape[0]
        first_line = "  ;1.0\n"
        #logger.info(first_line)
        fout.write(first_line)
        snp_list = snp_matrix.Position.tolist()
        snp_list = [x.split('|')[2] for x in snp_list]
        snp_list = " ;".join([str(x) for x in snp_list]) + " ;\n"
        #logger.info(snp_list)
        fout.write(snp_list)
        third_line = ("10;" * snp_number) + "\n"
        #logger.info(third_line)
        fout.write(third_line)
        transposed_snp_matrix = snp_matrix.T
        for index, row in transposed_snp_matrix.iloc[3:,:].iterrows():
            sample_header = ">"+ index+";1;;;;;;;\n"
            #logger.info(sample_header)
            fout.write(sample_header)
            snp_row = "".join([str(x) for x in row.tolist()]) + "\n"
            #logger.info(snp_row)
            fout.write(snp_row)
        ref_header = ">REF;1;;;;;;;\n"
        #logger.info(ref_header)
        fout.write(ref_header)
        ref_snp = "0" * snp_number
        #logger.info(ref_snp)
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
            common_snps = snp_matrix['Position'][snp_matrix.N == max_samples].astype(str).tolist()
            line = "\n".join(common_snps)
            fout.write("Position\n")
            fout.write(line)
    else:
        logger.info("No common SNPs were found")

def ddtb_compare(final_database, distance=0):

    database_file = os.path.abspath(final_database)
    check_file_exists(database_file)
    presence_ddbb = import_to_pandas(database_file, header=True)

    output_path = database_file.split(".")[0]

    logger.info("Output path is: " + output_path)


    logger.info(BLUE + BOLD + "Comparing all samples in " + database_file + END_FORMATTING)
    #prior_pairwise = datetime.datetime.now()
    #Calculate pairwise snp distance for all and save file
    #logger.info(CYAN + "Pairwise distance" + END_FORMATTING)
    
    #snp_distance_pairwise(presence_ddbb, pairwise_file)
    #after_pairwise = datetime.datetime.now()
    #logger.info("Done with pairwise in: %s" % (after_pairwise - prior_pairwise))

    #Calculate snp distance for all and save file
    logger.info(CYAN + "SNP distance" + END_FORMATTING)
    snp_dist_file = output_path + ".snp.tsv"
    pairwise_file = output_path + ".snp.pairwise.tsv"
    snp_distance_matrix(presence_ddbb, snp_dist_file, pairwise_file)

    #Calculate hamming distance for all and save file
    logger.info(CYAN + "Hamming distance" + END_FORMATTING)
    hmm_dist_file = output_path + ".hamming.tsv"
    hamming_distance_matrix(presence_ddbb, hmm_dist_file)
    """
    #Represent pairwise snp distance for all and save file
    logger.info(CYAN + "Drawing distance" + END_FORMATTING)
    prior_represent = datetime.datetime.now()
    png_dist_file = output_path + ".snp.distance.png"
    #clustermap_dataframe(presence_ddbb, png_dist_file)
    after_represent = datetime.datetime.now()
    logger.info("Done with distance drawing in: %s" % (after_represent - prior_represent))
    """
    #Represent dendrogram snp distance for all and save file
    logger.info(CYAN + "Drawing dendrogram" + END_FORMATTING)
    png_dend_file = output_path + ".snp.dendrogram.png"
    dendogram_dataframe(presence_ddbb, png_dend_file)

    #Output a Newick file distance for all and save file
    logger.info(CYAN + "Newick dendrogram" + END_FORMATTING)
    newick_file = output_path + ".nwk"
    linkage_to_newick(presence_ddbb, newick_file)

    #Output a binary snp matrix distance in rdf format
    logger.info(CYAN + "rdf format" + END_FORMATTING)
    rdf_file = output_path + ".rdf"
    matrix_to_rdf(presence_ddbb, rdf_file)

    #Output a list of all common snps in group compared
    logger.info(CYAN + "Common SNPs" + END_FORMATTING)
    common_file = output_path + ".common.txt"
    matrix_to_common(presence_ddbb, common_file)

    #Output files with group/cluster assigned to samples
    logger.info(CYAN + "Assigning clusters" + END_FORMATTING)
    matrix_to_cluster(pairwise_file, snp_dist_file, distance=distance)

    


if __name__ == '__main__':
    args = get_arguments()
    
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split('/')[-1]
    check_create_dir(output_dir)
    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = group_name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output_dir, 'Logs')
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


    logger.info("#################### COMPARE SNPS #########################")
    logger.info(args)

    group_compare = os.path.join(output_dir, group_name)
    compare_snp_matrix = group_compare + ".tsv"
    
    if args.only_compare == False:
        input_dir = os.path.abspath(args.input_dir)

        if args.recalibrate == False:
            ddtb_add(input_dir, group_compare, sample_filter=args.sample_list)
            ddtb_compare(compare_snp_matrix, distance=args.distance)
        else:
            coverage_dir = os.path.abspath(args.recalibrate)
            compare_snp_matrix_recal = group_compare + ".revised.final.tsv"
            compare_snp_matrix_recal_intermediate = group_compare + ".revised_intermediate.tsv"
            compare_snp_matrix_recal_mpileup = group_compare + ".revised_intermediate_mpileup.tsv"
            recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(input_dir, coverage_dir, min_freq_discard=0.1, min_alt_dp=10)
            recalibrated_snp_matrix_intermediate.to_csv(compare_snp_matrix_recal_intermediate, sep="\t", index=False)

            recalibrated_snp_matrix_mpileup = recalibrate_ddbb_vcf_intermediate(compare_snp_matrix_recal_intermediate, args.bam_folder, min_cov_low_freq = 10)
            recalibrated_snp_matrix_mpileup.to_csv(compare_snp_matrix_recal_mpileup, sep="\t", index=False)

            recalibrated_revised_df = revised_df(recalibrated_snp_matrix_mpileup, output_dir, min_freq_include=0.75, min_threshold_discard_uncov_sample=0.4, min_threshold_discard_uncov_pos=0.4, min_threshold_discard_htz_sample=0.6, min_threshold_discard_htz_pos=0.6, min_threshold_discard_all_pos=0.6, min_threshold_discard_all_sample=0.6, remove_faulty=True, drop_samples=True, drop_positions=True)
            recalibrated_revised_df.to_csv(compare_snp_matrix_recal, sep="\t", index=False)
            ddtb_compare(compare_snp_matrix_recal, distance=args.distance)
    else:
        compare_matrix = os.path.abspath(args.only_compare)
        ddtb_compare(compare_matrix, distance=args.distance)
