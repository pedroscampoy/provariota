#!/usr/bin/env python

import os
import sys
import logging
import argparse
import pandas as pd
import numpy as np
import re
import subprocess
from tabulate import tabulate
from misc import get_snpeff_path, check_create_dir
from vcf_process import bed_to_df,bed_to_df, add_bed_info, annotate_bed_s, obtain_output_dir, calculate_ALT_AD, calculate_true_ALT

logger = logging.getLogger()


##Import files containing annotation info and convert them to dictionary
#script_dir = os.path.dirname(os.path.realpath(__file__))


def extract_reference_vcf(input_vcf):
    """
    Read file until header ends and pick the first field corresponding to reference
    """
    with open(input_vcf, "r") as f:
        next_line = f.readline().strip()
        while next_line.startswith("#"):
            next_line = f.readline()
        
    reference = next_line.split()[0]
    return reference


def replace_reference(input_vcf, output, ref_old=False, ref_new="Chromosome" ):
    """
    This function replace all instances of a reference in a vcf file
    Depends on extract_reference_vcf
    190909 - Function now uses chromosome name in file and replaces it with term provided (default "Chromosome")
    """
    input_file = os.path.abspath(input_vcf)
    output_file = os.path.abspath(output)
    output_dir = os.path.dirname(output)

    check_create_dir(output_dir)

    if ref_old == False:
        ref_old = extract_reference_vcf(input_file)

    with open(input_file, 'r') as fi:
        with open(output_file, 'w') as fo:
            for line in fi:
                ref = ref_old + "\t"
                new = ref_new + "\t"
                line = line.replace(ref, new)
                fo.write(line)

def snpeff_annotation(args, vcf_file, database=False):
    #http://snpeff.sourceforge.net/SnpEff_manual.html
    # java -jar snpEff.jar -ud 0 -c path/to/snpEff/snpEff.config -stats SAMPLE.starts.annot.html
    #   Mycobacterium_tuberculosis_h37rv VCF_FILE > FILE.annot
   
    vcf_file = os.path.abspath(vcf_file)

    sample = os.path.basename(vcf_file).split(".")[0]
    file_name = (".").join(os.path.basename(vcf_file).split(".")[:-1])

    snpeff_config = get_snpeff_path()

    annotate_output_dir = obtain_output_dir(args, "Annotation")
    #output_dir = os.path.abspath(args.output)
    stat_name = sample + ".annot.html"
    annot_name = file_name + ".annot"

    stat_file = os.path.join(annotate_output_dir, stat_name)
    output_file = os.path.join(annotate_output_dir, annot_name)

    cmd = ["snpEff", "-ud", "0", "-c", snpeff_config, "-stats", stat_file, database, vcf_file]
    #execute_subprocess(cmd)
    with open(output_file, "w+") as outfile:
        #calculate coverage and save it in th eoutput file
        subprocess.run(cmd,
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)


def import_annot_to_pandas(vcf_file, sep='\t'):
    """
    Parse vcf outputted by snpEFF which adds the ANN field
    Dependences: calculate_ALT_AD
                calculate_true_ALT
    """
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
        
        ann_head = ["Allele","Annotation","Annotation_Impact","Gene_Name",
                    "Gene_ID","Feature_Type","Feature_ID","Transcript_BioType",
                    "Rank","HGVS.c","HGVS.p","cDNA.pos / cDNA.length",
                    "CDS.pos / CDS.length","AA.pos / AA.length","Distance",
                    "ERRORS / WARNINGS / INFO"]
        
        for index, data_row in dataframe.iterrows():
            info_fields = re.findall(r';*([a-zA-Z]{1,20})=', data_row.INFO)
            info_values = re.findall(r'-?\d+\.?\d*e?[+-]?\d{0,2}', data_row.INFO)
            ann_values_re = re.search(r'ANN=(.*)\|(.*)', data_row.INFO)
            all_ann_values = ann_values_re.group(1)
            ann_values = all_ann_values.split(",")[0].split("|")[:16]
            
            format_fields = data_row['FORMAT'].split(":")
            format_values = data_row['sample'].split(":")
                                    
            for ifield, ivalue in zip(info_fields,info_values):
                dataframe.loc[index,ifield] = ivalue
                
            for ffield, fvalue in zip(format_fields,format_values):
                dataframe.loc[index,ffield] = fvalue
            
            dataframe.loc[index,'ANN'] = all_ann_values
            
            for ann_field,ann_value in zip(ann_head, ann_values):
                dataframe.loc[index,ann_field] = ann_value
            
        dataframe.rename(columns={'AF':'af'}, inplace=True)
        
        dataframe['len_AD'] = dataframe['AD'].str.split(",").str.len()
        dataframe['REF_AD'] = dataframe['AD'].str.split(",").str[0]

        dataframe['ALT_AD'] = dataframe.apply(calculate_ALT_AD, axis=1)
        dataframe['ALT'] = dataframe.apply(calculate_true_ALT, axis=1)

        dataframe[['gt0','gt1']] = dataframe['GT'].str.split(r'[/|\|]', expand=True)

        dataframe['HGVS.c'] = dataframe['HGVS.c'].str.split(".").str[-1]
        dataframe['HGVS.p'] = dataframe['HGVS.p'].str.split(".").str[-1]
        dataframe['Gene length'] = dataframe['CDS.pos / CDS.length'].str.split("/").str[-1]
        dataframe['AA length'] = dataframe['AA.pos / AA.length'].str.split("/").str[-1]
        dataframe['AA length'] = dataframe['AA length'].replace('', 0)
        dataframe['HGVS.p'] = dataframe['HGVS.p'].replace('', 'None')
                
        to_float = ['QUAL', 'AC', 'af', 'AN', 'BaseQRankSum', 'DP', 'ExcessHet', 'FS',
       'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR','GQ','ALT_AD', 'REF_AD']
        
        to_int = ['POS', 'len_AD', 'gt0', 'gt1']
        
        to_str = ['#CHROM','REF','ALT', 'FILTER']
        
        for column in dataframe.columns:
            if column in to_float:
                dataframe[column] = dataframe[column].astype(float)
                
        for column in dataframe.columns:
            if column in to_int:
                dataframe[column] = dataframe[column].astype(int)
                
        for column in dataframe.columns:
            if column in to_str:
                dataframe[column] = dataframe[column].astype(str)
                
        dataframe['dp'] = (dataframe['REF_AD'] + dataframe['ALT_AD'])
        dataframe['aF'] = dataframe['REF_AD']/dataframe['dp']
        dataframe['AF'] = dataframe['ALT_AD']/dataframe['dp']
        

    else:
        print("This vcf file is not v4.2")
        sys.exit(1)
           
    return dataframe



def final_annotation(vcf_file_annot, *bed_files):
    """
    import annotated vcf with snpEff
    add Lineage info -> output final lineage to an external file
    add resistance info -> output final lineage to an external file
    """
    df_vcf = import_annot_to_pandas(vcf_file_annot, sep='\t')

    vcf_path = os.path.abspath(vcf_file_annot)
    output_dir = ("/").join(vcf_path.split("/")[:-1])
    vcf_name = vcf_path.split("/")[-1]

    tab_name = (".").join(vcf_name.split(".")[:-1])
    #extend_raw = ".raw.annot.tab"
    extend_final = ".annot.tsv"

    annotate_bed_s(df_vcf, *bed_files)

    #Add lineage info 
    #add_lineage_Coll(df_vcf)

    #Add resistance info

    #Retrieve only annotation fields
    #df_vcf_annot = df_vcf[['#CHROM', 'POS', 'ID', 'REF', 'ALT','Annotation',
    #   'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type',
    #   'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p',
    #   'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length', 'AA.pos / AA.length','Is_essential','Product', 'Lineage', 'Resistance']]
    
    #output all raw info into a file
    new_out_file = tab_name + extend_final
    output_raw_tab = os.path.join(output_dir, new_out_file)
    df_vcf.to_csv(output_raw_tab, sep='\t', index=False)
    

def get_reverse(nucleotyde):
    nucleotyde = str(nucleotyde)
    nucleotyde_rev = {'A' : 'T',
                     'T' : 'A',
                     'C' : 'G',
                     'G': 'C'}
    if len(nucleotyde) > 1:
        nucleotyde_str = nucleotyde[::-1] #Reverse nucleotide
        nucleotyde_str_fin = "".join([nucleotyde_rev[x] for x in nucleotyde_str]) #Complement nucleotide
        return nucleotyde_str_fin
    else:
        return nucleotyde_rev[nucleotyde]

css_report = """

    <style type="text/css">
    
    body {
        background-color: #cecccc;
        font: normal 20px Verdana, Arial, sans-serif;
        border: 1px solid black;
        border-radius: 5px;
        padding: auto;
        margin: auto;
        max-width: 1000px;
        min-width: 500px;
    }

    #center{
        background-color: white;
        border-radius: 5px;
        padding: 10px;
        display: block;
        text-align: center;
        height: 100%;
    }

    #info-div-center{
        display: block;
        text-align: center;
    }

    #info-text{
        display: inline-block;
        text-align: left;
    }

    p {
        font-size: 15px;
        font-weight: normal;
    }

    table {
        text-align: center;
        border-color: #000;
        border-spacing: 0px;
        border-style: solid;
        border-width: 1px;
        font-size: 0.75em;
        width: 100%;
    }

    th, td {
    border-bottom: 1px solid #ddd;
    padding:0 5px 0 5px;
    }

    th {
    background-color: rgb(76, 175, 170);
    }

    tr:nth-child(even) {background-color: #cecccc;}
    tr:hover {background-color:#7c7b7b;}

    footer{
        position: relative;
        width: 100%;
        display: flex;
        flex-flow: wrap;
    }
    footer p{
        padding-left: 10px;
        font-size: 0.75em;
        text-align: left;
    }
    @page {
        size: A4;
        margin: 0;
    }

    @media print {
    html, body {
        width: 210mm;
        height: 297mm;
        margin: 0px;
    }}

    </style>

    """

"""
#info-text::before{
        content: "";
        display: inline-block;
        height: 5px; 
        width: 20px;
        padding: 0.5px; 
        background-color: teal; 
        border-radius: 150px;  
        transform: rotate(200deg) skew(0deg); 
        -webkit-transform: rotate(200deg) skew(0deg); 
        border: 3px dashed #D8A3CA;
        margin-right: 0.5em;
        }
"""



def create_report(tab_annot, css=css_report, species="Mycobacterium tuberculosis", species_report="Main species: <i>Mycobacterium tuberculosis</i><br>"):
    #<div style="position: absolute; bottom: 5px; color: red; background-color: rgb(253, 253, 253)">
    #Text disclaimer 
    #</div>

    script_dir = os.path.dirname(os.path.realpath(__file__))
    annotation_dir = os.path.join(script_dir, "annotation/resistance")

    table_res = os.path.join(annotation_dir, "MTB_Resistance_Mediating.txt")
    df_res = pd.read_csv(table_res, sep="\t", header=0)
    df_res['High Confidence SNP'].fillna("no", inplace=True)


    output = os.path.dirname(tab_annot)
    sample = os.path.basename(tab_annot).split(".")[0]
    report_name = sample + ".annot.report.html"

    output_file = os.path.join(output, report_name)
    
    cummulative_report = ""

    with open(output_file, 'w+') as f:

        f.write(css)
        starter_div = """
        <div id = "center">
            <div id = "info-div-center">
                <div id = "info-text">
        """
        f.write(starter_div)
        line_sample = "Nombre de la muestra: " + sample + "<br> \
            <br>\n"
        f.write(line_sample)
        cummulative_report = cummulative_report + starter_div + line_sample

        line_species = "Especie: " + "<i>" + str(species) + "</i>" + "<br> \
            <br>\n"
        f.write(line_species)
        cummulative_report = cummulative_report + species_report + "<br>\n"

        df_annot = pd.read_csv(tab_annot, sep="\t", header=0)
        
        df_resistance = df_annot[df_annot.Resistance.notnull()]
        df_resistance_F = df_resistance[['POS', 'ALT', 'AF', 'Annotation', 'Gene_ID', 'Gene_Name', 'HGVS.c', 'HGVS.p', 'Resistance']]
        df_resistance_F.columns = ['Position', 'Alt. base', 'AF', 'Change', 'Gene ID', 'Gene Name', 'Codon change','AA change', 'Resistance']
        list_resistance = df_annot['Resistance'][df_annot.Resistance.notnull()].tolist()


        list_lineage = df_annot['Lineage'][df_annot.Lineage.notnull()].tolist()
        
        #Output Lineage info
        if len(list_lineage) > 0:
            list_lineage.sort(reverse=True)
            asterix = ""
            for sublineage_n in range(len(list_lineage)):
                if sublineage_n < (len(list_lineage) - 1):
                    if str(list_lineage[sublineage_n]).startswith(str(list_lineage[sublineage_n + 1])):
                        asterix = asterix + "*"
            final_lineage = str(list_lineage[0]) #+ " " + asterix
            line_lineage = "Esta cepa tiene los marcadores de Linaje: " + "<b>" + str(final_lineage) + "</b>" + "<br>\n \
                <br>\n"
            f.write(line_lineage)
            cummulative_report = cummulative_report + line_lineage
        else:
            line_lineage = "<br>\n"
                
            f.write(line_lineage)
            cummulative_report = cummulative_report + line_lineage
        
        #Output Resistance info
        if len(list_resistance) == 1:
            line_res_1 = "Esta cepa tiene " + str(len(list_resistance)) + " mutación asociada a resistencia:<br>\n \
                <br>\n \
                </div> \
                    </div>"
            f.write(line_res_1)
            cummulative_report = cummulative_report + line_res_1
            f.write(tabulate(df_resistance_F, headers='keys', tablefmt='html', showindex=False))
            table_res = tabulate(df_resistance_F, headers='keys', tablefmt='html', showindex=False)
            cummulative_report = cummulative_report + table_res + "\n"

            f.write("\n<br>\n \
            <footer> \
                <p>Los asteriscos (*) al final del campo 'Resistance' hacen referencia a posiciones de baja confianza. Se recomienda análisis fenotípico para los antibióticos señalados.</p>\n \
                <p>El campo AF indica la frecuencia de la mutación, siendo 1 el 100%</p>\n")

        elif len(list_resistance) > 1:
            line_res_1 = "Esta cepa tiene " + str(len(list_resistance)) + " mutaciones asociadas a resistencia:<br>\n \
                <br>\n \
                </div> \
                    </div>"
            f.write(line_res_1)
            cummulative_report = cummulative_report + line_res_1
            f.write(tabulate(df_resistance_F, headers='keys', tablefmt='html', showindex=False))
            table_res = tabulate(df_resistance_F, headers='keys', tablefmt='html', showindex=False)
            cummulative_report = cummulative_report + table_res + "\n"

            f.write("\n<br>\n \
            <footer> \
                <p>Los asteriscos (*) al final del campo 'Resistance' hacen referencia a posiciones de baja confianza. Se recomienda análisis fenotípico para los antibióticos señalados.</p>\n \
                <p>El campo AF indica la frecuencia de la mutación, siendo 1 el 100%</p>\n")
                
        else:
            f.write("No se han encontrado mutaciones asociadas a resistencia<br> \
                </div> \
                    </div> \
                    \n<br>\n \
                <footer>")
            cummulative_report = cummulative_report + "No Resistance positions were found<br>\n \
                </div> \
            </div> \
            </div>"
            

        f.write("\n<br>\n \
                <p>Este informe debe ser utilizado exclusivamente con fines de investigación. No utilizar ningún dato con fines asistenciales.</p>\n \
                <br>\n \
            </footer>\n \
            </div>")

    return cummulative_report

if __name__ == '__main__':
    print("#################### ANNOTATION #########################")