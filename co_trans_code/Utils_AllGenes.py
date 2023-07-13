#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import numpy as np
import pandas as pd
import os
import io


''' These Utils are related to the TCGA data of all protein-coding genes (data is on the power cluster at cancer_proj/gdc_db/data/filtered_feb_2021/AllGenes). 
The variant data is stored in maf files and an example for the data structure is  cancer_proj/gdc_db/data/filtered_feb_2021/AllGenes/ABCB1/BRCA.maf, meaning that the ABCB1 directory holds all the TCGA mutations of the ABCB1 gene and in the BRCA.maf file we can find the ABCB1 mutations of the BRCA patients. '''



''' The maf files have varying number of header rows. This function gets the number of header lines of a specific file, 
so we can input it to pd.read_csv and get the maf file as a table'''
def get_rows_to_skip(path):
    header_counter = 0
    with open(path, "r") as a_file:
        for line in a_file:
            if line[0]!= '#': #header rows start with the "#" sign
                break
            header_counter += 1
    return(header_counter)


''' Get all TCGA variant data for a specific gene'''
def get_cancer_muts_cur_gene(gene_name):
    cancerous_muts_path = f"../../../../../tamir1/cancer_proj/gdc_db/data/filtered_feb_2021/AllGenes/{gene_name}/"
    #iterate over all the MAF files in the folder of this gene and get all synonymous mutations
    cancer_muts_cur_gene = pd.DataFrame()
    for filename in os.listdir(cancerous_muts_path):
        if filename.endswith(".maf"):
            path_cur_cancer_type = f"{cancerous_muts_path}{filename}"
            num_rows_to_skip = get_rows_to_skip(path_cur_cancer_type) #the number of header notes is different for each file
            temp_df = pd.read_csv(path_cur_cancer_type , sep = "\t", skiprows = num_rows_to_skip)
            temp_df = temp_df[temp_df['Variant_Classification'] == 'Silent'] #keep only synonymous variants
            temp_df['Cancer_Type'] = filename.split(".")[0] #keep the data of the cancer type
            cancer_muts_cur_gene = pd.concat([cancer_muts_cur_gene, temp_df], sort=False)
    return(cancer_muts_cur_gene)

'''The list of TCGA mutations in a gene also contains mutations in the stop codon. However, our CDSs in this project do not 
contain the stop codon. So, we need to remove mutations in the stop codon from the list of mutations obtained from TCGA.'''

def remove_stop_codon_muts(cancer_muts_cur_gene:pd.DataFrame) -> pd.DataFrame:
    cancer_muts_cur_gene['CDS'] = cancer_muts_cur_gene['CDS_position'].apply(lambda x: int(x.split("/")[0]))
    len_CDS = int(cancer_muts_cur_gene.iloc[0]['CDS_position'].split("/")[1]) 
    cancer_muts_cur_gene['not_in_stop_codon'] = len_CDS - cancer_muts_cur_gene['CDS'] > 2 #"len CDS" is the cds length including the
    #stop codon. "cur_gene_syn_muts['CDS']" is the position of the mutatation. if len_CDS - cur_gene_syn_muts['CDS'] <= 2 the mutation
    #is in the stop codon. We remove these. 
    cancer_muts_cur_gene = cancer_muts_cur_gene[cancer_muts_cur_gene['not_in_stop_codon']]
    return(cancer_muts_cur_gene)


def read_vcf(path): #from https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})