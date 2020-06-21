#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 18:28:02 2019

@author: naim

This script simply takes the gtf format of the collapsed GENCODE genes 
downloaded from GTEx and reformats it as a custom dataframe
"""

import pandas as pd
import numpy as np

def rows_to_skip(afile, num_lines_to_check = 1000):
    num_header_lines = 0
    num_lines_checked = 0
    with open(afile, 'r') as f:
        nextline = f.readline()
        num_lines_checked += 1
        while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
            num_header_lines += 1
            nextline = f.readline()
            num_lines_checked += 1
    return num_header_lines


def convert_gtf(gencode):
    col8 = gencode.iloc[:,8]
    ENSG_name = [i.split(';')[0].strip().split(" ")[1].replace('\"','') for i in col8]
    name = [i.split(';')[3].strip().split(" ")[1].replace('\"','') for i in col8]
    chrom = gencode.iloc[:,0]
    if str(chrom[0])[0:3] != 'chr':
        chrom = ['chr' + str(i) for i in chrom]
    
    gencode = pd.concat([gencode, pd.Series(ENSG_name), pd.Series(name), pd.Series(chrom)], axis=1)
    gencode.columns = np.arange(12)    
    
    gene_grp = gencode.groupby(9)
    
    ENSG_name_col = []
    name_col = []
    chrom_col = []
    strand_col = []
    txStart_col = []
    txEnd_col = []
    exonStarts_col = []
    exonEnds_col = []
    
    for gene in gene_grp:
        gene_df = gene[1]
        txStart = int(gene_df.loc[gene_df[2]=="gene"][3])
        txEnd = int(gene_df.loc[gene_df[2]=="gene"][4])
        exonStarts = []
        exonEnds = []
        exons_df = gene_df.loc[gene_df[2]=="exon"]
        exonStarts = list(exons_df[3])
        exonEnds = list(exons_df[4])
        exonStarts_str = ",".join([str(x) for x in exonStarts])
        exonEnds_str = ",".join([str(x) for x in exonEnds])
        ENSG_name_col.append(list(gene_df[9])[0])
        name_col.append(list(gene_df[10])[0])
        chrom_col.append(list(gene_df[11])[0])
        strand_col.append(list(gene_df[6])[0])
        txStart_col.append(txStart)
        txEnd_col.append(txEnd)
        exonStarts_col.append(exonStarts_str)
        exonEnds_col.append(exonEnds_str)
    
    new_gencode_df = pd.DataFrame({
            'ENSG_name': ENSG_name_col,
            'name': name_col,
            'chrom': chrom_col,
            'strand': strand_col,
            'txStart': txStart_col,
            'txEnd': txEnd_col,
            'exonStarts': exonStarts_col,
            'exonEnds': exonEnds_col
            })
    
    return new_gencode_df


################################
# MAIN
################################
gencode19 = pd.read_csv("gencode.v19.genes.v7.patched_contigs.gtf", 
                      skiprows = rows_to_skip('gencode.v19.genes.v7.patched_contigs.gtf'), 
                      header=None, sep="\t", encoding="utf-8")
new_gencode_df19 = convert_gtf(gencode19)
new_gencode_df19.to_csv('collapsed_gencode_v19_hg19.gz', encoding='utf-8', index=False, sep="\t", compression='gzip')

gencode26 = pd.read_csv("gencode.v26.GRCh38.genes.gtf",
                        skiprows = rows_to_skip('gencode.v26.GRCh38.genes.gtf'),
                        header=None, sep="\t", encoding="utf-8")
new_gencode_df26 = convert_gtf(gencode26)
new_gencode_df26.to_csv('collapsed_gencode_v26_hg38.gz', encoding='utf-8', index=False, sep="\t", compression='gzip')
