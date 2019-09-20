#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 18:28:02 2019

@author: naim
"""

import pandas as pd
import numpy as np

gencode = pd.read_csv("gencode.v19.genes.v7.patched_contigs_mod.gtf", header=None, sep="\t", encoding="utf-8")
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

new_gencode_df.to_csv('collapsed_gencode_v19_hg19', encoding='utf-8', index=False, sep="\t")
