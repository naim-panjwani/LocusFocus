#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 09:59:02 2019

@author: naim
"""

import pandas as pd
import numpy as np
from pymongo import MongoClient
import pymongo
from datetime import datetime
import dask.dataframe as dd
import subprocess
import os

#passwd = pd.read_csv('.passwd', header=None, encoding='utf-8').iloc[0,0]

files_list = pd.read_csv('tissues.txt', header=None)
files_list = list(files_list.iloc[:,0])
files_list = [file.replace('.txt.gz','_fixed.txt.gz') for file in files_list]

# Initialize MongoDB database
conn = "mongodb://localhost:27017"
client = MongoClient(conn)
db = client.GTEx_V7

def push_variant_dict(gene_df):
    variant_id = [str(x).encode('utf-8') for x in list(gene_df['variant_id'])]
    pval = list(gene_df['pval_nominal'])
    beta = list(gene_df['slope'])
    se = list(gene_df['slope_se'])
    ma_samples = list(gene_df['ma_samples'])
    ma_count = list(gene_df['ma_count'])
    sample_maf = list(gene_df['maf'])
    geneid = gene_df.reset_index()['gene_id'][0]
    variants_list = []
    for row in np.arange(len(variant_id)):
        variants_list.append({
            'variant_id': variant_id[row].decode('utf-8')
            , 'pval': float(pval[row])
            , 'beta': float(beta[row])
            , 'se': float(se[row])
            , 'ma_samples': float(ma_samples[row])
            , 'ma_count': float(ma_count[row])
            , 'sample_maf': float(sample_maf[row])                     
            })
    gene_dict = {'gene_id': geneid, 'eqtl_variants': variants_list }
    collection.insert_one(gene_dict)

for file in files_list:
    tissue_name = file.split('.')[0].replace(' ','_')
    if tissue_name not in db.list_collection_names():
        collection = db[tissue_name]
        if file.endswith('gz') and os.path.isfile(file):
            print('Decompressing ' + file)
            subprocess.run(args=['gunzip','-f',file])
        print('Reading file ' + file)
        tissue_eqtls = dd.read_csv(file.replace('.gz',''), sep="\t", 
                                   usecols=['gene_id', 'variant_id', 'pval_nominal',
                                            'slope', 'slope_se', 'ma_samples', 'ma_count', 'maf'])
        tissue_eqtls = tissue_eqtls.set_index('gene_id')
        genes_dict_list = []
        #tissue_eqtls = tissue_eqtls.loc[:, ~tissue_eqtls.columns.isin(['chr','pos'])]
        print('Parsing file ' + file + ' and creating tissue collection')
        print(datetime.now().strftime('%c'))
        tissue_eqtls.groupby('gene_id').apply(push_variant_dict).compute()
#        for i in np.arange(len(result)):
#            genes_dict_list.append({
#                    'gene_id': result.keys()[i]
#                    ,'eqtl_variants': result[i]
#                    })
        print(tissue_name + ' collection created')
        print(datetime.now().strftime('%c'))
        print('Now indexing by gene_id')
        print(datetime.now().strftime('%c'))
        collection.create_index('gene_id')
        print('Indexing done')
        print(datetime.now().strftime('%c'))
        #print('Recompressing ' + file.replace('.gz',''))
        #subprocess.run(args=['gzip', file.replace('.gz','')])
        print('Deleting ' + file.replace('.gz',''))
        subprocess.run(args=['rm', '-f', file.replace('.gz','')])
        print('Done with tissue ' + tissue_name)
        print(datetime.now().strftime('%c'))

# Next, create the variant lookup table
collection = db['variant_table']
tbl_chunk = pd.read_csv('GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz',
                  sep="\t", chunksize=100, encoding='utf-8')
for tbl in tbl_chunk:
    tbl['chr'] = [int(str(x).replace('X','23')) for x in list(tbl['chr'])]
    tbl_dict = tbl.to_dict(orient='records')
    collection.insert_many(tbl_dict)
    
collection.create_index('variant_id')
collection.create_index([("chr", pymongo.ASCENDING),
                         ("variant_pos", pymongo.ASCENDING)])


#results = collection.find({'gene_id': gene})
#temp = list(results)
#len(temp)
#temp[0]['gene_id'] == gene
#temp[0]['eqtl_variants']

