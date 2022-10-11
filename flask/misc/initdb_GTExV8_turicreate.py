#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 18:12:05 2020

@author: naim
Does not work in Ubuntu - groupby requires large amount of memory; turicreate does not handle properly
"""

import turicreate as tc
import pandas as pd
import numpy as np
# import pymongo
from pymongo import MongoClient
# from datetime import datetime
import subprocess
import os
import gzip
import shutil

#passwd = pd.read_csv('.passwd', header=None, encoding='utf-8').iloc[0,0]

tissues = pd.read_csv(os.path.join('data', 'GTEx_v8_eQTL','tissues.txt'), header=None)
tissues = list(tissues.iloc[:,0])
#files_list = [ 'GTEx_Analysis_v8_eQTL_all_associations_' + tissue.replace(' ','_') + '.allpairs_fixed.txt.gz' for tissue in tissues ]
files_list = [ 'Pancreas.allpairs_fixed.txt.gz', 'Lung.allpairs_fixed.txt.gz' ]

# Initialize MongoDB database
conn = "mongodb://localhost:27017"
client = MongoClient(conn)
db = client.GTEx_V8

def push_variant_dict(gene_df):
    variant_id = [str(x).replace('chr','').encode('utf-8') for x in list(gene_df['variant_id'])]
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


def decompress(f):
    with gzip.open(f, 'rb') as f_in:
        with open(f.replace('.gz',''), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)




# for file in files_list:
file = files_list[1]
tissue_name = file.split('.')[0].replace(' ','_').replace('GTEx_Analysis_v8_eQTL_all_associations_','')
file = os.path.join('data','GTEx_v8_eQTL', file)
if tissue_name not in db.list_collection_names():
    collection = db[tissue_name]
    if file.endswith('gz') and os.path.isfile(file):
        print('Decompressing ' + file)
        subprocess.run(args=['gunzip','-f',file])
    print('Reading file ' + file)
    # tissue_eqtls = tc.SFrame.read_csv(file.replace('.gz',''), delimiter = '\t', header=True, 
    #                                   usecols=['gene_id', 'variant_id', 'pval_nominal',
    #                                             'slope', 'slope_se', 'ma_samples', 'ma_count', 'maf'])
    # tissue_eqtls.save(tissue_name + '.sframe')
    tissue_eqtls = tc.load_sframe(tissue_name + '.sframe')
    print('Parsing file ' + file + ' and creating tissue collection')       
    other_cols = ['variant_id', 'pval_nominal','slope', 'slope_se', 'ma_samples', 'ma_count', 'maf']
    agg_list = [tc.aggregate.CONCAT(i) for i in other_cols]
    gene_groups = tissue_eqtls.groupby('gene_id', dict(zip(other_cols, agg_list)))
    print(gene_groups)



