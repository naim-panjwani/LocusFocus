#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 21:20:01 2021

@author: naim
Adapted from: https://maxhalford.github.io/blog/pandas-streaming-groupby/
Works well in Ubuntu
"""

import itertools
# import multiprocessing as mp # Mongo does not support several write requests
import pandas as pd
import numpy as np
import os
import pymongo
from pymongo import MongoClient
import subprocess
from datetime import datetime


def stream_groupby_csv(path, key, agg, chunk_size=1e6, pool=None, **kwargs):

    # Make sure path is a list
    if not isinstance(path, list):
        path = [path]

    # Chain the chunks
    kwargs['chunksize'] = chunk_size
    chunks = itertools.chain(*[
        pd.read_csv(p, **kwargs)
        for p in path
    ])

    results = []
    orphans = pd.DataFrame()
        
    for chunk in chunks:

        # Add the previous orphans to the chunk
        chunk = pd.concat((orphans, chunk))

        # Determine which rows are orphans
        last_val = chunk[key].iloc[-1]
        is_orphan = chunk[key] == last_val

        # Put the new orphans aside
        chunk, orphans = chunk[~is_orphan], chunk[is_orphan]

        # If a pool is provided then we use apply_async
        if pool:
            results.append(pool.apply_async(agg, args=(chunk,)))
        else:
            results.append(agg(chunk))

    # If a pool is used then we have to wait for the results
    if pool:
        results = [r.get() for r in results]
    
    results.append(agg(orphans)) # ensure last chunk (gene) is pushed as well!

    return pd.concat(results)


def agg(chunk):
    """lambdas can't be serialized so we need to use a function"""
    chunk.set_index('gene_id', inplace=True)
    return chunk.groupby('gene_id').apply(push_variant_dict)


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


##########################################
# MAIN
##########################################

conn = "mongodb://localhost:27017"
client = MongoClient(conn)
db = client.GTEx_V7

tissues = pd.read_csv(os.path.join('data', 'GTEx_v7_eQTL','tissues4.txt'), header=None)
tissues = list(tissues.iloc[:,0])
files_list = [ tissue.replace(' ','_') + '.allpairs.txt.gz' for tissue in tissues ]
# files_list = [ 'Pancreas.allpairs_fixed.txt.gz', 'Lung.allpairs_fixed.txt.gz' ]


for file in files_list:
    #file = 'brain_test.allpairs.txt.gz'
    tissue_name = file.split('.')[0].replace(' ','_')
    file = os.path.join('data','GTEx_v7_eQTL', file)
    if tissue_name not in db.list_collection_names():
        collection = db[tissue_name]
        if file.endswith('gz') and os.path.isfile(file):
            print('Decompressing ' + file)
            subprocess.run(args=['gunzip','-f',file])
            #decompress(file)
        print('Parsing file ' + file + ' and creating tissue collection')
        results = stream_groupby_csv(
            path=[ file.replace('.gz','') ],
            key='gene_id',
            agg=agg,
            chunk_size=1e6,
            #pool=mp.Pool(processes=4), # problem with opening too many write requests to mongo
            sep="\t", 
            usecols=['gene_id', 'variant_id', 'pval_nominal',
                     'slope', 'slope_se', 'ma_samples', 'ma_count', 'maf']
        )
        print(tissue_name + ' collection created')
        print(datetime.now().strftime('%c'))
        print('Now indexing by gene_id')
        print(datetime.now().strftime('%c'))
        collection.create_index('gene_id')
        print('Indexing done')
        print(datetime.now().strftime('%c'))
        print('Recompressing ' + file.replace('.gz',''))
        subprocess.run(args=['gzip', file.replace('.gz','')])
        #print('Deleting ' + file.replace('.gz',''))
        #subprocess.run(args=['rm', '-f', file.replace('.gz','')])
        print('Done with tissue ' + tissue_name)
        print(datetime.now().strftime('%c'))
    #break

# Next, create the variant lookup table
print('Reading variant lookup file GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz')
collection = db['variant_table']
tbl_chunk = pd.read_csv(os.path.join('data','GTEx_v7_eQTL','GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz'), 
                        sep="\t", chunksize=1e5, encoding='utf-8')
print('Pushing variant information into variant_table collection by chunks')
for tbl in tbl_chunk:
    tbl['chr'] = [int(str(x).replace('chr','').replace('X','23')) for x in list(tbl['chr'])]
    tbl['variant_id'] = [x.replace('chr','') for x in list(tbl['variant_id'])]
    tbl_dict = tbl.to_dict(orient='records')
    collection.insert_many(tbl_dict)
print('Variant collection created')
print(datetime.now().strftime('%c'))

print('Now indexing by variant_id')
collection.create_index('variant_id')
print('Indexing by variant_id done')
print(datetime.now().strftime('%c'))
print('Now indexing by chr and pos (ascending) order')
collection.create_index([("chr", pymongo.ASCENDING),
                         ("variant_pos", pymongo.ASCENDING)])
print('Indexing by chr and pos done')
print(datetime.now().strftime('%c'))



# Test database results ok:
# tissue_name = 'test'
# collection = db[tissue_name]
# collection.estimated_document_count()
# gene = 'ENSG00000241860.2'
# results = collection.find({'gene_id': gene})
# temp = list(results)
# len(temp)
# temp[0]['gene_id'] == gene
# temp[0]['eqtl_variants']
# len(temp[0]['eqtl_variants'])
