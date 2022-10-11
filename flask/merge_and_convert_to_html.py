#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 17:43:52 2019

@author: naim
For usage instructions, simply run (python version 3):
python merge_and_convert_to_hmtl.py -h
"""

import argparse
import numpy as np
import pandas as pd
import gzip
from datetime import datetime
from tqdm import tqdm
import sys
#from bs4 import BeautifulSoup as bs

# Default column names:
CHROM = 'CHROM'
BP = 'BP'
SNP = 'SNP'
P = 'P'

genomicWindowLimit = 2e6

##############################
## Helper functions
##############################
def parseRegionText(regiontext):
    regiontext = regiontext.strip().replace(' ','')
    chrom = regiontext.split(':')[0].replace('chr','').replace('Chr','')
    pos = regiontext.split(':')[1]
    startbp = pos.split('-')[0].replace(',','')
    endbp = pos.split('-')[1].replace(',','')
    #chromLengths = pd.read_csv(os.path.join(MYDIR, 'data/hg19_chrom_lengths.txt'), sep="\t", encoding='utf-8')
    #chromLengths.set_index('sequence',inplace=True)
    if chrom == 'X':
        chrom = 23
        #maxChromLength = chromLengths.loc['chrX', 'length']
        try:
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise Exception("Invalid coordinates input")
    else:
        try:
            chrom = int(chrom)
            #maxChromLength = chromLengths.loc['chr'+str(chrom), 'length']
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise Exception("Invalid coordinates input")
    if chrom < 1 or chrom > 23:
        raise Exception('Chromosome input must be between 1 and 23')
    elif startbp > endbp:
        raise Exception('Starting chromosome basepair position is greater than ending basepair position')
#    elif startbp > maxChromLength or endbp > maxChromLength:
#        raise Exception('Start or end coordinates are out of range')
    else:
        return chrom, startbp, endbp


def checkColnames(filename, cols):
    '''
    Checks if the 'filename' has the column names in cols (a list of column names); 
    returns False if any column is absent
    '''
    num_lines_to_check = 1000
    num_lines_checked = 0
    if filename.endswith('gz'):
        with gzip.open(filename,'rb') as f:
            nextline = f.readline().decode('utf-8').strip()
            num_lines_checked += 1
            while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
                nextline = f.readline().decode('utf-8').strip()
                num_lines_checked += 1
            headerline = nextline.split('\t')
            if all(col in headerline for col in cols):
                return True
            else:
                return False
    else:
        with open(filename,'r') as f:
            nextline = f.readline().strip()
            num_lines_checked += 1
            while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
                nextline = f.readline().strip()
                num_lines_checked += 1
            headerline = nextline.split('\t')
            if all(col in headerline for col in cols):
                return True
            else:
                return False
        
def getNumHeaderLines(file, num_lines_to_check = 1000):
    num_header_lines = 0
    num_lines_checked = 0
    if file.endswith('gz'):
        with gzip.open(file, 'rb') as f:
            nextline = f.readline().decode('utf-8')
            num_lines_checked += 1
            while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
                num_header_lines += 1
                nextline = f.readline().decode('utf-8')
                num_lines_checked += 1
    else:
        with open(file, 'r') as f:
            nextline = f.readline()
            num_lines_checked += 1
            while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
                num_header_lines += 1
                nextline = f.readline()
                num_lines_checked += 1
    return num_header_lines

class Logger(object):
    def __init__(self, outfilename):
        self.terminal = sys.stdout
        self.log = open(outfilename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        pass    


##############################
## MAIN
##############################
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Merge several datasets together into HTML tables separated by <h3> title tags')
    parser.add_argument('filelist_filename', 
                        help='Filename containing the list of files to be merged together.\n \
                        The second column (tab-delimited) may contain descriptions of the datasets.\n \
                        The remaining columns specify the column names for chromosome, basepair position, SNP name, P-value (in that order).')
    parser.add_argument('coordinates', help="The region coordinates to subset from each file (e.g. 1:500,000-600,000")
    parser.add_argument('outfilename', help="Desired output filename for the merged file")
    args = parser.parse_args()  

    filelist_filename = str(args.filelist_filename)
    chrom, startbp, endbp = parseRegionText(args.coordinates)
    outfilename = str(args.outfilename.replace('.html','') + '.html')
    logfilename = str(outfilename.replace('.html','') + '.log')

#    filelist_filename = "files_test.txt"
#    outfilename = "slc9a3_gwas_mqtl.html"
#    coordinates = 'chr5:384,664-612803'
#    chrom, startbp, endbp = parseRegionText(coordinates)

    old_stdout = sys.stdout
    sys.stdout = Logger(logfilename)
    
    print('Start: ' + datetime.now().strftime('%c'))
    print('filelist_filename: ' + filelist_filename)
    print('Region: ' + str(chrom)+':'+str(startbp)+'-'+str(endbp))
    print('outfilename: ' + outfilename)
    print('Log file: ' + logfilename)

    print('Reading list of files file')
    filelist = pd.read_csv(filelist_filename, sep='\t', header=None)
    files = list(filelist.iloc[:,0])
    file_descriptions = list(filelist.iloc[:,1])
    chrom_colnames = list(filelist.iloc[:,2])
    bp_colnames = list(filelist.iloc[:,3])
    snp_colnames = list(filelist.iloc[:,4])
    pval_colnames = list(filelist.iloc[:,5])
    
    print('Verifying settings in filelist_filename')
    for i in np.arange(len(files)):
        if not checkColnames(files[i], list(filelist.iloc[i,2:6])):
            raise Exception('Column names given: ' + str(list(filelist.iloc[i,2:6])) + '\n Not all match for file ' + files[i])

    print('Merging files')
    final_merge = '<!DOCTYPE html>\n<html>'
    for i in tqdm(np.arange(len(files))):
        df = pd.read_csv(files[i], sep='\t', skiprows=getNumHeaderLines(files[i]))
        desired_cols = list(filelist.iloc[i,2:6])
        df = df[desired_cols].rename(columns={
                desired_cols[0]: CHROM
                ,desired_cols[1]: BP
                ,desired_cols[2]: SNP
                ,desired_cols[3]: P
                })
        if chrom == 23:
            df[CHROM] = np.repeat(chrom, df.shape[0])
        df[CHROM] = [ int(str(x).lower().replace('chr','').replace('chrom','')) for x in list(df[CHROM]) ]
        df = df.loc[ (df[CHROM]==chrom) & (df[BP]>=startbp) & (df[BP]<=endbp) ]
        h3tag = '<h3>'+file_descriptions[i]+'</h3>'
        df_html_str = df.to_html(index=False)
        final_merge += h3tag + df_html_str
    final_merge += '</html>'
    
    print('Writing merged output')
    with open(outfilename, 'w') as f_out:
        f_out.write(final_merge)
    
    # Quick check:
#    with open(outfilename, encoding='utf-8', errors='replace') as f:
#        html = f.read()
#    soup = bs(html, 'lxml')
#    table_titles = soup.find_all('h3')
#    table_titles = [x.text for x in table_titles]
#    tables = pd.read_html(outfilename)

    print('Done')
    print('End: ' + datetime.now().strftime('%c'))
    sys.stdout = old_stdout

    
