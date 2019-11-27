#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 17:43:52 2019

@author: naim
"""

import argparse
import numpy as np
import pandas as pd
import gzip
from datetime import datetime
from tqdm import tqdm
import sys


##############################
## Helper functions
##############################



##############################
## MAIN
##############################
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Merge several datasets together into HTML tables separated by <h3> title tags')
    parser.add_argument('filelist_filename', 
                        help='Filename containing the list of files to be merged together.\n \
                        The second column (tab-delimited) may contain descriptions of the datasets.\n \
                        The remaining columns specify the column names for chromosome, basepair position, SNP name, P-value (in that order)')
    parser.add_argument('outfilename', help="Desired output filename for the merged file")
    args = parser.parse_args()
    
    filelist_filename = str(args.file_list)
    outfilename = str(args.outfilename.replace('.html','') + '.html')
    logfilename = str(outfilename.replace('.html','') + '.log')

    old_stdout = sys.stdout
    log_file = open(logfilename, "w")
    sys.stdout = log_file
    
    print(datetime.now().strftime('%c'))
    print('filelist_filename: ' + filelist_filename)
    print('outfilename: ' + outfilename)
    print('Log file: ' + logfilename)

    print('Reading list of files file')
    filelist = pd.read_csv(filelist_filename, sep='\t', header=0)
    files = list(filelist.iloc[:,0])
    file_descriptions = list(filelist.iloc[:,1])    
    
    

    print(datetime.now().strftime('%c'))
    sys.stdout = old_stdout
    log_file.close()
    