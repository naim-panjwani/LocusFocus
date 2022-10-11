# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 19:59:32 2020

@author: Naim
"""

import pandas as pd

htmlpage = pd.read_html('https://www.ncbi.nlm.nih.gov/grc/human/data')
chromlengths = htmlpage[0].iloc[:,0:2]
chromlengths.columns = ['sequence', 'length']
chromlengths['sequence'] = [ 'chr' + str(x) for x in list(chromlengths['sequence'])]
chromlengths.to_csv('hg38_chrom_lengths.txt', index=False, encoding='utf-8', sep="\t")
