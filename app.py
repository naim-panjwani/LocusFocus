import requests
import json
import pandas as pd
import numpy as np
import tokens

def parseR2(urltext):
    temp = urltext.strip().split('\n')
    for t in temp:
        if t.strip().find('R2') != -1:
            return(float(t.strip().split(':')[1].strip()))
    return(np.nan)

lead_snp = 'rs7512462'
snp_list = ['rs61814953', 'rs1342063']
european_populations = ['CEU', 'TSI', 'FIN', 'GBR', 'IBS']
ld_type = 'r2'
base_url = 'https://ldlink.nci.nih.gov/LDlinkRest/ldpair?'

r2_pairs = []
for snp in snp_list:
    snp_query = 'var1=' + lead_snp + '&var2=' + snp
    population_query = '&pop=' + "%2B".join(european_populations) 
    #ld_query = '&r2_d=' + ld_type 
    token_query = '&token=' + tokens.token
    url = base_url + snp_query + population_query + token_query
    response = requests.get(url).text
    r2 = parseR2(response)
    #print(r2)
    r2_pairs.append(r2)
    
