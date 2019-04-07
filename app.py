import json
import requests
import pandas as pd
import numpy as np
import os
import math
import string
from tqdm import tqdm
import uuid
import tokens

import sqlalchemy as sa
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine, inspect, String, Integer

from flask import Flask, request, redirect, url_for, jsonify, render_template
from flask_sqlalchemy import SQLAlchemy
import pymysql
pymysql.install_as_MySQLdb()
#thepwd = open('pwd.txt').readline().replace('\n', '')

genomicWindowLimit = 2000000
fileSizeLimit = 200e6

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'upload'

####################################
# LD Querying from ldlink.nci.nih.gov
####################################

def parseR2(responseText):
    temp = responseText.strip().split('\n')
    for t in temp:
        if t.strip().find('R2') != -1:
            catchR2 = t.strip().split(':')[1].strip()
            try:
                result = float(catchR2)
            except:
                result = np.nan
            return result
    return np.nan

def parsePosition(responseText, snpname):
    temp = responseText.strip().split('\n')
    for t in temp:
        if t.strip().find(snpname) != -1:
            catchPOS = t.split(':')[1].replace(')','').strip()
            try:
                result = int(catchPOS)
            except:
                result = np.nan
            return result
    return np.nan

def queryLD(lead_snp, snp_list, populations=['CEU', 'TSI', 'FIN', 'GBR', 'IBS'], ld_type='r2'):
    base_url = 'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?'
    result_df = pd.DataFrame({'snp':[],'ld':[]})
    allowed_requests = 300 # GET request limit (POST not working; tried postman to verify as well)
    num_requests = math.ceil(len(snp_list) / (allowed_requests-1))
    for req in tqdm(np.arange(num_requests)):
        ld_values = []
        snps = []
        curr_snp_list = []
        start = req * allowed_requests
        end = min(start+allowed_requests-1, len(snp_list))
        curr_snp_list = snp_list[start:end]
        if lead_snp not in curr_snp_list:
            curr_snp_list.append(lead_snp)
        
        # Remove or solve problematic SNP names:
        allowed_chars = 'rs' + string.digits
        for i in np.arange(len(curr_snp_list)):
            snp = curr_snp_list[i]
            querysnp = snp.split(';')[0]
            curr_snp_list[i] = querysnp
            if any([char not in allowed_chars for char in querysnp]):
                curr_snp_list.pop(i)
        
        # Check if any SNPs are still left:
        if len(curr_snp_list) == 0:
            raise InvalidUsage('The provided file does not have any valid SNP names', status_code=410)
        
        ## POST requests (can do up to 1000 SNP batch requests) currently not working:
        # headers = {
        # 'Content-Type': 'application/json',
        # }
        # params = (
        #     ('token', tokens.token),
        # )
        # snpstring = "\\n".join(curr_snp_list)
        # data = f'{{"snps": {snpstring}, "pop": "CEU","r2_d": "d"}}'
        # response = requests.post('https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix', headers=headers, params=params, data=data, verify=False)

        # Do GET requests instead (max of 300 queries)
        snp_query = 'snps=' + '%0A'.join(curr_snp_list)
        population_query = '&pop=' + "%2B".join(populations)
        ld_query = '&r2_d=' + ld_type
        token_query = '&token=' + tokens.token
        url = base_url + snp_query + population_query + ld_query + token_query
        response = requests.get(url)
        if(response):
            data = response.text.strip().split('\n')[1:]
            tempSNPs = response.text.strip().split('\n')[0].split('\t')[1:]
            if lead_snp in tempSNPs:
                lead_snp_col = tempSNPs.index(lead_snp)
                data.pop(lead_snp_col)
            snps = [datum.strip().split('\t')[0] for datum in data if datum != '']
            ld_values = [datum.split('\t')[lead_snp_col+1] for datum in data if datum != '']
            for i in np.arange(len(ld_values)):
                try:
                    ld_values[i] = float(ld_values[i])
                except:
                    ld_values[i] = np.nan
            result_df = result_df.append(pd.DataFrame({'snp': snps, 'ld': ld_values}), ignore_index=True)
        snp_df = pd.DataFrame({'snp':snp_list})
        merged_df = snp_df.reset_index().merge(result_df, how='left', on='snp', sort=False).sort_values('index').drop(columns=['index']) # to preserve order of input snps
    return merged_df

# Very slow queryLD:
# def queryLD(lead_snp, snp_list, populations=['CEU', 'TSI', 'FIN', 'GBR', 'IBS'], ld_type='r2'):
#     base_url = 'https://ldlink.nci.nih.gov/LDlinkRest/ldpair?'
#     ld_values = []
#     positions = []
#     for snp in tqdm(snp_list):
#         snp_query = 'var1=' + lead_snp + '&var2=' + snp
#         population_query = '&pop=' + "%2B".join(populations) 
#         ld_query = '&r2_d=' + ld_type
#         token_query = '&token=' + tokens.token
#         url = base_url + snp_query + population_query + ld_query + token_query
#         response = requests.get(url).text
#         r2 = parseR2(response)
#         pos = parsePosition(response, snp)
#         ld_values.append(r2)
#         positions.append(pos)
#     return ld_values, positions


# from timeit import default_timer as timer
# start = timer()
# ld = queryLD(lead_snp, snp_list)
# end = timer()
# print(end - start)


####################################
# HG19 positions' querying to UCSC's snp151 table
####################################
def getSNPPositions(snp_list):
    pos = []
    allowed_chars = 'rs' + string.digits
    engine = sa.create_engine('mysql://genome@genome-mysql.cse.ucsc.edu:3306/hg19')
    for snp in tqdm(snp_list):
        querysnp = snp.split(';')[0]
        if any([char not in allowed_chars for char in querysnp]):
            pos.append(np.nan)
        else:
            pos.append(engine.execute(f"select chromEnd from snp151 where name='{querysnp}'").fetchall()[0][0])
    return pos


####################################
# Helper functions
####################################
def parseRegionText(regiontext):
    chr = regiontext.split(':')[0].replace('chr','')
    pos = regiontext.split(':')[1]
    startbp = pos.split('-')[0]
    endbp = pos.split('-')[1]
    chromLengths = pd.read_csv('data/hg19_chrom_lengths.txt', sep="\t", encoding='utf-8')
    chromLengths.set_index('sequence',inplace=True)
    if chr == 'X':
        chr = 23
        maxChromLength = chromLengths.loc['chrX', 'length']
    else:
        try:
            chr = int(chr)
            maxChromLength = chromLengths.loc['chr'+str(chr), 'length']
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage("Invalid coordinates input", status_code=410)
    if chr < 1 or chr > 23:
        raise InvalidUsage('Chromosome input must be between 1 and 23', status_code=410)
    elif startbp > endbp:
        raise InvalidUsage('Starting chromosome basepair position is greater than ending basepair position', status_code=410)
    elif startbp > maxChromLength or endbp > maxChromLength:
        raise InvalidUsage('Start or end coordinates are out of range', status_code=410)
    else:
        return chr, startbp, endbp


#####################################
# API Routes
#####################################
class InvalidUsage(Exception):
    status_code = 400

    def __init__(self, message, status_code=None, payload=None):
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        rv = dict(self.payload or ())
        rv['message'] = self.message
        return rv

@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    response = jsonify(error.to_dict())
    response.status_code = error.status_code
    return response

@app.route("/populations")
def get1KGPopulations():
    populations = pd.read_csv('data/populations.tsv', sep='\t')
    return jsonify(populations.to_dict(orient='list'))

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    data = {"success": False}
    if request.method == 'POST':
        if request.files.get('file'):
            # read the file
            file = request.files['file']

            # read the filename
            filename = file.filename

            # create a path to the uploads folder
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)

            file.save(filepath)

            # Load
            print('Loading file')
            gwas_data = pd.read_csv(filepath, sep="\t", encoding='utf-8')

            chromcol = request.form['chrom-col']
            if chromcol=='': chromcol='#CHROM'
            poscol = request.form['pos-col']
            if poscol=='': poscol='BP'
            snpcol = request.form['snp-col']
            if snpcol=='': snpcol='SNP'
            pcol = request.form['pval-col']
            if pcol=='': pcol='P'
            lead_snp = request.form['leadsnp']
            if lead_snp=='': lead_snp = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ]['SNP'])[0]
            regiontext = request.form['locus']
            print('regiontext',regiontext)
            if regiontext == "": raise InvalidUsage("Must enter coordinates", status_code=410)
            print('Parsing region text')
            chr, startbp, endbp = parseRegionText(regiontext)
            print(chr,startbp,endbp)
            print('Subsetting GWAS data to entered region')
            gwas_data = gwas_data.loc[ (gwas_data[chromcol] == chr) & (gwas_data[poscol] >= startbp) & (gwas_data[poscol] <= endbp) ]
            if gwas_data.shape[0] == 0: InvalidUsage('No data found for entered region', status_code=410)
            pops = request.form.getlist('LD-populations')
            if len(pops) == 0: pops = ['CEU','TSI','FIN','GBR','IBS']
            print('Populations:', pops)
            ld_type = request.form['ld-type']
            print(ld_type)
            gtex_tissues = request.form.getlist('GTEx-tissues')
            print('GTEx tissues:',gtex_tissues)
            if len(gtex_tissues) == 0: raise InvalidUsage('Select at least one GTEx tissue', status_code=410)
            gene = request.form['gencodeID']
            if gene=='': raise InvalidUsage('Please enter a gene of interest (eg. ENSG00000174502)')

            # Omit any rows with missing values:
            gwas_data = gwas_data[[snpcol,pcol]]
            gwas_data.dropna(inplace=True)

            # Get LD via API queries to LDlink:
            print('Gathering LD information from LDlink')
            snp_list = list(gwas_data[snpcol])
            ld_df = queryLD(lead_snp, snp_list, pops, ld_type)

            # Get SNP positions via queries to UCSC Genome's MySQL snp151 table:
            print('Gathering HG19 positions from UCSC dbSNP151 database')
            positions = getSNPPositions(snp_list)

            # Make a data dictionary to return as JSON for javascript plot:
            data['snps'] = snp_list
            data['pvalues'] = list(gwas_data[pcol])
            data['lead_snp'] = lead_snp
            data['ld_values'] = ld_df['ld']
            data['positions'] = positions
            data['chr'] = chr
            data['startbp'] = startbp
            data['endbp'] = endbp
            data['ld_populations'] = pops
            data['gtex_tissues'] = gtex_tissues

            # Get GTEx data for the tissues and SNPs selected:
            ensembl_eqtl_base_url = 'http://rest.ensembl.org/eqtl/variant_name/homo_sapiens/'
            gene_query = f'stable_id={gene};'
            query_suffix = 'content-type=application/json'
            for tissue in tqdm(gtex_tissues):
                tissue_query = f'tissue={tissue};'
                data[tissue] = []
                for snp in tqdm(snp_list):
                    querysnp = snp.split(';')[0]
                    snp_query = f'{querysnp}?statistic=p-value;'
                    url = ensembl_eqtl_base_url + snp_query + gene_query + tissue_query + query_suffix
                    response = requests.get(url, headers={ "Content-Type" : "application/json"})
                    if response:
                        eqtl = response.json()
                        try:
                            data[tissue].append(eqtl[0]['minus_log10_p_value'])
                        except:
                            data[tissue].append(np.nan)
                    else:
                        data[tissue].append(np.nan)

            # indicate that the request was a success
            data['success'] = True
            print('Loading a success')

            # Save data in JSON format for plotting
            my_id = uuid.uuid1()
            json.dump(data, open(f'session_data/form_data{my_id}.json', 'w'))

        return redirect("/", code=302)

    return render_template("index.html")
    

if __name__ == "__main__":
    app.run()

