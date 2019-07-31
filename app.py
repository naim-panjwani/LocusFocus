import json
import requests
import pandas as pd
import numpy as np
import os
import math
import string
from tqdm import tqdm
import uuid
from pprint import pprint
import subprocess

import sqlalchemy as sa
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine, inspect, String, Integer

from flask import Flask, request, redirect, url_for, jsonify, render_template, flash
from werkzeug.utils import secure_filename
from flask_sqlalchemy import SQLAlchemy
import pymysql
pymysql.install_as_MySQLdb()
from pymongo import MongoClient
#thepwd = open('pwd.txt').readline().replace('\n', '')

genomicWindowLimit = 2000000
one_sided_SS_window_size = 100000 # (100 kb on either side of the lead SNP)
fileSizeLimit = 500 # in KB

MYDIR = os.path.dirname(__file__)

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'static/upload'
app.config['MAX_CONTENT_LENGTH'] = fileSizeLimit * 1024
ALLOWED_EXTENSIONS = set(['txt', 'tsv'])

# token = ""
# with open('tokens.txt') as f:
#     token = f.read().replace('\n','')

collapsed_genes_df = pd.read_csv(os.path.join(MYDIR, 'data/collapsed_gencode_v19_hg19.gz'), compression='gzip', sep='\t', encoding='utf-8')
# ensg_to_genename_map = pd.read_csv(os.path.join(MYDIR, 'data/gencode_ensg_to_name_map.txt'), sep='\t', encoding='utf-8', header=None)
# ensg_to_genename_map.columns = ['ENSG_name', 'name']
ld_mat_diag_constant = 1e-6

conn = "mongodb://localhost:27017"
client = MongoClient(conn)
db = client.GTEx_V7

####################################
# Helper functions
####################################
def parseRegionText(regiontext):
    regiontext = regiontext.strip().replace(' ','')
    chrom = regiontext.split(':')[0].replace('chr','').replace('Chr','')
    pos = regiontext.split(':')[1]
    startbp = pos.split('-')[0].replace(',','')
    endbp = pos.split('-')[1].replace(',','')
    chromLengths = pd.read_csv(os.path.join(MYDIR, 'data/hg19_chrom_lengths.txt'), sep="\t", encoding='utf-8')
    chromLengths.set_index('sequence',inplace=True)
    if chrom == 'X':
        chrom = 23
        maxChromLength = chromLengths.loc['chrX', 'length']
        try:
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage("Invalid coordinates input", status_code=410)
    else:
        try:
            chrom = int(chrom)
            maxChromLength = chromLengths.loc['chr'+str(chrom), 'length']
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage("Invalid coordinates input", status_code=410)
    if chrom < 1 or chrom > 23:
        raise InvalidUsage('Chromosome input must be between 1 and 23', status_code=410)
    elif startbp > endbp:
        raise InvalidUsage('Starting chromosome basepair position is greater than ending basepair position', status_code=410)
    elif startbp > maxChromLength or endbp > maxChromLength:
        raise InvalidUsage('Start or end coordinates are out of range', status_code=410)
    else:
        return chrom, startbp, endbp

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def writeList(alist, filename):
    with open(filename, 'w') as f:
        for item in alist:
            f.write("%s\n" % item)

def writeMat(aMat, filename):
    aMat = np.matrix(aMat)
    with open(filename, 'w') as f:
        for row in np.arange(aMat.shape[0]):
            for col in np.arange(aMat.shape[1] - 1):
                f.write("%s\t" % str(aMat[row,col]))
            f.write("%s\n" % str(aMat[row,-1]))

####################################
# LD Calculation from 1KG using PLINK
####################################

def plink_ldmat(pop, chrom, snp_positions, outfilename):
    # positions must be in hg19 coordinates
    if chrom == 'X': chrom = 23
    try:
        chrom = int(chrom)
    except:
        raise InvalidUsage("Invalid chromosome", status_code=410)
    if chrom not in np.arange(1,24):
        raise InvalidUsage("Invalid chromosome", status_code=410)
    plink_filepath = ""
    if chrom == 23:
        plink_filepath = os.path.join(MYDIR, "data", pop, "chrX")
    else:
        plink_filepath = os.path.join(MYDIR, "data", pop, f"chr{chrom}")
    # make snps file to extract:
    snps = [f"chr{str(int(chrom))}:{str(int(position))}" for position in snp_positions]
    writeList(snps, outfilename + "_snps.txt")
    plink_path = subprocess.run(args=["which","plink"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
    plinkrun = subprocess.run(args=[
        "./plink", '--bfile', plink_filepath
        , "--chr", str(chrom)
        , "--extract", outfilename + "_snps.txt"
        , "--from-bp", str(min(snp_positions))
        , "--to-bp", str(max(snp_positions))
        , "--r2", "square"
        , "--make-bed"
        , "--threads", "1"
        , "--out", outfilename
        ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if plinkrun.returncode != 0:
        raise InvalidUsage(plinkrun.stdout.decode('utf-8'), status_code=410)
    ld_snps = list(pd.read_csv(outfilename + ".bim", sep="\t", header=None).iloc[:,1])
    ldmat = np.matrix(pd.read_csv(outfilename + ".ld", sep="\t", header=None))
    return ld_snps, ldmat

def plink_ld_pairwise(lead_snp_position, pop, chrom, snp_positions, outfilename):
    # positions must be in hg19 coordinates
    # returns NaN for SNPs not in 1KG LD file; preserves order of input snp_positions
    if chrom == 'X': chrom = 23
    try:
        chrom = int(chrom)
    except:
        raise InvalidUsage("Invalid chromosome", status_code=410)
    if chrom not in np.arange(1,24):
        raise InvalidUsage("Invalid chromosome", status_code=410)
    plink_filepath = ""
    if chrom == 23:
        plink_filepath = os.path.join(MYDIR, "data", pop, "chrX")
    else:
        plink_filepath = os.path.join(MYDIR, "data", pop, f"chr{chrom}")
    # make snps file to extract:
    snps = [f"chr{str(int(chrom))}:{str(int(position))}" for position in snp_positions]
    writeList(snps, outfilename + "_snps.txt")
    plink_path = subprocess.run(args=["which","plink"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
    plinkrun = subprocess.run(args=[
        "./plink", '--bfile', plink_filepath
        , "--chr", str(chrom)
        , "--extract", outfilename + "_snps.txt"
        , "--from-bp", str(min(snp_positions))
        , "--to-bp", str(max(snp_positions))
        , "--ld-snp", f"chr{str(int(chrom))}:{str(int(lead_snp_position))}"
        , "--r2"
        , "--ld-window-r2", "0"
        , "--ld-window", "999999"
        , "--ld-window-kb", "200000"
        , "--make-bed"
        , "--threads", "1"
        , "--out", outfilename
        ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if plinkrun.returncode != 0:
        raise InvalidUsage(plinkrun.stdout.decode('utf-8'), status_code=410)
    ld_results = pd.read_csv(outfilename + ".ld", delim_whitespace=True)
    available_r2_positions = ld_results[['BP_B', 'R2']]
    pos_df = pd.DataFrame({'pos': snp_positions})
    merged_df = pd.merge(pos_df, available_r2_positions, how='left', left_on="pos", right_on="BP_B", sort=False)[['pos', 'R2']]
    merged_df.fillna(-1, inplace=True)
    return merged_df

def get_gtex_v7(tissue, gene_id):
    tissue = tissue.title().replace(' ','_')
    gene_id = gene_id.upper()
    ensg_name = ""
    if tissue not in db.list_collection_names():
        raise InvalidUsage('Tissue not found', status_code=410)
    collection = db[tissue]
    if gene_id.startswith('ENSG'):
        i = list(collapsed_genes_df['ENSG_name']).index(gene_id)
        ensg_name = list(collapsed_genes_df['ENSG_name'])[i]
    elif gene_id in list(collapsed_genes_df['name']):
        i = list(collapsed_genes_df['name']).index(gene_id)
        ensg_name = list(collapsed_genes_df['ENSG_name'])[i]
    else:
        raise InvalidUsage(f'Gene name {gene_id} not found', status_code=410)
    results = list(collection.find({'gene_id': ensg_name}))
    response = []
    try:
        response = results[0]['eqtl_variants']
    except:
        return pd.DataFrame([{'error': f'No eQTL data for {gene_id} in {tissue}'}])
    results_df = pd.DataFrame(response)
    chrom = int(list(results_df['variant_id'])[0].split('_')[0].replace('X','23'))
    positions = [ int(x.split('_')[1]) for x in list(results_df['variant_id']) ]
    variants_query = db.variant_table.aggregate([
        { '$match': { '$and': [ 
            { 'chr': chrom }, 
            { 'variant_pos': { '$gte': min(positions), '$lte': max(positions) } } 
            ] 
            } 
        }
    ])
    variants_df = pd.DataFrame(list(variants_query)).drop(['_id'], axis=1)
    x = pd.merge(results_df, variants_df, on='variant_id')
    x.rename(columns={'rs_id_dbSNP147_GRCh37p13': 'rs_id'}, inplace=True)
    return x


def get_gtex_data(tissue, gene, snp_list, raiseErrors = False):
    gtex_data = []
    rsids = True
    if snp_list[0].startswith('rs'):
        rsids = True
    elif snp_list[0].endswith('_b37'):
        rsids = False
    else:
        raise InvalidUsage('Variant naming format not supported; ensure all are rs ID\'s or formatted as chrom_pos_ref_alt_b37 eg. 1_205720483_G_A_b37')
    ensg_gene = gene
    if gene in list(collapsed_genes_df['name']):
        ensg_gene = collapsed_genes_df['ENSG_name'][list(collapsed_genes_df['name']).index(gene)]
    if gene in list(collapsed_genes_df['ENSG_name']):
        gene = collapsed_genes_df['name'][list(collapsed_genes_df['ENSG_name']).index(gene)]
    print(f'Gathering eQTL data for {gene} ({ensg_gene}) in {tissue}')
    response_df = get_gtex_v7(tissue, gene)
    if 'error' not in response_df.columns:
        eqtl = response_df
        if rsids:
            snp_df = pd.DataFrame(snp_list, columns=['rs_id'])
            gtex_data = snp_df.reset_index().merge(eqtl, on='rs_id', how='left', sort=False).sort_values('index')
        else:
            snp_df = pd.DataFrame(snp_list, columns=['variant_id'])
            gtex_data = snp_df.reset_index().merge(eqtl, on='variant_id', how='left', sort=False).sort_values('index')
    else:
        try:
            error_message = list(response_df['error'])[0]
            if raiseErrors:
                raise InvalidUsage(error_message, status_code=410)
        except:
            if raiseErrors:
                raise InvalidUsage("No response for tissue " + tissue.replace("_"," ") + " and gene " + gene + " ( " + ensg_gene + " )", status_code=410)
    return gtex_data


def get_gtex_data_pvalues(eqtl_data, snp_list):
    rsids = True
    if snp_list[0].startswith('rs'):
        rsids = True
    elif snp_list[0].endswith('_b37'):
        rsids = False
    else:
        raise InvalidUsage('Variant naming format not supported; ensure all are rs ID\'s or formatted as chrom_pos_ref_alt_b37 eg. 1_205720483_G_A_b37')    
    if rsids:
        gtex_data = pd.merge(eqtl_data, pd.DataFrame(snp_list, columns=['rs_id']), on='rs_id', how='right')
    else:
        gtex_data = pd.merge(eqtl_data, pd.DataFrame(snp_list, columns=['variant_id']), on='variant_id', how='right')
    return list(gtex_data['pval'])

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
    populations = pd.read_csv(os.path.join(MYDIR, 'data/populations.tsv'), sep='\t')
    return jsonify(populations.to_dict(orient='list'))

@app.route("/genenames")
def getGeneNames():
    return jsonify(list(collapsed_genes_df['name']))

@app.route("/gtex_v7/tissues_list")
def list_tissues():
    tissues = list(db.list_collection_names())
    tissues.remove('variant_table')
    return jsonify(tissues)

@app.route("/gtex_v7/<tissue>/<gene_id>")
def get_gtex(tissue, gene_id):
    x = get_gtex_v7(tissue, gene_id)
    return jsonify(x.to_dict(orient='records'))

@app.route("/gtex_v7/<tissue>/<gene_id>/<variant>")
def get_gtex_variant(tissue, gene_id, variant):
    x = get_gtex_v7(tissue, gene_id)
    response_df = x
    result = []
    if variant.startswith("rs"):
        result = response_df.loc[ response_df['rs_id'] == variant ]
    elif variant.endswith("_b37"):
        result = response_df.loc[ response_df['variant_id'] == variant ]
    else:
        raise InvalidUsage(f'variant name {variant} not found', status_code=410)
    if result.shape[0] == 0:
        raise InvalidUsage(f'variant name {variant} not found', status_code=410)
    return jsonify(result.to_dict(orient='records'))

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    data = {"success": False}
    ldmat_file_supplied = False
    if request.method == 'POST':
        if request.files.get('file'):
            # read the file
            file = request.files['file']
            # read the filename
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename) # more secure
                # create a path to the uploads folder
                filepath = os.path.join(MYDIR, app.config['UPLOAD_FOLDER'], filename)
                file.save(filepath)
            else:
                raise InvalidUsage('GWAS summary statistics file type not allowed', status_code=410)
            try:
                ldmat_file = request.files['ldmat']
                # read the filename
                if ldmat_file and allowed_file(ldmat_file.filename):
                    ldmat_filename = secure_filename(ldmat_file.filename) # more secure
                    # create a path to the uploads folder
                    ldmat_filepath = os.path.join(MYDIR, app.config['UPLOAD_FOLDER'], ldmat_filename)
                    file.save(ldmat_filepath)
                    ldmat_file_supplied = True
                else:
                    raise InvalidUsage('LD matrix file type not allowed', status_code=410)
            except:
                pass
            # Load
            my_session_id = uuid.uuid4()
            print(f'Session ID: {my_session_id}')
            print('Loading file')
            gwas_data = pd.read_csv(filepath, sep="\t", encoding='utf-8')
            chromcol = request.form['chrom-col']
            if chromcol=='': chromcol='#CHROM'
            if chromcol not in gwas_data.columns:
                raise InvalidUsage('Chromosome column not found', status_code=410)
            poscol = request.form['pos-col']
            if poscol=='': poscol='BP'
            if poscol not in gwas_data.columns:
                raise InvalidUsage('Basepair position column not found', status_code=410)
            snpcol = request.form['snp-col']
            if snpcol=='': snpcol='SNP'
            if snpcol not in gwas_data.columns:
                raise InvalidUsage('SNP column not found', status_code=410)
            pcol = request.form['pval-col']
            if pcol=='': pcol='P'
            if pcol not in gwas_data.columns:
                raise InvalidUsage('Basepair position column not found', status_code=410)
            lead_snp = request.form['leadsnp']
            snp_list = list(gwas_data[snpcol])
            snp_list = [asnp.split(';')[0] for asnp in snp_list] # cleaning up the SNP names a bit
            if lead_snp=='': lead_snp = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ][snpcol])[0].split(';')[0]
            if lead_snp not in snp_list:
                raise InvalidUsage('Lead SNP not found', status_code=410)
            lead_snp_position = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ][poscol])[0]
            regiontext = request.form['locus']
            print('regiontext',regiontext)
            if regiontext == "": regiontext = "1:205500000-206000000"
            print('Parsing region text')
            chrom, startbp, endbp = parseRegionText(regiontext)
            if (endbp - startbp) > genomicWindowLimit:
                raise InvalidUsage(f'Entered region size is larger than {genomicWindowLimit} bp', status_code=410)
            print(chrom,startbp,endbp)
            print('Subsetting GWAS data to entered region')            
            if chrom == 23 and list(gwas_data[chromcol])[0] == "X":
                gwas_data = gwas_data.loc[ (gwas_data[chromcol] == "X") & (gwas_data[poscol] >= startbp) & (gwas_data[poscol] <= endbp) ]
            else:
                gwas_data = gwas_data.loc[ (gwas_data[chromcol] == chrom) & (gwas_data[poscol] >= startbp) & (gwas_data[poscol] <= endbp) ]
            if gwas_data.shape[0] == 0:
                raise InvalidUsage('No data found for entered region', status_code=410)
            gwas_data.sort_values(by=[ poscol ], inplace=True)
            #pops = request.form.getlist('LD-populations')
            #if len(pops) == 0: pops = ['CEU','TSI','FIN','GBR','IBS']
            pops = request.form['LD-populations']
            if len(pops) == 0: pops = 'EUR'
            print('Populations:', pops)
            # ld_type = request.form['ld-type']
            # print(ld_type)
            gtex_tissues = request.form.getlist('GTEx-tissues')
            print('GTEx tissues:',gtex_tissues)
            if len(gtex_tissues) == 0: raise InvalidUsage('Select at least one GTEx tissue', status_code=410)            
            gene = request.form['gencodeID']
            if gene=='': 
                gene='ENSG00000174502.14'
            elif not (str(gene).upper().startswith('ENSG')):
                try:
                    gene = str(list(collapsed_genes_df.loc[ collapsed_genes_df['name'] == str(gene).upper() ]['ENSG_name'])[0])
                except:
                    raise InvalidUsage('Gene name not recognized', status_code=410)
            # Omit any rows with missing values:
            gwas_data = gwas_data[[ chromcol, poscol, snpcol, pcol ]]
            gwas_data.dropna(inplace=True)
            snp_list = list(gwas_data[snpcol])
            snp_list = [asnp.split(';')[0] for asnp in snp_list] # cleaning up the SNP names a bit
            # Get LD:
            print('Calculating pairwise LD using PLINK')
            positions = list(gwas_data[poscol])
            #ld_df = queryLD(lead_snp, snp_list, pops, ld_type)
            ld_df = plink_ld_pairwise(lead_snp_position, pops, chrom, positions, os.path.join(MYDIR, "static", "session_data", f"ld-{my_session_id}"))
            data = {}
            data['snps'] = snp_list
            data['pvalues'] = list(gwas_data[pcol])
            data['lead_snp'] = lead_snp
            data['ld_values'] = list(ld_df['R2'])
            data['positions'] = positions
            data['chrom'] = chrom
            data['startbp'] = startbp
            data['endbp'] = endbp
            data['ld_populations'] = pops
            data['gtex_tissues'] = gtex_tissues
            # Get GTEx data for the tissues and SNPs selected:
            print('Gathering GTEx data')
            gtex_data = {}
            for tissue in tqdm(gtex_tissues):
                gtex_data[tissue] = get_gtex_data(tissue, gene, snp_list, raiseErrors=True).to_dict(orient='records')
            data.update(gtex_data)
            
            # Determine the region to calculate the Simple Sum (SS):
            SS_start = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ][poscol])[0] - one_sided_SS_window_size
            SS_end = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ][poscol])[0] + one_sided_SS_window_size
            data['SS_region'] = [SS_start, SS_end]
            
            # Obtain any genes to be plotted in the region:
            print('Summarizing genes to be plotted in this region')
            genes_to_draw = collapsed_genes_df.loc[ (collapsed_genes_df['chrom'] == ('chr' + str(chrom).replace('23','X'))) &
                                                    ( ((collapsed_genes_df['txStart'] >= startbp) & (collapsed_genes_df['txStart'] <= endbp)) | 
                                                      ((collapsed_genes_df['txEnd'] >= startbp  ) & (collapsed_genes_df['txEnd'] <= endbp  )) | 
                                                      ((collapsed_genes_df['txStart'] <= startbp) & (collapsed_genes_df['txEnd'] >= endbp  )) )]
            genes_data = []
            for i in np.arange(genes_to_draw.shape[0]):
                genes_data.append({
                    'name': list(genes_to_draw['name'])[i]
                    ,'txStart': list(genes_to_draw['txStart'])[i]
                    ,'txEnd': list(genes_to_draw['txEnd'])[i]
                    ,'exonStarts': [int(bp) for bp in list(genes_to_draw['exonStarts'])[i].split(',')]
                    ,'exonEnds': [int(bp) for bp in list(genes_to_draw['exonEnds'])[i].split(',')]
                })

            # Indicate that the request was a success
            data['success'] = True
            print('Loading a success')

            # Save data in JSON format for plotting
            sessionfile = f'session_data/form_data-{my_session_id}.json'
            sessionfilepath = os.path.join(MYDIR, 'static', sessionfile)
            json.dump(data, open(sessionfilepath, 'w'))
            genes_sessionfile = f'session_data/genes_data-{my_session_id}.json'
            genes_sessionfilepath = os.path.join(MYDIR, 'static', genes_sessionfile) 
            json.dump(genes_data, open(genes_sessionfilepath, 'w'))

            # # Getting Simple Sum P-values
            # 2. Subset the region (step 1 was determining the region to do the SS calculation on - see above SS_start and SS_end variables):
            print('SS_start: ' + str(SS_start))
            print('SS_end:' + str(SS_end))
            chromList = [('chr' + str(chrom).replace('23','X')), str(chrom).replace('23','X')]
            gwas_chrom_col = pd.Series([str(x) for x in list(gwas_data[chromcol])])
            SS_chrom_bool = [x for x in gwas_chrom_col.isin(chromList) if x == True]
            SS_gwas_data = gwas_data.loc[ SS_chrom_bool & (gwas_data[poscol] >= SS_start) & (gwas_data[poscol] <= SS_end) ]
            #print(gwas_data.shape)
            #print(SS_gwas_data.shape)
            #print(SS_gwas_data)
            if SS_gwas_data.shape[0] == 0: InvalidUsage('No data points found for entered Simple Sum region', status_code=410)
            PvaluesMat = [list(SS_gwas_data[pcol])]
            SS_snp_list = list(SS_gwas_data[snpcol])
            SS_snp_list = [asnp.split(';')[0] for asnp in SS_snp_list] # cleaning up the SNP names a bit
            SS_positions = list(SS_gwas_data[poscol])
            # Extra file written:
            gwas_df = pd.DataFrame({
                'Position': SS_positions,
                'SNP': SS_snp_list,
                'P': list(SS_gwas_data[pcol])
            })
            gwas_df.to_csv(os.path.join(MYDIR, 'static', f'session_data/gwas_df-{my_session_id}.txt'), index=False, encoding='utf-8', sep="\t")
            # 3. Determine the genes to query
            query_genes = list(genes_to_draw['name'])
            # 4. Query and extract the eQTL p-values for all tissues & genes from GTEx
            for tissue in gtex_tissues:
                for gene in query_genes:
                    gtex_eqtl_df = get_gtex_data(tissue, gene, SS_snp_list)
                    print(len(gtex_eqtl_df))
                    if len(gtex_eqtl_df) > 0:
                        pvalues = list(gtex_eqtl_df['pval'])
                    else:
                        pvalues = np.repeat(np.nan, len(SS_snp_list))
                    PvaluesMat.append(pvalues)
                    print(f'tissue: {tissue}, gene: {gene}, len(pvalues): {len(pvalues)}')
                    print(f'len(SS_positions): {len(SS_positions)}, len(SS_snp_list): {len(SS_snp_list)}')
                    # Extra files written:
                    eqtl_df = pd.DataFrame({
                        'Position': SS_positions,
                        'SNP': SS_snp_list,
                        'P': pvalues
                    })
                    eqtl_df.to_csv(os.path.join(MYDIR, 'static', f'session_data/eqtl_df-{tissue}-{gene}-{my_session_id}.txt'), index=False, encoding='utf-8', sep="\t")
            # 5. Get the LD matrix via PLINK subprocess call:
            plink_outfilename = f'session_data/ld-{my_session_id}'
            plink_outfilepath = os.path.join(MYDIR, 'static', plink_outfilename)
            if ldmat_file_supplied:
                ld_mat_snps = [f'chr{chrom}:{SS_pos}' for SS_pos in SS_positions]
                ld_mat = np.matrix(pd.read_csv(ldmat_filepath, sep="\t", encoding='utf-8', header=None))
                ld_mat_positions = [int(snp.split(":")[1]) for snp in ld_mat_snps]
            else:
                ld_mat_snps, ld_mat = plink_ldmat(pops, chrom, SS_positions, plink_outfilepath)
                ld_mat_positions = [int(snp.split(":")[1]) for snp in ld_mat_snps]
            np.fill_diagonal(ld_mat, np.diag(ld_mat) + ld_mat_diag_constant)
            # 6. Shrink the P-values matrix to include only the SNPs available in the LD matrix:
            PvaluesMat = np.matrix(PvaluesMat)
            Pmat_indices = [i for i, e in enumerate(SS_positions) if e in ld_mat_positions]
            PvaluesMat = PvaluesMat[:, Pmat_indices]
            # 7. Write the p-values and LD matrix into session_data
            Pvalues_file = f'session_data/Pvalues-{my_session_id}.txt'
            ldmatrix_file = f'session_data/ldmat-{my_session_id}.txt'
            Pvalues_filepath = os.path.join(MYDIR, 'static', Pvalues_file)
            ldmatrix_filepath = os.path.join(MYDIR, 'static', ldmatrix_file)
            writeMat(PvaluesMat, Pvalues_filepath)
            writeMat(ld_mat, ldmatrix_filepath)
            #### Extra files written for LD matrix:
            writeList(ld_mat_snps, os.path.join(MYDIR,'static', f'session_data/ldmat_snps-{my_session_id}.txt'))
            writeList(ld_mat_positions, os.path.join(MYDIR,'static', f'session_data/ldmat_positions-{my_session_id}.txt'))
            ####
            Rscript_code_path = os.path.join(MYDIR, 'getSimpleSumStats.R')
            Rscript_path = subprocess.run(args=["which","Rscript"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
            RscriptRun = subprocess.run(args=[Rscript_path, Rscript_code_path, Pvalues_filepath, ldmatrix_filepath], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            if RscriptRun.returncode != 0:
                raise InvalidUsage(RscriptRun.stdout, status_code=410)
            SSPvalues = RscriptRun.stdout.replace('\n',' ').split(' ')
            SSPvalues = [float(SSP) for SSP in SSPvalues if SSP!='']
            for i in np.arange(len(SSPvalues)):
                if SSPvalues[i] != -1:
                    SSPvalues[i] = np.format_float_scientific((-np.log10(SSPvalues[i])), precision=2)
            SSPvaluesMat = np.array(SSPvalues).reshape(len(gtex_tissues), len(query_genes))
            SSPvalues_dict = {
                'Genes': query_genes
                ,'Tissues': gtex_tissues
                ,'SSPvalues': SSPvaluesMat.tolist()
            }
            SSPvalues_file = f'session_data/SSPvalues-{my_session_id}.json'
            SSPvalues_filepath = os.path.join(MYDIR, 'static', SSPvalues_file)
            json.dump(SSPvalues_dict, open(SSPvalues_filepath, 'w'))
            return render_template("plot.html", sessionfile = sessionfile, genesfile = genes_sessionfile, SSPvalues_file = SSPvalues_file)
        return render_template("invalid_input.html")
    return render_template("index.html")

if __name__ == "__main__":
    app.run()

