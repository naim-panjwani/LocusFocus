import json
#import requests
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import uuid
import subprocess
from datetime import datetime
from bs4 import BeautifulSoup as bs
import re
import pysam
import secrets

from flask import Flask, request, redirect, url_for, jsonify, render_template, flash, send_file
from werkzeug.utils import secure_filename
from flask_sitemap import Sitemap
from flask_uploads import UploadSet, configure_uploads, DATA

from pymongo import MongoClient

#import getSimpleSumStats

genomicWindowLimit = 2000000
one_sided_SS_window_size = 100000 # (100 kb on either side of the lead SNP)
fileSizeLimit = 100 * 1024 * 1024 # in Bytes

MYDIR = os.path.dirname(__file__)
APP_STATIC = os.path.join(MYDIR, 'static')

##################
# Default settings
##################

default_region = "1:205500000-206000000"
default_chromname = "#CHROM"
default_posname = "POS"
default_snpname = "ID"
default_refname = "REF"
default_altname = "ALT"
default_pname = "P"
default_betaname = "BETA"
default_stderrname = "SE"
default_nname = "N"
default_mafname = "MAF"

# Default column names for secondary datasets:
CHROM = 'CHROM'
BP = 'BP'
SNP = 'SNP'
P = 'P'

coloc2colnames = ['CHR','POS','SNPID','A2','A1','BETA','SE','PVAL','MAF', 'N']
coloc2eqtlcolnames = coloc2colnames + ['ProbeID']
coloc2gwascolnames = coloc2colnames + ['type']

################
################

app = Flask(__name__)
ext = Sitemap(app=app)
app.config['UPLOAD_FOLDER'] = os.path.join(MYDIR, 'static/upload/')
app.config['UPLOADED_FILES_DEST'] = os.path.join(MYDIR, 'static/upload/')
app.config['MAX_CONTENT_LENGTH'] = fileSizeLimit
ALLOWED_EXTENSIONS = set(['txt', 'tsv', 'ld', 'html'])
app.config['UPLOADED_FILES_ALLOW'] = ALLOWED_EXTENSIONS
app.secret_key = secrets.mysecret
files = UploadSet('files', DATA)
configure_uploads(app, files)

collapsed_genes_df_hg19 = pd.read_csv(os.path.join(MYDIR, 'data/collapsed_gencode_v19_hg19.gz'), compression='gzip', sep='\t', encoding='utf-8')
collapsed_genes_df_hg38 = pd.read_csv(os.path.join(MYDIR, 'data/collapsed_gencode_v26_hg38.gz'), compression='gzip', sep='\t', encoding='utf-8')

collapsed_genes_df = collapsed_genes_df_hg19 # For now
ld_mat_diag_constant = 1e-6

conn = "mongodb://localhost:27017"
client = MongoClient(conn)
db = client.GTEx_V7 # For now

available_gtex_versions = ["V7", "V8"]
valid_populations = ["EUR", "AFR","EAS", "SAS", "AMR", "ASN"]

####################################
# Helper functions
####################################
def parseRegionText(regiontext, build):
    if build not in ['hg19', 'hg38']:
        raise InvalidUsage(f'Unrecognized build: {build}', status_code=410)
    regiontext = regiontext.strip().replace(' ','').replace(',','').replace('chr','')
    if not re.search("^\d+:\d+-\d+$", regiontext.replace('X','23').replace('x','23')):
       raise InvalidUsage(f'Invalid coordinate format. {regiontext} e.g. 1:205,000,000-206,000,000', status_code=410)
    chrom = regiontext.split(':')[0].lower().replace('chr','').upper()
    pos = regiontext.split(':')[1]
    startbp = pos.split('-')[0].replace(',','')
    endbp = pos.split('-')[1].replace(',','')
    chromLengths = pd.read_csv(os.path.join(MYDIR, 'data', build + '_chrom_lengths.txt'), sep="\t", encoding='utf-8')
    chromLengths.set_index('sequence',inplace=True)
    if chrom in ['X','x'] or chrom == '23':
        chrom = 23
        maxChromLength = chromLengths.loc['chrX', 'length']
        try:
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage(f"Invalid coordinates input: {regiontext}", status_code=410)
    else:
        try:
            chrom = int(chrom)
            if chrom == 23:
                maxChromLength = chromLengths.loc['chrX', 'length']
            else:
                maxChromLength = chromLengths.loc['chr'+str(chrom), 'length']
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage(f"Invalid coordinates input {regiontext}", status_code=410)
    if chrom < 1 or chrom > 23:
        raise InvalidUsage('Chromosome input must be between 1 and 23', status_code=410)
    elif startbp > endbp:
        raise InvalidUsage('Starting chromosome basepair position is greater than ending basepair position', status_code=410)
    elif startbp > maxChromLength or endbp > maxChromLength:
        raise InvalidUsage('Start or end coordinates are out of range', status_code=410)
    elif (endbp - startbp) > genomicWindowLimit:
        raise InvalidUsage(f'Entered region size is larger than {genomicWindowLimit/10**6} Mbp', status_code=410)
    else:
        return chrom, startbp, endbp

def allowed_file(filenames):
    if type(filenames) == type('str'):
        return '.' in filenames and filenames.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
    for filename in filenames:
        if not ('.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS):
            return False
    return True

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


def genenames(genename, build):
    # Given either ENSG gene name or HUGO gene name, returns both HUGO and ENSG names
    ensg_gene = genename
    if build.lower() in ["hg19","grch37"]:
        collapsed_genes_df = collapsed_genes_df_hg19
    elif build.lower() in ["hg38", "grch38"]:
        collapsed_genes_df = collapsed_genes_df_hg38
    if genename in list(collapsed_genes_df['name']):
        ensg_gene = collapsed_genes_df['ENSG_name'][list(collapsed_genes_df['name']).index(genename)]
    if genename in list(collapsed_genes_df['ENSG_name']):
        genename = collapsed_genes_df['name'][list(collapsed_genes_df['ENSG_name']).index(genename)]
    return genename, ensg_gene

def classify_files(filenames):
    gwas_filepath = ''
    ldmat_filepath = ''
    html_filepath = ''
    extensions = []
    for file in filenames:
        filename = secure_filename(file.filename)
        extension = filename.split('.')[-1]
        if extension not in extensions:
            if extension in ['txt', 'tsv']:
                extensions.extend(['txt','tsv'])
            else:
                extensions.append(extension)
        else:
            raise InvalidUsage('Please upload up to 3 different file types as described', status_code=410)
        if extension in ['txt', 'tsv']:
            gwas_filepath = os.path.join(MYDIR, app.config['UPLOAD_FOLDER'], filename)
        elif extension in ['ld']:
            ldmat_filepath = os.path.join(MYDIR, app.config['UPLOAD_FOLDER'], filename)
        elif extension in ['html']:
            html_filepath = os.path.join(MYDIR, app.config['UPLOAD_FOLDER'], filename)
    return gwas_filepath, ldmat_filepath, html_filepath


def isSorted(l):
    # l is a list
    # returns True if l is sorted, False otherwise
    return all(l[i] <= l[i+1] for i in range(len(l)-1))

def Xto23(l):
    newl = []
    for x in l:
        if str(str(x).strip().lower().replace('chr','').upper()) == "X":
            newl.append(23)
        elif str(str(x).strip().lower().replace('chr','')) in [str(i) for i in list(np.arange(1,24))]:
            newl.append(int(str(x).strip().lower().replace('chr','')))
        else:
            raise InvalidUsage('Chromosome unrecognized', status_code=410)
    return newl

def verifycol(formname, defaultname, filecolnames, error_message_):
    """
    Checks if the user-entered column name (formname)
    (or the default column name if no column name was entered - defaultname)
    can be found in the dataset column names (ie. filecolnames list).
    If not, the error_message_ is output and program halted with 410 status
    """
    theformname = formname
    if formname=='': theformname=str(defaultname)
    if theformname not in filecolnames:
        raise InvalidUsage(error_message_, status_code=410)
    return theformname


def buildSNPlist(df, chromcol, poscol, refcol, altcol, build):
    snplist = []
    if build.lower() in ["hg38","grch38"]:
        build = 'b38'
    else:
        build = 'b37'
    for i in np.arange(df.shape[0]):
        chrom = list(df[chromcol])[i]
        pos = list(df[poscol])[i]
        ref = list(df[refcol])[i]
        alt = list(df[altcol])[i]
        try:
            snplist.append(str(chrom)+"_"+str(pos)+"_"+str(ref)+"_"+str(alt)+"_"+str(build))
        except:
            raise InvalidUsage(f'Could not convert marker at row {str(i)}')
    return snplist


def fetchSNV(chrom, bp, ref, build):
    variantid = '.'
    
    if ref is None or ref=='.':
        ref=''
    
    # Ensure valid region:
    try:
        regiontxt = str(chrom) + ":" + str(bp) + "-" + str(int(bp)+1)
    except:
        raise InvalidUsage(f'Invalid input for str(chrom):str(bp)')
    chrom, startbp, endbp = parseRegionText(regiontxt, build)
    chrom = str(chrom).replace('chr','').replace('23',"X")
    
    # Load dbSNP151 SNP names from region indicated
    dbsnp_filepath = ''
    if build.lower() in ["hg38", "grch38"]:
        suffix = 'b38'
        dbsnp_filepath = os.path.join(MYDIR, 'data', 'dbSNP151', 'GRCh38p7', 'All_20180418.vcf.gz')
    else:
        suffix = 'b37'
        dbsnp_filepath = os.path.join(MYDIR, 'data', 'dbSNP151', 'GRCh37p13', 'All_20180423.vcf.gz')

    # Load variant info from dbSNP151
    tbx = pysam.TabixFile(dbsnp_filepath)
    varlist = []
    for row in tbx.fetch(str(chrom), bp-1, bp):
        rowlist = str(row).split('\t')
        chromi = rowlist[0].replace('chr','')
        posi = rowlist[1]
        idi = rowlist[2]
        refi = rowlist[3]
        alti = rowlist[4]
        varstr = '_'.join([chromi, posi, refi, alti, suffix])
        varlist.append(varstr)

    # Check if there is a match to an SNV with the provided info
    if len(varlist) == 1:
        variantid = varstr
    elif len(varlist) > 1 and ref != '':
        for v in varlist:
            if v.split('_')[2] == ref:
                variantid = v
                break
    return variantid


def mapAndCleanSNPs(variantlist, regiontxt, build):
    """
    Input: Variant names in any of these formats: rsid, chrom_pos_ref_alt, chrom:pos_ref_alt, chrom:pos_ref_alt_b37/b38 
    Output: chrom_pos_ref_alt_b37/b38 variant ID format.
    In the case of multi-allelic variants (e.g. rs2211330(T/A,C)), formats such as 1_205001063_T_A,C_b37 are accepted
    If variant ID format is chr:pos, and the chr:pos has a unique biallelic SNV, then it will be assigned that variant
    """
    
    # Ensure valid region:
    chrom, startbp, endbp = parseRegionText(regiontxt, build)
    chrom = str(chrom).replace('23',"X")

    # Load dbSNP151 SNP names from region indicated
    dbsnp_filepath = ''
    if build.lower() in ["hg38", "grch38"]:
        suffix = 'b38'
        dbsnp_filepath = os.path.join(MYDIR, 'data', 'dbSNP151', 'GRCh38p7', 'All_20180418.vcf.gz')
    else:
        suffix = 'b37'
        dbsnp_filepath = os.path.join(MYDIR, 'data', 'dbSNP151', 'GRCh37p13', 'All_20180423.vcf.gz')
    
    
    # Load dbSNP file
    #delayeddf = delayed(pd.read_csv)(dbsnp_filepath,skiprows=getNumHeaderLines(dbsnp_filepath),sep='\t')
    #dbsnp = dd.from_delayed(delayeddf)
    tbx = pysam.TabixFile(dbsnp_filepath)
    print('Compiling list of known variants in the region from dbSNP151')
    chromcol = []
    poscol = []
    idcol = []
    refcol = []
    altcol = []
    variantid = [] # in chr_pos_ref_alt_build format
    rsids = dict({}) # a multi-allelic variant rsid (key) can be represented in several variantid formats (values)
    for row in tbx.fetch(str(chrom), startbp, endbp):
        rowlist = str(row).split('\t')
        chromi = rowlist[0].replace('chr','')
        posi = rowlist[1]
        idi = rowlist[2]
        refi = rowlist[3]
        alti = rowlist[4]
        varstr = '_'.join([chromi, posi, refi, alti, suffix])
        chromcol.append(chromi)
        poscol.append(posi)
        idcol.append(idi)
        refcol.append(refi)
        altcol.append(alti)
        variantid.append(varstr)
        rsids[idi] = [varstr]
        altalleles = alti.split(',') # could have more than one alt allele (multi-allelic)
        if len(altalleles)>1:
            varstr = '_'.join([chromi, posi, refi, altalleles[0], suffix])
            rsids[idi].append(varstr)
            for i in np.arange(len(altalleles)-1):
                varstr = '_'.join([chromi, posi, refi, altalleles[i+1], suffix])
                rsids[idi].append(varstr)
    
    print('Cleaning and mapping list of variants')
    variantlist = [asnp.split(';')[0].replace(':','_').replace('.','') for asnp in variantlist] # cleaning up the SNP names a bit
    stdvariantlist = []
    for variant in variantlist:
        variantstr = variant.replace('chr','')
        if variantstr.startswith('rs'):
            try:
                stdvariantlist.append(rsids[variantstr][0])
            except:
                stdvariantlist.append('.')
        elif re.search("^\d+_\d+_[A,T,G,C]+_[A,T,C,G]+,*", variantstr.replace('X','23')):
            strlist = variantstr.split('_')
            strlist = list(filter(None, strlist)) # remove empty strings
            try:
                achr, astart, aend = parseRegionText(strlist[0]+":"+strlist[1]+"-"+str(int(strlist[1])+1), build)
                if str(achr).replace('23','X') == str(chrom) and astart >= startbp and astart <= endbp:
                    variantstr = variantstr.replace("_"+str(suffix),"") + "_"+str(suffix)
                    if len(variantstr.split('_')) == 5:
                        stdvariantlist.append(variantstr)
                    else:
                        raise InvalidUsage(f'Variant format not recognizable: {variant}', status_code=410)
                else:
                    stdvariantlist.append('.')
            except:
                raise InvalidUsage(f'Problem with variant {variant}')
        elif re.search("^\d+_\d+_*[A,T,G,C]*", variantstr.replace('X','23')):
            strlist = variantstr.split('_')
            strlist = list(filter(None, strlist)) # remove empty strings
            try:
                achr, astart, aend = parseRegionText(strlist[0]+":"+strlist[1]+"-"+str(int(strlist[1])+1), build)
                if str(achr).replace('23','X') == str(chrom) and astart >= startbp and astart <= endbp:
                    if len(strlist)==3:
                        aref=strlist[2]
                    else:
                        aref=''
                    stdvariantlist.append(fetchSNV(achr, astart, aref, build))
                else:
                    stdvariantlist.append('.')
            except:
                raise InvalidUsage(f'Problem with variant {variant}')
        else:
            raise InvalidUsage(f'Variant format not recognized: {variant}', status_code=410)
    return stdvariantlist


def decomposeVariant(variant_list):
    """
    Parameters
    ----------
    variantid_list : list
        list of str standardized variants in chr_pos_ref_alt_build format
        
    Returns
    -------
    A pandas.dataframe with chromosome, pos, reference and alternate alleles columns
    """
    chromlist = [x.split('_')[0] if len(x.split('_'))==5 else x for x in variant_list]
    chromlist = [int(x) for x in chromlist if x!="X"]
    poslist = [int(x.split('_')[1]) if len(x.split('_'))==5 else x for x in variant_list]
    reflist = [x.split('_')[2] if len(x.split('_'))==5 else x for x in variant_list]
    altlist = [x.split('_')[3] if len(x.split('_'))==5 else x for x in variant_list]
    df = pd.DataFrame({
        default_chromname: chromlist
        ,default_posname: poslist
        ,default_refname: reflist
        ,default_altname: altlist
        })
    return df


def addVariantID(gwas_data, chromcol, poscol, refcol, altcol, build):
    """
    
    Parameters
    ----------
    gwas_data : pandas.DataFrame
        Has a minimum of chromosome, position, reference and alternate allele columns.
    chromcol : str
        chromosome column name in gwas_data
    poscol : str
        position column name in gwas_data
    refcol : str
        reference allele column name in gwas_data
    altcol : str
        alternate allele column name in gwas_data

    Returns
    -------
    pandas.dataframe with list of standardized variant ID's in chrom_pos_ref_alt_build format added to gwas_data

    """
    varlist = []
    buildstr = 'b37'
    if build == 'hg38':
        buildstr = 'b38'
    chromlist = list(gwas_data[chromcol])
    poslist = list(gwas_data[poscol])
    reflist = list(gwas_data[refcol])
    altlist = list(gwas_data[altcol])
    for i in np.arange(gwas_data.shape[0]):
        chrom = chromlist[i]
        pos = poslist[i]
        ref = reflist[i]
        alt = altlist[i]
        varlist.append('_'.join([chrom,pos,ref,alt,buildstr]))
    gwas_data[default_snpname] = varlist
    return gwas_data



def subsetLocus(build, summaryStats, regiontext, chromcol, poscol, pcol):
    # regiontext format example: "1:205500000-206000000"
    if regiontext == "": regiontext = default_region
    print('Parsing region text')
    chrom, startbp, endbp = parseRegionText(regiontext, build)
    print(chrom,startbp,endbp)
    print('Subsetting GWAS data to entered region')            
    bool1 = [x == chrom for x in Xto23(list(summaryStats[chromcol]))]
    bool2 = [x>=startbp and x<=endbp for x in list(summaryStats[poscol])]
    bool3 = [not x for x in list(summaryStats.isnull().any(axis=1))]
    gwas_indices_kept = [ (x and y) and z for x,y,z in zip(bool1,bool2,bool3)]
    summaryStats = summaryStats.loc[ gwas_indices_kept ].copy()
    summaryStats.sort_values(by=[ poscol ], inplace=True)
    chromcolnum = list(summaryStats.columns).index(chromcol)
    summaryStats.reset_index(drop=True, inplace=True)
    summaryStats.iloc[:,chromcolnum] = Xto23(list(summaryStats[chromcol]))
    if summaryStats.shape[0] == 0:
        raise InvalidUsage('No data found for entered region', status_code=410)
    # Check for invalid p=0 rows:
    zero_p = [x for x in list(summaryStats[pcol]) if x==0]
    if len(zero_p)>0:
        raise InvalidUsage('P-values of zero detected; please replace with a non-zero p-value')
    return summaryStats, gwas_indices_kept


def getLeadSNPindex(leadsnpname, summaryStats, snpcol, pcol):
    lead_snp = leadsnpname
    snp_list = list(summaryStats.loc[:,snpcol])
    snp_list = [asnp.split(';')[0] for asnp in snp_list] # cleaning up the SNP names a bit
    if lead_snp=='': lead_snp = list(summaryStats.loc[ summaryStats.loc[:,pcol] == min(summaryStats.loc[:,pcol]) ].loc[:,snpcol])[0].split(';')[0]
    if lead_snp not in snp_list:
        raise InvalidUsage('Lead SNP not found', status_code=410)
    lead_snp_position_index = list(summaryStats.loc[:,snpcol]).index(lead_snp)
    return lead_snp_position_index

####################################
# LD Calculation from 1KG using PLINK (on-the-fly)
####################################

def resolve_plink_filepath(build, pop, chrom):
    """
    Returns the file path of the binary plink file
    """
    if chrom == 'X': chrom = 23
    try:
        chrom = int(chrom)
    except:
        raise InvalidUsage(f"Invalid chromosome {str(chrom)}", status_code=410)
    if chrom not in np.arange(1,24):
        raise InvalidUsage(f"Invalid chromosome {str(chrom)}", status_code=410)
    if pop not in valid_populations:
        raise InvalidUsage(f"{str(pop)} is not a recognized population", status_code=410)
    plink_filepath = ""
    if build.lower() in ["hg19","grch37"]:
        if chrom == 23:
            plink_filepath = os.path.join(MYDIR, "data", "1000Genomes_GRCh37", pop, "chrX")
        else:
            plink_filepath = os.path.join(MYDIR, "data", "1000Genomes_GRCh37", pop, f"chr{chrom}")
    elif build.lower() in ["hg38","grch38"]:
        if chrom == 23:
            plink_filepath = os.path.join(MYDIR, "data", "1000Genomes_GRCh38", "chrX")
        else:
            plink_filepath = os.path.join(MYDIR, "data", "1000Genomes_GRCh38", f"chr{chrom}")
    else:
        raise InvalidUsage(f'{str(build)} is not a recognized genome build')
    return plink_filepath

def plink_ldmat(build, pop, chrom, snp_positions, outfilename):
    plink_filepath = resolve_plink_filepath(build, pop, chrom)
    # make snps file to extract:
    snps = [f"chr{str(int(chrom))}:{str(int(position))}" for position in snp_positions]
    writeList(snps, outfilename + "_snps.txt")
    #plink_path = subprocess.run(args=["which","plink"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
    if build.lower() in ["hg19","grch37"]:
        if os.name == 'nt':
            plinkrun = subprocess.run(args=[
                "./plink.exe", '--bfile', plink_filepath
                , "--chr", str(chrom)
                , "--extract", outfilename + "_snps.txt"
                , "--from-bp", str(min(snp_positions))
                , "--to-bp", str(max(snp_positions))
                , "--r2", "square"
                , "--make-bed"
                , "--threads", "1"
                , "--out", outfilename
                ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        else:    
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
    elif build.lower() in ["hg38","grch38"]:
        popfile = os.path.join(MYDIR, 'data', '1000Genomes_GRCh38', str(pop)+'.txt')
        if os.name == 'nt':
            plinkrun = subprocess.run(args=[
                "./plink.exe", '--bfile', plink_filepath
                , "--keep", popfile # this is the difference in running GRCh38
                , "--chr", str(chrom)
                , "--extract", outfilename + "_snps.txt"
                , "--from-bp", str(min(snp_positions))
                , "--to-bp", str(max(snp_positions))
                , "--r2", "square"
                , "--make-bed"
                , "--threads", "1"
                , "--out", outfilename
                ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        else:    
            plinkrun = subprocess.run(args=[
                "./plink", '--bfile', plink_filepath
                , "--keep", popfile # this is the difference in running GRCh38
                , "--chr", str(chrom)
                , "--extract", outfilename + "_snps.txt"
                , "--from-bp", str(min(snp_positions))
                , "--to-bp", str(max(snp_positions))
                , "--r2", "square"
                , "--make-bed"
                , "--threads", "1"
                , "--out", outfilename
                ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)            
    else:
        raise InvalidUsage(f'{str(build)} is not a recognized genome build')
    if plinkrun.returncode != 0:
        raise InvalidUsage(plinkrun.stdout.decode('utf-8'), status_code=410)
    ld_snps = list(pd.read_csv(outfilename + ".bim", sep="\t", header=None).iloc[:,1])
    ldmat = np.matrix(pd.read_csv(outfilename + ".ld", sep="\t", header=None))
    return ld_snps, ldmat

def plink_ld_pairwise(build, lead_snp_position, pop, chrom, snp_positions, snp_pvalues, outfilename):
    # positions must be in hg19 coordinates
    # returns NaN for SNPs not in 1KG LD file; preserves order of input snp_positions
    plink_filepath = resolve_plink_filepath(build, pop, chrom)
    # make snps file to extract:
    snps = [f"chr{str(int(chrom))}:{str(int(position))}" for position in snp_positions]
    writeList(snps, outfilename + "_snps.txt")
    
    # Ensure lead snp is also present in 1KG; if not, choose next best lead SNP
    lead_snp = f"chr{str(int(chrom))}:{str(int(lead_snp_position))}"
    the1kg_snps = list(pd.read_csv(plink_filepath + ".bim", sep="\t", header=None).iloc[:,1])
    new_lead_snp = lead_snp
    new_lead_snp_position = int(lead_snp_position)
    while (new_lead_snp not in the1kg_snps) and (len(snp_positions) != 1):
        print(new_lead_snp + ' not in 1KG ' + str(len(snp_positions)) + ' SNPs left ')
        lead_snp_index = snp_positions.index(new_lead_snp_position)
        snp_positions.remove(new_lead_snp_position)
        del snp_pvalues[lead_snp_index]
        new_lead_snp_position = snp_positions[ snp_pvalues.index(min(snp_pvalues)) ]
        new_lead_snp = f"chr{str(int(chrom))}:{str(int(new_lead_snp_position))}"
    if len(snp_positions) == 0:
        raise InvalidUsage('No alternative lead SNP found in the 1000 Genomes', status_code=410)
    lead_snp = new_lead_snp
    lead_snp_position = new_lead_snp_position
    print('Lead SNP in use: ' + lead_snp)

    #plink_path = subprocess.run(args=["which","plink"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
    if build.lower() in ["hg19","grch37"]:
        if os.name == 'nt':
            plinkrun = subprocess.run(args=[
                "./plink.exe", '--bfile', plink_filepath
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
        else:
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
    elif build.lower() in ["hg38","grch38"]:
        popfile = os.path.join(MYDIR, 'data', '1000Genomes_GRCh38', str(pop)+'.txt')
        if os.name == 'nt':
            plinkrun = subprocess.run(args=[
                "./plink.exe", '--bfile', plink_filepath
                , "--keep", popfile # this is the difference in running GRCh38
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
        else:
            plinkrun = subprocess.run(args=[
                "./plink", '--bfile', plink_filepath
                , "--keep", popfile # this is the difference in running GRCh38
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
    else:
        raise InvalidUsage(f'{str(build)} is not a recognized genome build')
    if plinkrun.returncode != 0:
        raise InvalidUsage(plinkrun.stdout.decode('utf-8'), status_code=410)
    ld_results = pd.read_csv(outfilename + ".ld", delim_whitespace=True)
    available_r2_positions = ld_results[['BP_B', 'R2']]
    pos_df = pd.DataFrame({'pos': snp_positions})
    merged_df = pd.merge(pos_df, available_r2_positions, how='left', left_on="pos", right_on="BP_B", sort=False)[['pos', 'R2']]
    merged_df.fillna(-1, inplace=True)
    return merged_df, new_lead_snp_position


####################################
# Getting GTEx Data from Local MongoDB Database
####################################

# This is the main function to extract the data for a tissue and gene_id:
def get_gtex(version, tissue, gene_id):
    if version.upper() == "V8":
        db = client.GTEx_V8
        collapsed_genes_df = collapsed_genes_df_hg38
    elif version.upper() == "V7":
        db = client.GTEx_V7
        collapsed_genes_df = collapsed_genes_df_hg19
    
    tissue = tissue.replace(' ','_')
    #gene_id = gene_id.upper()
    ensg_name = ""
    if tissue not in db.list_collection_names():
        raise InvalidUsage(f'Tissue {tissue} not found', status_code=410)
    collection = db[tissue]
    if gene_id.startswith('ENSG'):
        i = list(collapsed_genes_df['ENSG_name']).index(gene_id)
        ensg_name = list(collapsed_genes_df['ENSG_name'])[i]
    elif gene_id in list(collapsed_genes_df['name']):
        i = list(collapsed_genes_df['name']).index(gene_id)
        ensg_name = list(collapsed_genes_df['ENSG_name'])[i]
    else:
        raise InvalidUsage(f'Gene name {gene_id} not found', status_code=410)
    #t1=datetime.now()
    results = list(collection.find({'gene_id': ensg_name}))
    #t2=datetime.now()
    #print(f'Time to retrieve gene collection: {t2-t1}')
    #t1=datetime.now()
    response = []
    try:
        response = results[0]['eqtl_variants']
    except:
        return pd.DataFrame([{'error': f'No eQTL data for {gene_id} in {tissue}'}])
    results_df = pd.DataFrame(response)
    chrom = int(list(results_df['variant_id'])[0].split('_')[0].replace('X','23'))
    positions = [ int(x.split('_')[1]) for x in list(results_df['variant_id']) ]
    #t2=datetime.now()
    #print(f'Reformatting time: {t2-t1}')
    #t1=datetime.now()
    variants_query = db.variant_table.find(
        { '$and': [ 
            { 'chr': chrom }, 
            { 'variant_pos': { '$gte': min(positions), '$lte': max(positions) } } 
        ]} 
    )
    #t2=datetime.now()
    #print(f'Time for variant table query: {t2-t1}')
    #t1=datetime.now()
    t1=datetime.now()
    variants_list = list(variants_query)
    t2=datetime.now()
    print(f'Time for variant query conversion to list: {t2-t1}')
    t1=datetime.now()
    variants_df = pd.DataFrame(variants_list)
    t2=datetime.now()
    print(f'Time to convert the variant list to dataframe: {t2-t1}')
    #print(variants_df.shape)
    #t2=datetime.now()
    #print(f'Reformatting time: {t2-t1}')
    #t1=datetime.now()
    variants_df = variants_df.drop(['_id'], axis=1)
    #t2=datetime.now()
    #print(f'_id dropping time: {t2-t1}')
    #t1=datetime.now()
    x = pd.merge(results_df, variants_df, on='variant_id')
    #t2=datetime.now()
    #print(f'Merging time: {t2-t1}')
    if version.upper() == "V7":
        x.rename(columns={'rs_id_dbSNP147_GRCh37p13': 'rs_id'}, inplace=True)
    elif version.upper() == "V8":
        x.rename(columns={'rs_id_dbSNP151_GRCh38p7': 'rs_id'}, inplace=True)
    return x

# Function to merge the GTEx data with a particular snp_list
def get_gtex_data(version, tissue, gene, snp_list, raiseErrors = False):
    build = "hg19"
    if version.upper() == "V8":
        build = "hg38"
    gtex_data = []
    rsids = True
    rsid_snps = [ x for x in snp_list if x.startswith('rs') ]
    b37_snps = [ x for x in snp_list if x.endswith('_b37') ]
    b38_snps = [ x for x in snp_list if x.endswith('_b38') ]
    if len(rsid_snps) > 0:
        rsids = True
    elif len(b37_snps) or len(b38_snps) > 0:
        rsids = False
    else:
        raise InvalidUsage('Variant naming format not supported; ensure all are rs ID\'s are formatted as chrom_pos_ref_alt_b37 eg. 1_205720483_G_A_b37')
    hugo_gene, ensg_gene = genenames(gene, build)
    print(f'Gathering eQTL data for {hugo_gene} ({ensg_gene}) in {tissue}')
    response_df = pd.DataFrame({})
    if version.upper() == "V7":
        response_df = get_gtex("V7", tissue, gene)
    elif version.upper() == "V8":
        response_df = get_gtex("V8", tissue, gene)
    if 'error' not in response_df.columns:
        eqtl = response_df
        if rsids:
            snp_df = pd.DataFrame(snp_list, columns=['rs_id'])
            #idx = pd.Index(list(snp_df['rs_id']))
            idx2 = pd.Index(list(eqtl['rs_id']))
            #snp_df = snp_df[~idx.duplicated()]
            eqtl = eqtl[~idx2.duplicated()]
            # print('snp_df.shape' + str(snp_df.shape))
            gtex_data = snp_df.reset_index().merge(eqtl, on='rs_id', how='left', sort=False).sort_values('index')
            # print('gtex_data.shape' + str(gtex_data.shape))
            # print(gtex_data)
        else:
            snp_df = pd.DataFrame(snp_list, columns=['variant_id'])
            gtex_data = snp_df.reset_index().merge(eqtl, on='variant_id', how='left', sort=False).sort_values('index')
    else:
        try:
            error_message = list(response_df['error'])[0]
            gtex_data = pd.DataFrame({})
        except:
            if raiseErrors:
                raise InvalidUsage("No response for tissue " + tissue.replace("_"," ") + " and gene " + hugo_gene + " ( " + ensg_gene + " )", status_code=410)
    return gtex_data


# This function simply merges the eqtl_data extracted with the snp_list, 
# then returns a list of the eQTL pvalues for snp_list (if available)
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

@app.errorhandler(413)
def request_entity_too_large(error):
    return 'File Too Large', 413

@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    response = jsonify(error.to_dict())
    response.status_code = error.status_code
    return response

@app.route("/populations")
def get1KGPopulations():
    populations = pd.read_csv(os.path.join(MYDIR, 'data/populations.tsv'), sep='\t')
    return jsonify(populations.to_dict(orient='list'))

@app.route("/genenames/<build>")
def getGeneNames(build):
    if build.lower() == "hg38":
        collapsed_genes_df = collapsed_genes_df_hg38
    elif build.lower() == "hg19":
        collapsed_genes_df = collapsed_genes_df_hg19
    return jsonify(list(collapsed_genes_df['name']))

@app.route("/genenames/<build>/<chrom>/<startbp>/<endbp>")
def getGenesInRange(build, chrom, startbp, endbp):
    collapsed_genes_df = collapsed_genes_df_hg19
    if build.lower() == "hg38":
        collapsed_genes_df = collapsed_genes_df_hg38
    regiontext = str(chrom) + ":" + startbp + "-" + endbp
    chrom, startbp, endbp = parseRegionText(regiontext, build)
    genes_to_draw = collapsed_genes_df.loc[ (collapsed_genes_df['chrom'] == ('chr' + str(chrom).replace('23','X'))) &
                                                    ( ((collapsed_genes_df['txStart'] >= startbp) & (collapsed_genes_df['txStart'] <= endbp)) | 
                                                      ((collapsed_genes_df['txEnd'] >= startbp  ) & (collapsed_genes_df['txEnd'] <= endbp  )) | 
                                                      ((collapsed_genes_df['txStart'] <= startbp) & (collapsed_genes_df['txEnd'] >= endbp  )) )]
    return jsonify(list(genes_to_draw['name']))

@app.route("/gtex/<version>/tissues_list")
def list_tissues(version):
    if version.upper() == "V8":
        db = client.GTEx_V8
    elif version.upper() == "V7":
        db = client.GTEx_V7
    tissues = list(db.list_collection_names())
    tissues.remove('variant_table')
    return jsonify(tissues)

@app.route("/gtex/<version>/<tissue>/<gene_id>")
def get_gtex_route(version, tissue, gene_id):
    x = get_gtex(version, tissue, gene_id)
    x = x.fillna(-1)
    return jsonify(x.to_dict(orient='records'))

@app.route("/gtex/<version>/<tissue>/<gene_id>/<variant>")
def get_gtex_variant(version, tissue, gene_id, variant):
    x = get_gtex(version, tissue, gene_id)
    response_df = x
    result = []
    if variant.startswith("rs"):
        result = response_df.loc[ response_df['rs_id'] == variant ]
    elif variant.endswith("_b37") or variant.endswith("_b38"):
        result = response_df.loc[ response_df['variant_id'] == variant ]
    else:
        raise InvalidUsage(f'variant name {variant} not found', status_code=410)
    if result.shape[0] == 0:
        raise InvalidUsage(f'variant name {variant} not found', status_code=410)
    return jsonify(result.to_dict(orient='records'))

@app.route("/previous_session", methods=['GET', 'POST'])
def prev_session():
    if request.method == 'POST':
        old_session_id = request.form['session-id']
        if old_session_id != '':
            my_session_id = old_session_id
            sessionfile =  f'session_data/form_data-{my_session_id}.json'
            genes_sessionfile = f'session_data/genes_data-{my_session_id}.json'
            SSPvalues_file = f'session_data/SSPvalues-{my_session_id}.json'
            sessionfilepath = os.path.join(MYDIR, 'static', sessionfile)
            genes_sessionfilepath = os.path.join(MYDIR, 'static', genes_sessionfile)
            SSPvalues_filepath = os.path.join(MYDIR, 'static', SSPvalues_file)
        else: # blank input
            raise InvalidUsage('Invalid input')
        # print(f'Session filepath: {sessionfilepath} is {str(os.path.isfile(sessionfilepath))}')
        # print(f'Genes filepath: {genes_sessionfilepath} is {str(os.path.isfile(genes_sessionfilepath))}')
        # print(f'SSPvalues filepath: {SSPvalues_filepath} is {str(os.path.isfile(SSPvalues_filepath))}')
        if not (os.path.isfile(sessionfilepath) and os.path.isfile(genes_sessionfilepath) and os.path.isfile(SSPvalues_filepath)):
            raise InvalidUsage(f'Could not locate session {my_session_id}')
        return render_template("plot.html", sessionfile = sessionfile, genesfile = genes_sessionfile, SSPvalues_file = SSPvalues_file, coloc2_file = coloc2_file, sessionid = my_session_id)
    return render_template('session_form.html')

@app.route("/session_id/<old_session_id>")
def prev_session_input(old_session_id):
    if old_session_id != '':
        my_session_id = old_session_id
        sessionfile =  f'session_data/form_data-{my_session_id}.json'
        genes_sessionfile = f'session_data/genes_data-{my_session_id}.json'
        SSPvalues_file = f'session_data/SSPvalues-{my_session_id}.json'
        sessionfilepath = os.path.join(MYDIR, 'static', sessionfile)
        genes_sessionfilepath = os.path.join(MYDIR, 'static', genes_sessionfile)
        SSPvalues_filepath = os.path.join(MYDIR, 'static', SSPvalues_file)
    else: # blank input
        raise InvalidUsage('Invalid input')
    print(f'Session filepath: {sessionfilepath} is {str(os.path.isfile(sessionfilepath))}')
    print(f'Genes filepath: {genes_sessionfilepath} is {str(os.path.isfile(genes_sessionfilepath))}')
    print(f'SSPvalues filepath: {SSPvalues_filepath} is {str(os.path.isfile(SSPvalues_filepath))}')
    if not (os.path.isfile(sessionfilepath) and os.path.isfile(genes_sessionfilepath) and os.path.isfile(SSPvalues_filepath)):
        raise InvalidUsage(f'Could not locate session {my_session_id}')
    return render_template("plot.html", sessionfile = sessionfile, genesfile = genes_sessionfile, SSPvalues_file = SSPvalues_file, coloc2_file = coloc2_file, sessionid = my_session_id)
    

@app.route("/update/<session_id>/<newgene>")
def update_colocalizing_gene(session_id, newgene):
    sessionfile = f'session_data/form_data-{session_id}.json'
    sessionfilepath = os.path.join(APP_STATIC, sessionfile)
    data = json.load(open(sessionfilepath, 'r'))
    gtex_tissues = data['gtex_tissues']
    snp_list = data['snps']
    gtex_version = data['gtex_version']
    if gtex_version.upper() not in available_gtex_versions:
        gtex_version = "V7"
    # gtex_data = {}
    for tissue in tqdm(gtex_tissues):
        data[tissue] = pd.DataFrame({})
        eqtl_df = get_gtex_data(gtex_version, tissue, newgene, snp_list)
        #eqtl_filepath = os.path.join(APP_STATIC, f'session_data/eqtl_df-{tissue}-{newgene}-{session_id}.txt')
        # if os.path.isfile(eqtl_filepath):
        if len(eqtl_df) > 0:
            eqtl_df.fillna(-1, inplace=True)
        data[tissue] = eqtl_df.to_dict(orient='records')
    # data.update(gtex_data)
    # json.dump(data, open(sessionfilepath, 'w'))

    return jsonify(data)


@app.route("/regionCheck/<build>/<regiontext>")
def regionCheck(build, regiontext):
    message = dict({'response': "OK"})
    if build not in ['hg19', 'hg38']:
        message['response'] = f'Unrecognized build: {build}'
        return jsonify(message)
    regiontext = regiontext.strip().replace(' ','').replace(',','').replace('chr','')
    if not re.search("^\d+:\d+-\d+$", regiontext.replace('X','23').replace('x','23')):
        message['response'] = 'Invalid coordinate format. e.g. 1:205,000,000-206,000,000'
        return jsonify(message)
    chrom = regiontext.split(':')[0].lower().replace('chr','').upper()
    pos = regiontext.split(':')[1]
    startbp = pos.split('-')[0].replace(',','')
    endbp = pos.split('-')[1].replace(',','')
    chromLengths = pd.read_csv(os.path.join(MYDIR, 'data', build + '_chrom_lengths.txt'), sep="\t", encoding='utf-8')
    chromLengths.set_index('sequence',inplace=True)
    if chrom == 'X' or chrom == '23':
        chrom = 23
        maxChromLength = chromLengths.loc['chrX', 'length']
        try:
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            message['response'] = "Invalid coordinate input"
            return jsonify(message)
    else:
        try:
            chrom = int(chrom)
            if chrom == 23:
                maxChromLength = chromLengths.loc['chrX', 'length']
            else:
                maxChromLength = chromLengths.loc['chr'+str(chrom), 'length']
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            message['response'] = "Invalid coordinate input"
            return jsonify(message)
    if chrom < 1 or chrom > 23:
        message['response'] = 'Chromosome input must be between 1 and 23'
    elif startbp > endbp:
        message['response'] = 'Starting chromosome basepair position is greater than ending basepair position'
    elif startbp > maxChromLength or endbp > maxChromLength:
        message['response'] = 'Start or end coordinates are out of range'
    elif (endbp - startbp) > genomicWindowLimit:
        message['response'] = f'Entered region size is larger than {genomicWindowLimit/10**6} Mbp'
        return jsonify(message)
    else:
        return jsonify(message)
    return jsonify(message)



@app.route('/', methods=['GET', 'POST'])
def index():
    data = {"success": False}
    
    # Initializing timing variables:
    t1_total = np.nan
    file_size = np.nan
    ldmat_file_size = np.nan
    upload_time = np.nan
    ldmat_upload_time = np.nan
    gwas_load_time = np.nan
    ld_pairwise_time = np.nan
    user_ld_load_time = np.nan
    gtex_one_gene_time = np.nan
    gene_list_time = np.nan
    SS_region_subsetting_time = np.nan
    gtex_all_queries_time = np.nan
    ldmat_time = np.nan
    ldmat_subsetting_time = np.nan
    SS_time = np.nan

    # coloc2-specific secondary dataset columns:
    BETA = 'BETA'
    SE = 'SE'
    ALT = 'A1'
    REF = 'A2'
    MAF = 'MAF'
    ProbeID = 'ProbeID'
    N = 'N'
    
    #######################################################
    # Uploading files
    #######################################################
    if request.method == 'POST':
        t1_total = datetime.now()
        if request.files.get('files[]'):
            t1 = datetime.now () # timer to get total upload time
            if 'files[]' in request.files:
                filenames = request.files.getlist('files[]')
                for file in filenames:
                    filename = secure_filename(file.filename)
                    filepath = os.path.join(MYDIR, app.config['UPLOAD_FOLDER'], filename)
                    file.save(filepath)
                gwas_filepath, ldmat_filepath, html_filepath = classify_files(filenames)
            upload_time = datetime.now() - t1

            #######################################################
            # Checking form input parameters
            #######################################################
            my_session_id = uuid.uuid4()
            coordinate = request.form['coordinate']
            gtex_version = "V7"
            collapsed_genes_df = collapsed_genes_df_hg19
            if coordinate.lower() == "hg38":
                gtex_version = "V8"
                collapsed_genes_df = collapsed_genes_df_hg38

            print(f'Session ID: {my_session_id}')
            print(f'Coordinate: {coordinate}')
            print(f'GTEx version: {gtex_version}')
            print('Loading file')
            t1 = datetime.now() # timing started for GWAS loading/subsetting/cleaning
            gwas_data = pd.read_csv(gwas_filepath, sep="\t", encoding='utf-8')
            
            inferVariant = request.form.get('markerCheckbox')
            chromcol, poscol, refcol, altcol = ('','','','')
            snpcol = ''
            columnnames = []
            if inferVariant:
                print('User would like variant locations inferred')
                snpcol = verifycol(formname = request.form['snp-col'], defaultname = default_snpname, filecolnames = gwas_data.columns, error_message_='Variant ID column not found')
                columnnames = [ snpcol ]
            else:
                chromcol = verifycol(formname = request.form['chrom-col'], defaultname = default_chromname, filecolnames = gwas_data.columns, error_message_=f"Chromosome column ({request.form['chrom-col']}) not found")
                poscol = verifycol(formname = request.form['pos-col'], defaultname = default_posname, filecolnames = gwas_data.columns, error_message_=f"Basepair position column ({request.form['pos-col']}) not found")
                refcol = verifycol(formname = request.form['ref-col'], defaultname = default_refname, filecolnames = gwas_data.columns, error_message_=f"Reference allele column ({request.form['ref-col']}) not found")
                altcol = verifycol(formname = request.form['alt-col'], defaultname = default_altname, filecolnames = gwas_data.columns, error_message_=f"Alternate allele column ({request.form['alt-col']}) not found")
                snpcol = request.form['snp-col'] # optional input in this case
                if snpcol != '':
                    snpcol = verifycol(formname = snpcol, defaultname = default_snpname, filecolnames = gwas_data.columns, error_message_='Variant ID column not found')
                    columnnames = [ chromcol, poscol, snpcol, refcol, altcol ]
                else:
                    columnnames = [ chromcol, poscol, refcol, altcol ]
                    print('No SNP ID column provided')
                # Check whether data types are ok:
                if not all(isinstance(x, int) for x in Xto23(list(gwas_data[chromcol]))):
                    raise InvalidUsage(f'Chromosome column ({chromcol}) contains unrecognizable values')
                if not all(isinstance(x, int) for x in list(gwas_data[poscol])):
                    raise InvalidUsage(f'Position column ({poscol}) has non-integer entries')
            pcol = verifycol(formname = request.form['pval-col'], defaultname = default_pname, filecolnames = gwas_data.columns, error_message_='P-value column not found')
            columnnames.append(pcol)
            if not all(isinstance(x, float) for x in list(gwas_data[pcol])):
                raise InvalidUsage(f'P-value column ({pcol}) has non-numeric entries')
            runcoloc2 = request.form.get('coloc2check')
            if runcoloc2:
                print('User would like COLOC2 results')
                betacol = verifycol(formname = request.form['beta-col'], defaultname = default_betaname, filecolnames = gwas_data.columns, error_message_='Beta column not found')
                stderrcol = verifycol(formname = request.form['stderr-col'], defaultname = default_betaname, filecolnames = gwas_data.columns, error_message_='Stderr column not found')
                numsamplescol = verifycol(formname = request.form['numsamples-col'], defaultname = default_nname, filecolnames = gwas_data.columns, error_message_='Number of samples column not found')
                mafcol = verifycol(formname = request.form['maf-col'], defaultname = default_mafname, filecolnames = gwas_data.columns, error_message_='MAF column not found')
                columnnames.extend([ betacol, stderrcol, numsamplescol, mafcol ])
                studytype = request.form['studytype']
                studytypedf = pd.DataFrame({'type': np.repeat(studytype,gwas_data.shape[0]).tolist()})
                if 'type' not in gwas_data.columns:
                    gwas_data = pd.concat([gwas_data, studytypedf], axis=1)
                columnnames.append('type')
                if studytype == 'cc': 
                    coloc2gwascolnames.append('Ncases')
                    numcases = request.form['numcases']
                    if not str(numcases).isdigit(): raise InvalidUsage('Number of cases entered must be an integer', status_code=410)
                    numcasesdf = pd.DataFrame({'Ncases': np.repeat(int(numcases), gwas_data.shape[0]).tolist()})
                    if 'Ncases' not in gwas_data.columns:
                        gwas_data = pd.concat([gwas_data, numcasesdf], axis=1)
                    columnnames.append('Ncases')
                if not all(isinstance(x, float) for x in list(gwas_data[betacol])):
                    raise InvalidUsage(f'Beta column ({betacol}) has non-numeric entries')
                if not all(isinstance(x, float) for x in list(gwas_data[stderrcol])):
                    raise InvalidUsage(f'Standard error column ({stderrcol}) has non-numeric entries')
                if not all(isinstance(x, int) for x in list(gwas_data[numsamplescol])):
                    raise InvalidUsage(f'Number of samples column ({numsamplescol}) has non-integer entries')
                if not all(isinstance(x, float) for x in list(gwas_data[mafcol])):
                    raise InvalidUsage(f'MAF column ({mafcol}) has non-numeric entries')
                


            # Further check column names provided:
            if len(set(columnnames)) != len(columnnames):
                raise InvalidUsage(f'Duplicate column names provided: {columnnames}')
            gwas_data = gwas_data[ columnnames ]
            if snpcol == '':
                gwas_data = addVariantID(gwas_data, chromcol, poscol, refcol, altcol)
                snpcol = default_snpname

            # LD:
            pops = request.form['LD-populations']
            if len(pops) == 0: pops = 'EUR'
            print('Populations:', pops)
            
            # GTEx tissues and genes:
            gtex_tissues = request.form.getlist('GTEx-tissues')
            print('GTEx tissues:',gtex_tissues)
            gtex_genes = request.form.getlist('region-genes')
            if len(gtex_tissues) > 0 and len(gtex_genes) == 0:
                raise InvalidUsage('Please select one or more genes to complement your GTEx tissue(s) selection')
            elif len(gtex_genes) > 0 and len(gtex_tissues) == 0:
                raise InvalidUsage('Please select one or more tissues to complement your GTEx gene(s) selection')
            # old code remnant:
            # if gtex_version=="V7": gene='ENSG00000174502.14'
            # if gtex_version=="V8": gene='ENSG00000174502.18'
            gene = gtex_genes[0]



            # Set-based P override:            
            setbasedP = request.form['setbasedP']
            if setbasedP=='':
                setbasedP = 'default'
            else:
                try:
                    setbasedP = float(setbasedP)
                    if setbasedP < 0 or setbasedP > 1:
                        raise InvalidUsage('Set-based p-value threshold given is not between 0 and 1')
                except:
                    raise InvalidUsage('Invalid value provided for the set-based p-value threshold. Value must be numeric between 0 and 1.')

            # Ensure custom LD matrix and GWAS files are sorted for accurate matching:
            if ldmat_filepath != '' and poscol != '' and not isSorted(list(gwas_data[poscol])):
                raise InvalidUsage('GWAS data input is not sorted and may not match with the LD matrix', status_code=410)

            regionstr = request.form['locus']
            if regionstr == "": regionstr = default_region
            leadsnpname = request.form['leadsnp']

            #######################################################
            # Standardizing variant ID's to chrom_pos_ref_alt_build format
            #######################################################
            if inferVariant:
                # standardize variant id's:
                variant_list = mapAndCleanSNPs(list(gwas_data[snpcol]), regionstr, coordinate)
                if all(x=='.' for x in variant_list):
                    raise InvalidUsage(f'None of the variants provided could be mapped to {regiontext}!', status_code=410)
                # get the chrom, pos, ref, alt info from the standardized variant_list
                vardf = decomposeVariant(variant_list)
                gwas_data = pd.concat([vardf, gwas_data], axis=1)
                chromcol = default_chromname
                poscol = default_posname
                refcol = default_refname
                altcol = default_altname
            
            #######################################################
            # Subsetting GWAS file
            #######################################################
            gwas_data, gwas_indices_kept = subsetLocus(coordinate, gwas_data, regionstr, chromcol, poscol, pcol)
            lead_snp_position_index = getLeadSNPindex(leadsnpname, gwas_data, snpcol, pcol)
            lead_snp_position = gwas_data.iloc[lead_snp_position_index,:][poscol]
            positions = list(gwas_data[poscol])
            snp_list = list(gwas_data[snpcol])
            snp_list = [asnp.split(';')[0] for asnp in snp_list] # cleaning up the SNP names a bit
            lead_snp = snp_list[ lead_snp_position_index ]
            pvals = list(gwas_data[pcol])
            chrom, startbp, endbp = parseRegionText(regionstr, coordinate)

            gwas_load_time = datetime.now() - t1

            ####################################################################################################
            # Get LD:
            if ldmat_filepath == '':
                t1 = datetime.now() # timing started for pairwise LD
                print('Calculating pairwise LD using PLINK')
                #ld_df = queryLD(lead_snp, snp_list, pops, ld_type)
                ld_df, new_lead_snp_position = plink_ld_pairwise(coordinate, lead_snp_position, pops, chrom, positions, pvals, os.path.join(MYDIR, "static", "session_data", f"ld-{my_session_id}"))
                if new_lead_snp_position != lead_snp_position:
                    lead_snp_position_index = list(gwas_data[poscol]).index(new_lead_snp_position)
                    lead_snp = snp_list[ lead_snp_position_index ]
                    lead_snp_position = new_lead_snp_position
                r2 = list(ld_df['R2'])
                ld_pairwise_time = datetime.now() - t1
            else:
                print('---------------------------------')
                print('Loading user-supplied LD matrix')
                print('---------------------------------')
                t1 = datetime.now() # timer started for loading user-defined LD matrix
                ld_mat = pd.read_csv(ldmat_filepath, sep="\t", encoding='utf-8', header=None)
                gwas_data_numRows = gwas_data.shape[0]
                ld_mat = ld_mat.loc[ gwas_indices_kept, gwas_indices_kept ]
                r2 = list(ld_mat.iloc[:, lead_snp_position_index])
                ld_mat = np.matrix(ld_mat)
                if not ((ld_mat.shape[0] == ld_mat.shape[1]) and (ld_mat.shape[0] == gwas_data.shape[0])):
                    raise InvalidUsage('GWAS and LD matrix input have different dimensions', status_code=410)
                user_ld_load_time = datetime.now() - t1
            
            data = {}
            data['snps'] = snp_list
            data['pvalues'] = list(gwas_data[pcol])
            data['lead_snp'] = lead_snp
            data['ld_values'] = r2
            data['positions'] = positions
            data['chrom'] = chrom
            data['startbp'] = startbp
            data['endbp'] = endbp
            data['ld_populations'] = pops
            data['gtex_tissues'] = gtex_tissues
            data['gene'] = genenames(gene, coordinate)[0]
            data['coordinate'] = coordinate
            data['gtex_version'] = gtex_version
            data['set_based_p'] = setbasedP
            SSlocustext = request.form['SSlocus']
            # SSlocus defined below
            
            
            #######################################################
            # Loading any secondary datasets uploaded
            #######################################################
            secondary_datasets = {}
            table_titles = []
            if html_filepath != '':
                print('Loading secondary datasets provided')
                with open(html_filepath, encoding='utf-8', errors='replace') as f:
                    html = f.read()
                    if (not html.startswith('<h3>')) and (not html.startswith('<html>')) and (not html.startswith('<table>') and (not html.startswith('<!DOCTYPE html>'))):
                        raise InvalidUsage('Secondary dataset(s) provided are not formatted correctly. Please use the merge_and_convert_to_html.py script for formatting.', status_code=410)
                soup = bs(html, 'lxml')
                table_titles = soup.find_all('h3')
                table_titles = [x.text for x in table_titles]
                tables = pd.read_html(html_filepath)
                for i in np.arange(len(tables)):
                    table = tables[i]
                    secondary_datasets[table_titles[i]] = table.fillna(-1).to_dict(orient='records')
            data['secondary_dataset_titles'] = table_titles
            if runcoloc2:
                data['secondary_dataset_colnames'] = ['CHR', 'POS', 'SNPID', 'PVAL', BETA, SE, 'N', ALT, REF, MAF, ProbeID]
            else:
                data['secondary_dataset_colnames'] = [CHROM, BP, SNP, P]
            data.update(secondary_datasets)

            ####################################################################################################
            t1 = datetime.now() # set timer for extracting GTEx data for selected gene:
            # Get GTEx data for the tissues and SNPs selected:
            gtex_data = {}
            if len(gtex_tissues)>0:
                print('Gathering GTEx data')
                for tissue in tqdm(gtex_tissues):
                    eqtl_df = get_gtex_data(gtex_version, tissue, gene, snp_list, raiseErrors=True) # for the full region (not just the SS region)
                    if len(eqtl_df) > 0:
                        eqtl_df.fillna(-1, inplace=True)
                    gtex_data[tissue] = eqtl_df.to_dict(orient='records')
            data.update(gtex_data)
            gtex_one_gene_time = datetime.now() - t1
            
            ####################################################################################################
            # Checking that there is at least one secondary dataset for colocalization
            ####################################################################################################
            if len(gtex_tissues)==0 and html_filepath == '':
                raise InvalidUsage('Please provide at least one secondary dataset or select at least one GTEx tissue for colocalization analysis')
            
            ####################################################################################################
            t1 = datetime.now() # timer for determining the gene list            
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
            gene_list_time = datetime.now() - t1

            ####################################################################################################
            # 1. Determine the region to calculate the Simple Sum (SS):
            if SSlocustext != '':
                SSchrom, SS_start, SS_end = parseRegionText(SSlocustext, coordinate)
            else:
                SS_start = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ][poscol])[0] - one_sided_SS_window_size
                SS_end = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ][poscol])[0] + one_sided_SS_window_size
            data['SS_region'] = [SS_start, SS_end]
            
            ####################################################################################################
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

            ####################################################################################################
            # # Getting Simple Sum P-values
            # 2. Subset the region (step 1 was determining the region to do the SS calculation on - see above SS_start and SS_end variables):
            t1 = datetime.now() # timer for subsetting SS region
            print('SS_start: ' + str(SS_start))
            print('SS_end:' + str(SS_end))
            chromList = [('chr' + str(chrom).replace('23','X')), str(chrom).replace('23','X')]
            gwas_chrom_col = pd.Series([str(x) for x in list(gwas_data[chromcol])])
            SS_chrom_bool = [x for x in gwas_chrom_col.isin(chromList) if x == True]
            SS_indices = SS_chrom_bool & (gwas_data[poscol] >= SS_start) & (gwas_data[poscol] <= SS_end)
            SS_gwas_data = gwas_data.loc[ SS_indices ]
            if runcoloc2:
                coloc2_gwasdf = SS_gwas_data.rename(columns={
                    chromcol: 'CHR'
                    ,poscol: 'POS'
                    ,snpcol: 'SNPID'
                    ,pcol: 'PVAL'
                    ,refcol: REF
                    ,altcol: ALT
                    ,betacol: BETA
                    ,stderrcol: SE
                    ,mafcol: MAF
                    ,numsamplescol: numsamplescol
                })
                coloc2_gwasdf = coloc2_gwasdf.reindex(columns=coloc2gwascolnames)
            #print(gwas_data.shape)
            #print(SS_gwas_data.shape)
            #print(SS_gwas_data)
            if SS_gwas_data.shape[0] == 0: InvalidUsage('No data points found for entered Simple Sum region', status_code=410)
            PvaluesMat = [list(SS_gwas_data[pcol])]
            SS_snp_list = list(SS_gwas_data[snpcol])
            # SS_snp_list = [asnp.split(';')[0] for asnp in SS_snp_list] # cleaning up the SNP names a bit
            if SSlocustext == '':
                SSlocustext = str(chrom) + ":" + str(SS_start) + "-" + str(SS_end)
            SS_snp_list = mapAndCleanSNPs(SS_snp_list, SSlocustext, coordinate)
            SS_positions = list(SS_gwas_data[poscol])
            if len(SS_positions) != len(set(SS_positions)):
                dups = set([x for x in SS_positions if SS_positions.count(x) > 1])
                raise InvalidUsage('Duplicate chromosome positions detected at: ' + str(dups))
            # Extra file written:
            gwas_df = pd.DataFrame({
                'Position': SS_positions,
                'SNP': SS_snp_list,
                'P': list(SS_gwas_data[pcol])
            })
            gwas_df.to_csv(os.path.join(MYDIR, 'static', f'session_data/gwas_df-{my_session_id}.txt'), index=False, encoding='utf-8', sep="\t")
            SS_region_subsetting_time = datetime.now() - t1
            
            ####################################################################################################
            # 3. Determine the genes to query
            #query_genes = list(genes_to_draw['name'])
            query_genes = gtex_genes
            coloc2eqtl_df = pd.DataFrame({})
            # 4. Query and extract the eQTL p-values for all tissues & genes from GTEx
            t1 = datetime.now() # timer set to check how long data extraction from Mongo takes
            if len(gtex_tissues)>0:
                print('Obtaining eQTL p-values for selected tissues and surrounding genes')
                for tissue in gtex_tissues:
                    for agene in query_genes:
                        gtex_eqtl_df = get_gtex_data(gtex_version, tissue, agene, SS_snp_list)
#                        print('len(gtex_eqtl_df) '+ str(len(gtex_eqtl_df)))
#                        print('gtex_eqtl_df.shape ' + str(gtex_eqtl_df.shape))
#                        print('len(SS_snp_list) ' + str(len(SS_snp_list)))
                        if len(gtex_eqtl_df) > 0:
                            #gtex_eqtl_df.fillna(-1, inplace=True)
                            pvalues = list(gtex_eqtl_df['pval'])
                            if runcoloc2:
                                tempdf = gtex_eqtl_df.rename(columns={
                                    'rs_id': 'SNPID'
                                    ,'pval': 'PVAL'
                                    ,'beta': BETA
                                    ,'se': SE
                                    ,'sample_maf': MAF
                                    ,'chr': 'CHR'
                                    ,'variant_pos': 'POS'
                                    ,'ref': REF
                                    ,'alt': ALT
                                })
                                numsamples = round(tempdf['ma_count'][0] / tempdf[MAF][0])
                                numsampleslist = np.repeat(numsamples, tempdf.shape[0]).tolist()
                                tempdf = pd.concat([tempdf,pd.Series(numsampleslist,name='N')],axis=1)
                                probeid = str(tissue) + ':' + str(agene)
                                probeidlist = np.repeat(probeid, tempdf.shape[0]).tolist()
                                tempdf = pd.concat([tempdf, pd.Series(probeidlist,name='ProbeID')],axis=1)
                                tempdf = tempdf.reindex(columns = coloc2eqtlcolnames)
                                coloc2eqtl_df = pd.concat([coloc2eqtl_df, tempdf], axis=0)
                        else:
                            pvalues = np.repeat(np.nan, len(SS_snp_list))
                        PvaluesMat.append(pvalues)
#                        print(f'tissue: {tissue}, gene: {gene}, len(pvalues): {len(pvalues)}')
#                        print(f'len(SS_positions): {len(SS_positions)}, len(SS_snp_list): {len(SS_snp_list)}')
                        
                        # Extra files written:
#                        eqtl_df = pd.DataFrame({
#                            'Position': SS_positions,
#                            'SNP': SS_snp_list,
#                            'P': pvalues
#                        })
#                        eqtl_df.to_csv(os.path.join(MYDIR, 'static', f'session_data/eqtl_df-{tissue}-{agene}-{my_session_id}.txt'), index=False, encoding='utf-8', sep="\t")
#                        print(f'Time to extract eQTLs for {tissue} and {agene}:' + str(datetime.now()-t1))

            gtex_all_queries_time = datetime.now() - t1

            ####################################################################################################
            # 4.2 Extract user's secondary datasets' p-values
            ####################################################################################################
            if len(secondary_datasets)>0:
                if runcoloc2:
                    print('Saving uploaded secondary datasets for coloc2 run')
                    for i in np.arange(len(secondary_datasets)):
                        secondary_dataset = tables[i]
                        secondary_dataset.set_index('SNPID', inplace=True)
                        secondary_dataset = secondary_dataset[~secondary_dataset.index.duplicated()]
                        # merge to keep only SNPs already present in the GWAS/primary dataset (SS subset):
                        snp_df = pd.DataFrame(SS_snp_list, columns=['SNPID'])
                        secondary_data = snp_df.reset_index().merge(secondary_dataset, left_on=SNP, right_on='SNPID', how='left', sort=False).sort_values('index')
                        pvalues = list(secondary_dataset['PVAL'])
                        PvaluesMat.append(pvalues)
                        coloc2eqtl_df = pd.concat([coloc2eqtl_df, secondary_data.reindex(columns = coloc2eqtlcolnames)], axis=0)
                else:
                    print('Obtaining p-values for uploaded secondary dataset(s)')
                    for i in np.arange(len(secondary_datasets)):
                        secondary_dataset = tables[i]
                        # remove duplicate SNPs
                        idx = pd.Index(list(secondary_dataset[SNP]))
                        secondary_dataset = secondary_dataset[~idx.duplicated()]
                        # merge to keep only SNPs already present in the GWAS/primary dataset (SS subset):
                        snp_df = pd.DataFrame(SS_snp_list, columns=[SNP])
                        secondary_data = snp_df.reset_index().merge(secondary_dataset, on=SNP, how='left', sort=False).sort_values('index')
                        pvalues = list(secondary_data[P])
                        PvaluesMat.append(pvalues)

            ####################################################################################################
            print('Extracting LD matrix')
            # 5. Get the LD matrix via PLINK subprocess call or use user-provided LD matrix:
            t1 = datetime.now() # timer for calculating the LD matrix
            plink_outfilename = f'session_data/ld-{my_session_id}'
            plink_outfilepath = os.path.join(MYDIR, 'static', plink_outfilename)
            if ldmat_filepath != '':
                ld_mat_snps = [f'chr{chrom}:{SS_pos}' for SS_pos in SS_positions]
                ld_mat_positions = [int(snp.split(":")[1]) for snp in ld_mat_snps]
                ld_mat = ld_mat[SS_indices ][:,SS_indices]
            else:
                ld_mat_snps, ld_mat = plink_ldmat(coordinate, pops, chrom, SS_positions, plink_outfilepath)
                ld_mat_positions = [int(snp.split(":")[1]) for snp in ld_mat_snps]
            np.fill_diagonal(ld_mat, np.diag(ld_mat) + ld_mat_diag_constant)
            ldmat_time = datetime.now() - t1
            
            ####################################################################################################
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
            ldmat_subsetting_time = datetime.now() - t1
            ####
            
            ####################################################################################################
            print('Calculating Simple Sum stats')
            t1 = datetime.now() # timer for Simple Sum calculation time
            Rscript_code_path = os.path.join(MYDIR, 'getSimpleSumStats.R')
            # Rscript_path = subprocess.run(args=["which","Rscript"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
            SSresult_path = os.path.join(MYDIR, 'static', f'session_data/SSPvalues-{my_session_id}.txt')
                        
            RscriptRun = subprocess.run(args=['Rscript', Rscript_code_path, Pvalues_filepath, ldmatrix_filepath, '--set_based_p', str(setbasedP), '--outfilename', SSresult_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            if RscriptRun.returncode != 0:
                raise InvalidUsage(RscriptRun.stdout, status_code=410)
            SSdf = pd.read_csv(SSresult_path, sep='\t', encoding='utf-8')

            # SSPvalues = RscriptRun.stdout.replace('\n',' ').split(' ')
            # SSPvalues = [float(SSP) for SSP in SSPvalues if SSP!='']
            
            #SSPvalues, num_SNP_used_for_SS, comp_used = getSimpleSumStats.get_simple_sum_p(np.asarray(PvaluesMat), np.asarray(ld_mat))
            
            SSPvalues = SSdf['Pss'].tolist()
            num_SNP_used_for_SS = SSdf['n'].tolist()
            comp_used = SSdf['comp_used'].tolist()
            
            for i in np.arange(len(SSPvalues)):
                if SSPvalues[i] > 0:
                    SSPvalues[i] = np.format_float_scientific((-np.log10(SSPvalues[i])), precision=2)
            SSPvaluesMatGTEx = np.empty(0)
            SSPvaluesSecondary = []
            if len(gtex_tissues)>0:
                SSPvaluesMatGTEx = np.array(SSPvalues[0:(len(gtex_tissues) * len(query_genes))]).reshape(len(gtex_tissues), len(query_genes))
            if len(SSPvalues) > len(gtex_tissues) * len(query_genes):
                SSPvaluesSecondary = SSPvalues[(len(gtex_tissues) * len(query_genes)) : (len(SSPvalues))]
            SSPvalues_dict = {
                'Genes': query_genes
                ,'Tissues': gtex_tissues
                ,'Secondary_dataset_titles': table_titles
                ,'SSPvalues': SSPvaluesMatGTEx.tolist() # GTEx pvalues
                ,'SSPvalues_secondary': SSPvaluesSecondary
                ,'Num_SNPs_Used_for_SS': [int(x) for x in num_SNP_used_for_SS]
                ,'Computation method': comp_used
            }
            print(SSPvalues_dict)
            SSPvalues_file = f'session_data/SSPvalues-{my_session_id}.json'
            SSPvalues_filepath = os.path.join(MYDIR, 'static', SSPvalues_file)
            json.dump(SSPvalues_dict, open(SSPvalues_filepath, 'w'))
            SS_time = datetime.now() - t1
            t2_total = datetime.now() - t1_total

            ####################################################################################################
            if runcoloc2:
                print('Calculating COLOC2 stats')
                coloc2gwasfilepath = os.path.join(MYDIR, 'static', f'session_data/coloc2gwas_df-{my_session_id}.txt')
                coloc2_gwasdf.dropna().to_csv(coloc2gwasfilepath, index=False, encoding='utf-8', sep="\t")
                coloc2eqtlfilepath = os.path.join(MYDIR, 'static', f'session_data/coloc2eqtl_df-{my_session_id}.txt')
                coloc2eqtl_df.dropna().to_csv(coloc2eqtlfilepath, index=False, encoding='utf-8', sep="\t")
                Rscript_code_path = os.path.join(MYDIR, 'coloc2', 'run_coloc2.R')
                coloc2result_path = os.path.join(MYDIR, 'static', f'session_data/coloc2result_df-{my_session_id}.txt')
                coloc2RscriptRun = subprocess.run(args=['Rscript', Rscript_code_path, coloc2gwasfilepath, coloc2eqtlfilepath, '--outfilename', coloc2result_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                if coloc2RscriptRun.returncode != 0:
                    raise InvalidUsage(coloc2RscriptRun.stdout, status_code=410)
                coloc2df = pd.read_csv(coloc2result_path, sep='\t', encoding='utf-8').fillna(-1)
                # save as json:
                coloc2_dict = {
                    'ProbeID': coloc2df['ProbeID'].tolist()
                    ,'PPH4abf': coloc2df['PPH4abf'].tolist()
                }
            else:
                coloc2_dict = {
                    'ProbeID': []
                    ,'PPH4abf': []
                }
            coloc2_file = f'session_data/coloc2result-{my_session_id}.json'
            coloc2_filepath = os.path.join(MYDIR, 'static', coloc2_file)
            json.dump(coloc2_dict, open(coloc2_filepath,'w'))

            ####################################################################################################
            

            timing_file = f'session_data/times-{my_session_id}.txt'
            timing_file_path = os.path.join(MYDIR, 'static', timing_file)
            with open(timing_file_path, 'w') as f:
                f.write('-----------------------------------------------------------\n')
                f.write(' Times Report\n')
                f.write('-----------------------------------------------------------\n')
                f.write(f'File size: {file_size/1000:.0f} KB\n')
                f.write(f'Upload time: {upload_time}\n')
                if not np.isnan(ldmat_upload_time):
                    f.write(f'LD matrix file size: {ldmat_file_size/1000} KB\n')
                    f.write(f'LD matrix upload time: {ldmat_upload_time}\n')
                    f.write(f'LD matrix loading and subsetting time: {user_ld_load_time}\n')
                f.write(f'GWAS load time: {gwas_load_time}\n')
                f.write(f'Pairwise LD calculation time: {ld_pairwise_time}\n')
                f.write(f'Extracting GTEx eQTLs for user-specified gene: {gtex_one_gene_time}\n')
                f.write(f'Finding all genes to draw and query time: {gene_list_time}\n')
                f.write(f'Number of genes found in the region: {genes_to_draw.shape[0]}\n')
                f.write(f'Time to subset Simple Sum region: {SS_region_subsetting_time}\n')
                f.write(f'Time to extract all eQTL data from {len(gtex_tissues)} tissues and {len(query_genes)} genes: {gtex_all_queries_time}\n')
                f.write(f'Time for calculating the LD matrix: {ldmat_time}\n')
                f.write(f'Time for subsetting the LD matrix: {ldmat_subsetting_time}\n')
                num_nmiss_tissues = -1 # because first row are the GWAS pvalues
                for i in np.arange(len(PvaluesMat.tolist())):
                    if not np.isnan(PvaluesMat.tolist()[i][0]):
                        num_nmiss_tissues += 1
                f.write(f'Time for calculating the Simple Sum P-values: {SS_time}\n')
                f.write(f'For {num_nmiss_tissues} pairwise calculations out of {PvaluesMat.shape[0]-1}\n')
                if num_nmiss_tissues != 0: f.write(f'Time per Mongo query: {gtex_all_queries_time/num_nmiss_tissues}\n')
                if num_nmiss_tissues != 0: f.write(f'Time per SS calculation: {SS_time/num_nmiss_tissues}\n')
                f.write(f'Total time: {t2_total}\n')
            
            # Compress session data files for easy download:
            print('Compressing data for downloading')
            downloadfile = f'session_data/LocusFocus_session_data-{my_session_id}.tar.gz'
            downloadfilepath = os.path.join(MYDIR, 'static', downloadfile)
            files_to_compress = f'session_data/*{my_session_id}*'
            files_to_compress_path = os.path.join(MYDIR, 'static', files_to_compress)
            compressrun = subprocess.run(args=['tar', 'zcvf', downloadfilepath, files_to_compress_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            # if compressrun.returncode != 0:
            #     raise InvalidUsage(compressrun.stdout.decode('utf-8'), status_code=410)

            return render_template("plot.html", sessionfile = sessionfile, genesfile = genes_sessionfile, SSPvalues_file = SSPvalues_file, coloc2_file = coloc2_file, sessionid = my_session_id)
        return render_template("invalid_input.html")
    return render_template("index.html")


@app.route('/downloaddata/<my_session_id>')
def downloaddata(my_session_id):
     downloadfile = f'session_data/LocusFocus_session_data-{my_session_id}.tar.gz'
     downloadfilepath = os.path.join(MYDIR, 'static', downloadfile)
     return send_file(downloadfilepath, as_attachment=True)


app.config['SITEMAP_URL_SCHEME'] = 'https'


@ext.register_generator
def index():
    # Not needed if you set SITEMAP_INCLUDE_RULES_WITHOUT_PARAMS=True
    # yield 'index', {}
    urls = ['locusfocus.research.sickkids.ca',
        'https://locusfocus.research.sickkids.ca/session_id/00dfdb4d-c86a-423b-adc7-4740b7b43695',
        'https://locusfocus.research.sickkids.ca/previous_session']
    return urls



if __name__ == "__main__":
    app.run()

