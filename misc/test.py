import os
import pandas as pd
import numpy as np
import re
import gzip
import dask.dataframe as dd
from dask.delayed import delayed
import subprocess
import pysam


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



def parseRegionText(regiontext, build):
    if build not in ['hg19', 'hg38']:
        raise InvalidUsage(f'Unrecognized build: {build}', status_code=410)
    regiontext = regiontext.strip().replace(' ','').replace(',','').replace('chr','')
    if not re.search("^\d+:\d+-\d+$", regiontext.replace('X','23').replace('x','23')):
       raise InvalidUsage('Invalid coordinate format. e.g. 1:205,000,000-206,000,000', status_code=410)
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
            raise InvalidUsage("Invalid coordinates input", status_code=410)
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
            raise InvalidUsage("Invalid coordinates input", status_code=410)
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


def subsetLocus(build, summaryStats, regiontext, columnnames):
    # regiontext format example: "1:205500000-206000000"
    if regiontext == "": regiontext = default_region
    chromcol, poscol, snpcol, pcol = columnnames
    print('Parsing region text')
    chrom, startbp, endbp = parseRegionText(regiontext, build)
    if (endbp - startbp) > genomicWindowLimit:
        raise InvalidUsage(f'Entered region size is larger than {genomicWindowLimit} bp', status_code=410)
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



def getNumHeaderLines(vcf_filename, num_lines_to_check = 1000):
    num_header_lines = 0
    num_lines_checked = 0
    with gzip.open(vcf_filename, 'rb') as f:
        nextline = f.readline().decode('utf-8')
        num_lines_checked += 1
        while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
            num_header_lines += 1
            nextline = f.readline().decode('utf-8')
            num_lines_checked += 1
    return num_header_lines


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

#### Testing functions above
## Some global variables:
default_chromname = "#CHROM"
default_posname = "POS"
default_snpname = "ID"
default_refname = "REF"
default_altname = "ALT"

MYDIR = os.getcwd()
coordinates = 'hg38'
gwasdata = pd.read_csv(os.path.join(os.getcwd(), 'data', 'sample_datasets', 'MI_GWAS_2019_1_205500-206000kbp_hg38.tsv'), sep='\t', encoding='utf-8')
genomicWindowLimit = 2e6
default_region = "1:205500000-206000000"

## Test mapAndCleanSNPs function
regiontext = '1:205000000-206000000'
variantlist = ['rs2211330', '1_205001063_T_A,C', '1:205001063_T_C', 'rs7512462', '1_205930467']
variantlist.extend(['X_115000731','rs62599779', 'X_115012777_CAA_C', 'rs782251282', 
                    'X_115012777_C_CA,CAA,CAAA','X_115012777_C_CAA', 'X_115013906_C_T_b38', 
                    'X_115013906_C_T_b38', 'rs17095917'])
chrom, startbp, endbp = parseRegionText(regiontext, coordinates)
mapAndCleanSNPs(variantlist, regiontext, coordinates)
regiontext = 'X:114000000-116000000'
chrom, startbp, endbp = parseRegionText(regiontext, coordinates)
mapAndCleanSNPs(variantlist, regiontext, coordinates)

regiontext = '1:205000000-206000000'
variantlist = list(gwasdata['SNP'])
stdvarlist = mapAndCleanSNPs(variantlist, regiontext, coordinates)
if all(x=='.' for x in stdvarlist):
    print(f'None of the variants provided could be mapped to {regiontext}!')
vardf = decomposeVariant(stdvarlist)

gwasdata.drop(['chr','variant_pos','ref','alt'], axis=1, inplace=True)
gwasdata2 = pd.concat([vardf, gwasdata], axis=1)

gwasdata3 = addVariantID(gwasdata2, default_chromname, default_posname, default_refname, default_altname)

## Test subsetLocus function
columnnames = ['chr', 'variant_pos', 'SNP', 'ref', 'alt', 'P']
gwas_data, gwas_indices_kept = subsetLocus(build = coordinates, summaryStats = gwasdata, regiontext = regiontext, columnnames = columnnames)

