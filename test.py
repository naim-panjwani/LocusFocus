import pandas as pd
data = pd.read_csv('data/populations.tsv',sep="\t", encoding='utf-8')
data['lead_snp'] = 'rs662702'
subdata = data[['Population Code','lead_snp']]
subdata.head()

import os
MYDIR = os.getcwd()

test1 = "1"
print(test1.replace('chr',''))
test2 = 'chr1'
print(test2.replace('chr',''))

chroms = pd.read_csv('data/hg19_chrom_lengths.txt', sep="\t", encoding='utf-8')
# chroms.head()
# chroms[1]
# chroms[1] = [int(chr.strip().replace(',','')) for chr in chroms[1]]
# chroms.head()
# chroms[0] = [chr.strip() for chr in chroms[0]]

# chroms.to_csv('data/hg19_chrom_lengths.txt', index=False, encoding='utf-8', sep="\t")

gwas_data = pd.read_csv('data/MI_GWAS_2019_1_205500-206000kbp.tsv', sep="\t")
chromcol='#CHROM'
poscol='BP'
snpcol='SNP'
pcol='P'
regiontext = "1:205500000-206000000"
lead_snp = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ]['SNP'])[0]
collapsed_genes_df = pd.read_csv('data/collapsed_gencode_v19_hg19.gz', compression='gzip', sep='\t', encoding='utf-8')
gene='ENSG00000174502'
chrom=1
startbp=205500000
endbp=206000000
gtex_tissues = ['Pancreas']
tissue = gtex_tissues[0]
gwas_data = gwas_data.loc[ (gwas_data[chromcol] == chrom) & (gwas_data[poscol] >= startbp) & (gwas_data[poscol] <= endbp) ]
pops = 'EUR'
gwas_data = gwas_data[[ chromcol, poscol, snpcol, pcol ]]
gwas_data.dropna(inplace=True)
snp_list = list(gwas_data[snpcol])
snp_list = [asnp.split(';')[0] for asnp in snp_list] # cleaning up the SNP names a bit
positions = list(gwas_data[poscol])

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
    gtex_data[tissue] = get_gtex_data(tissue, gene, snp_list, positions, raiseErrors=True)
data.update(gtex_data)
# Obtain any genes to be plotted in the region:
print('Summarizing genes to be plotted in this region')
genes_to_draw = collapsed_genes_df.loc[ (collapsed_genes_df['chrom'] == ('chr' + str(chrom).replace('23','X'))) &
                                        ( ((collapsed_genes_df['txStart'] >= startbp) & (collapsed_genes_df['txStart'] <= endbp)) | 
                                            ((collapsed_genes_df['txEnd'] >= startbp  ) & (collapsed_genes_df['txEnd'] <= endbp  )) ) ]
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






engine = sa.create_engine('mysql://root:root@localhost:3306/snp151')
chunks = pd.read_csv('data/snp151.gz', compression='gzip', sep="\t", chunksize=10000)
for chunk in chunks:
    chunk = chunk.rename(columns={'#chrom':'chrom', 'chromEnd':'position'}).set_index(['name','chrom','position'])
    chunk.head()
    chunk.to_sql(name="snp151_hg19", if_exists='append', index=True, index_label=['name','chrom','position'], 
                        con=engine, dtype={'name': String(50), 'chrom': String(100), 'position': Integer})
    engine.execute("ALTER TABLE `snp151`.`snp151_hg19` add primary key(name(50));")

############### Testing overlaps of genes #################
genes = pd.read_csv('data/gencode_v19_hg19.gz', compression='gzip', sep="\t")
transcript_mapping_file = pd.read_csv('data/ENST_ENSG_mapping_file.txt', sep='\t')
chrom=1
startbp=205500000
endbp=206000000
genes.loc[refseq['#chrom'] == ('chr'+str(chrom))].head()
# chr_groups = refseq.groupby('#chrom',sort=False, group_keys=False)
# sorted_grouped_refseq = refseq.apply(lambda x: x.sort_values(['#chrom' 'txStart','name2'])).groupby('name2')
transcripts = [gene.split('.')[0] for gene in genes['#name']]
rowIndices = [ list(transcript_mapping_file['Transcript stable ID']).index(transcript) for transcript in tqdm(transcripts) ]
ENSG_list = list(transcript_mapping_file.iloc[ rowIndices,: ]['Gene stable ID'])
genes['name3'] = ENSG_list
genes.to_csv('data/GencodeGenes_v19_hg19.gz', compression='gzip', sep="\t", encoding='utf-8', index=False)



import requests

headers = {
    'Content-Type': 'application/json',
}

params = (
    ('token', '6e3b9da087ac'),
)

data = '{"snps": "rs3\nrs4", "pop":"CEU","r2_d": "d"}'

requests.post('https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix', headers=headers, params=params, data=json.dumps(data))



from pymongo import MongoClient
import pandas as pd
conn = "mongodb+srv://naimesca:YF2izjv!@cluster0-p1stg.mongodb.net/test?retryWrites=true"
client = MongoClient(conn)
db = client.GWAS_QTL_app_db
chromLengths_collection = db.hg19_chromLengths.find()
chromLengths = pd.DataFrame(list(chromLengths_collection)).drop(["_id"], axis=1)
chromLengths.set_index('sequence',inplace=True)
print(chromLengths)

