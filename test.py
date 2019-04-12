import pandas as pd
data = pd.read_csv('data/populations.tsv',sep="\t", encoding='utf-8')
data['lead_snp'] = 'rs662702'
subdata = data[['Population Code','lead_snp']]
subdata.head()

test1 = "1"
print(test1.replace('chr',''))
test2 = 'chr1'
print(test2.replace('chr',''))

chroms = pd.read_csv('data/hg19_chrom_lengths.txt', sep="\t", encoding='utf-8')
chroms.head()
chroms[1]
chroms[1] = [int(chr.strip().replace(',','')) for chr in chroms[1]]
chroms.head()
chroms[0] = [chr.strip() for chr in chroms[0]]

chroms.to_csv('data/hg19_chrom_lengths.txt', index=False, encoding='utf-8', sep="\t")

gwas_data = pd.read_csv('data/MI_GWAS_2019_1_205500-206000kbp.tsv', sep="\t")
#gwas_data = gwas_data[['SNP','P']]
gwas_data.dropna(inplace=True)
gwas_data.head()
snpcol='SNP'
pcol='P'
lead_snp = list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ]['SNP'])[0]
snp_list = list(gwas_data[snpcol])
positions = list(gwas_data['BP'])


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
