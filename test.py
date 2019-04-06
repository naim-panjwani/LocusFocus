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
gwas_data = gwas_data[['SNP','P']]
gwas_data.dropna(inplace=True)
gwas_data.head()
snpcol='SNP'
pcol='P'
list(gwas_data.loc[ gwas_data[pcol] == min(gwas_data[pcol]) ]['SNP'])[0]
snp_list = list(gwas_data[snpcol])

