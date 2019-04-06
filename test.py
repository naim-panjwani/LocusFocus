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
