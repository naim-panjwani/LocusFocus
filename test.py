import pandas as pd
data = pd.read_csv('data/populations.tsv',sep="\t", encoding='utf-8')
data['lead_snp'] = 'rs662702'
