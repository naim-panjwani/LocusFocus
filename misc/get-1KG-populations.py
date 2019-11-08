import pandas as pd

url = 'http://www.internationalgenome.org/category/population/'
populations = pd.read_html(url)[0]
populations = populations[['Population Code', 'Population Description', 'Super Population Code']]

populations.to_csv('data/populations.tsv', sep='\t', index=False, encoding='utf-8')
