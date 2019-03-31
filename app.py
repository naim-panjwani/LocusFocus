import requests
import json
import pandas as pd
import numpy as np
import tokens

from flask import Flask, jsonify, render_template
from flask_sqlalchemy import SQLAlchemy
import pymysql
pymysql.install_as_MySQLdb()
#thepwd = open('pwd.txt').readline().replace('\n', '')

app = Flask(__name__)
genomicWindowLimit = 2000000
fileSizeLimit = 200e6


####################################
# LD Querying
####################################

def parseR2(urltext):
    temp = urltext.strip().split('\n')
    for t in temp:
        if t.strip().find('R2') != -1:
            return(float(t.strip().split(':')[1].strip()))
    return(np.nan)

lead_snp = 'rs7512462'
snp_list = ['rs61814953', 'rs1342063']
european_populations = ['CEU', 'TSI', 'FIN', 'GBR', 'IBS']
ld_type = 'r2'
base_url = 'https://ldlink.nci.nih.gov/LDlinkRest/ldpair?'

r2_pairs = []
for snp in snp_list:
    snp_query = 'var1=' + lead_snp + '&var2=' + snp
    population_query = '&pop=' + "%2B".join(european_populations) 
    #ld_query = '&r2_d=' + ld_type 
    token_query = '&token=' + tokens.token
    url = base_url + snp_query + population_query + token_query
    response = requests.get(url).text
    r2 = parseR2(response)
    #print(r2)
    r2_pairs.append(r2)
    

#####################################
# API Routes
#####################################

@app.route("/")
def index():
    """Return the homepage."""
    return render_template("index.html")

@app.route("/populations")
def get1KGPopulations():
    populations = pd.read_csv('data/populations.tsv', sep='\t')
    return jsonify(populations.to_dict(orient='list'))



if __name__ == "__main__":
    app.run()

