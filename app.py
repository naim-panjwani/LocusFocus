import json
import requests
import pandas as pd
import numpy as np
import tokens
import os

from flask import Flask, request, redirect, url_for, jsonify, render_template
from flask_sqlalchemy import SQLAlchemy
import pymysql
pymysql.install_as_MySQLdb()
#thepwd = open('pwd.txt').readline().replace('\n', '')

genomicWindowLimit = 2000000
fileSizeLimit = 200e6

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'upload'

####################################
# LD Querying
####################################

def parseR2(urltext):
    temp = urltext.strip().split('\n')
    for t in temp:
        if t.strip().find('R2') != -1:
            return float(t.strip().split(':')[1].strip())
    return np.nan

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

# @app.route("/")
# def index():
#     """Return the homepage."""
#     return render_template("index.html")

@app.route("/populations")
def get1KGPopulations():
    populations = pd.read_csv('data/populations.tsv', sep='\t')
    return jsonify(populations.to_dict(orient='list'))

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    data = {"success": False}
    if request.method == 'POST':
        if request.files.get('file'):
            # read the file
            file = request.files['file']

            # read the filename
            filename = file.filename

            # create a path to the uploads folder
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)

            file.save(filepath)

            # Load
            data = pd.read_csv(filepath, sep="\t", encoding='utf-8')
            #data.head()

            snpcol = request.form['snp-col']
            pcol = request.form['pval-col']
            lead_snp = request.form['leadsnp']
            regiontext = request.form['locus']

            data['lead_snp'] = lead_snp

            # indicate that the request was a success
            success = True

        return jsonify(data.to_dict(orient='list'))

    return render_template("index.html")
    

if __name__ == "__main__":
    app.run()

