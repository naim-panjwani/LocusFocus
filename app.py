import json
import requests
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

@app.route("/")
def index():
    """Return the homepage."""
    return render_template("index.html")

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

            # Load the saved image using Keras and resize it to the Xception
            # format of 299x299 pixels
            image_size = (299, 299)
            im = keras.preprocessing.image.load_img(filepath,
                                                    target_size=image_size,
                                                    grayscale=False)

            # preprocess the image and prepare it for classification
            image = prepare_image(im)

            global graph
            with graph.as_default():
                preds = model.predict(image)
                results = decode_predictions(preds)
                # print the results
                print(results)

                data["predictions"] = []

                # loop over the results and add them to the list of
                # returned predictions
                for (imagenetID, label, prob) in results[0]:
                    r = {"label": label, "probability": float(prob)}
                    data["predictions"].append(r)

                # indicate that the request was a success
                data["success"] = True

        return jsonify(data)

    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    '''




if __name__ == "__main__":
    app.run()

