#!/bin/bash

# Setup conda environment:
conda create --name locusfocusR --file environments/spec-file-old.txt
conda activate locusfocusR
python -m pip install lxml
python -m pip install blinker
python -m pip install flask_sitemap
python -m pip install flask-uploads
python -m pip install turicreate

# Create missing folders:


# Download dbSNP151:


# Download 1000 Genomes data for LD - optional:


## Filter for biallelic SNVs and convert 1000 Genomes data to binary PLINK format:



# Install mongodb - optional:
curl -fsSL https://www.mongodb.org/static/pgp/server-4.4.asc | sudo apt-key add -
echo "deb [ arch=amd64,arm64 ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/4.4 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-4.4.list
sudo apt update
sudo apt install mongodb-org
sudo systemctl start mongod.service
sudo systemctl enable mongod

# Create GTEx collections - optional:
## Download and fix the files:


## Push data to a MongoDB database
python misc/initdb_GTExV8.py


