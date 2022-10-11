#!/bin/bash
conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda create -n locusfocusR python r-base spyder pandas numpy scipy bs4 lxml flask pymongo python-dateutil tqdm libuuid dask pysam sphinx Werkzeug=0.16.1 r-base r-r.utils rstudio r-argparser r-compquadform r-data.table 
pip install blinker==1.4
pip install flask_sitemap==0.3.0
pip install flask-uploads==0.2.1
