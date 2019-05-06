#!/bin/bash
ERROR_FILE="/home/ec2-user/GWAS-QTL-Explore/gunicorn_error.log"
LOG_FILE="/home/ec2-user/GWAS-QTL-Explore/gunicorn_log.log"
if [[ -e ${ERROR_FILE} ]]; then sudo rm -f ${ERROR_FILE}; fi
if [[ -e ${LOG_FILE} ]]; then sudo rm -f ${LOG_FILE}; fi
sudo /usr/local/bin/gunicorn --timeout 6000 app:app --bind=127.0.0.1:8000 --name GWAS-QTL-Explore --daemon --error-logfile $ERROR_FILE, --log-file $LOG_FILE --workers=3

sudo /etc/init.d/nginx start
