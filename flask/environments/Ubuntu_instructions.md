Here are the steps to follow to end up with an environment with a local LocusFocus web server in Ubuntu 20.04.1 LTS.

```bash
conda create --name locusfocusR --file spec-file-old.txt
conda activate locusfocusR
python -m pip install lxml
python -m pip install blinker
python -m pip install flask_sitemap
python -m pip install flask-uploads
```

This will allow one to use LocusFocus to test colocalization if all 3 files are available (primary GWAS, secondary colocalization dataset, and LD matrix) by simply turning on the webserver and issuing the command:

```bash
python app.py
```

and then navigating to the local URL given (usually http://127.0.0.1:5000/).

