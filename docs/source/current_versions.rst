.. _current_versions:

##################
Current versions
##################

***************************
Version History
***************************

- v0.0.1 (released Sept. 6, 2019)
- v0.0.2 (released Nov. 8, 2019)
   - Fixed file upload issue; up to 3 files may be uploaded at once and file types are auto-detected
   - Increased file size limit to 100 MB
   - Added ability to change the eQTL gene within the plot output page
   - Added interactive table of Simple Sum -log10 P-value results, with ability to download the table in various formats
   - Added t\ :sup:`2` test for testing whether the given secondary datasets have a significant signal above a Bonferroni adjustment
   - Added ability to upload SNP names in chr_pos_ref_alt_b37 format for cases when the rs ID is not available
- v1.0 alpha (released Dec. 5, 2019)
   - Enabled ability to upload secondary datasets as a merged HTML file in addition to any selected GTEx tissues
- v1.0.1 alpha (released Dec. 18, 2019)
   - Added ability to transpose table and adjust plot figure drawing parameters  
- v1.1.0 alpha (released Apr. 1, 2020)
   - Internal change to calculate the Simple Sum using an R script instead of Python. 
   - This change enables the use of the app in Windows as the rpy2 package is no longer a requirement
- v1.3.0 alpha (released Jul. 31, 2020)
   - Addition of hg38 coordinate support
   - Added the latest GTEx version 8 (hg38) eQTL analyses for use as secondary datasets for colocalization testing
   - Added GRCh38 re-aligned 1000 Genomes (phase 3) as option for LD matrix
   - Using GENCODE v26 for hg38 gene track
   - Added support for COLOC2 colocalization testing
- v1.4.0 alpha (released Aug. 6, 2020)
   - Added ability to export images as svg vector format
   - Bug fix for the merge_and_convert_to_html_coloc2.py script
   - More sample datasets added that are compatible for COLOC2 runs
- v1.4.1 alpha (Oct. 6, 2020)
   - Bug fix in identifying lead SNP
- v1.4.2 alpha (Oct. 9, 2020)
   - Fixed issue where SS was not being computed due to non-singular matrix error
   - Improved initial plotting time of colocalization plot by filtering out GWAS p-values less than 0.1 by default
- v1.4.3 alpha (Nov. 10, 2020)
   - Fixed issue where the Simple Sum calculation was not being performed in the case where no missing data was present
- v1.4.4 alpha (Nov. 16, 2020)
   - Fixed rsid mapping bug
   - Fixed bug when matching the top SNP with secondary datasets
- v1.4.5 alpha (Feb. 04, 2020)
   - Fixed issue where the SS window did not follow the user-specified lead SNP
   - Added a function to clean up uploaded dataset of common file input reading issues


******************
Datasets
******************

- Human reference: hg19 (GRCh37.p13) and hg38 (GRCh38.p7)
- `GTEx: versions 7 (hg19) and 8 (hg38) <https://gtexportal.org/home/>`_
- GENCODE: `version 19 (hg19) <https://github.com/naim-panjwani/LocusFocus/blob/master/data/collapsed_gencode_v19_hg19.gz>`_ (the transcript models were collapsed into a single gene model)
- GENCODE: `version 26 (hg38) <https://github.com/naim-panjwani/LocusFocus/blob/master/data/collapsed_gencode_v26_hg38.gz>`_ (the transcript models were collapsed into a single gene model)
- `1000 Genomes (phase 3) aligned to GRCh37 biallelic SNV call set  <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>`_
- `1000 Genomes (phase 3) biallelic SNV call set re-aligned to GRCh38 <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/>`_
- `dbSNP151 GRCh37.p13 <ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13>`_
- `dbSNP151 GRCh38.p7 <ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/>`_

******************
Programs
******************

All required programs and versions are specified in the `yml file <https://github.com/naim-panjwani/LocusFocus/blob/master/environment.yml>`_ 
or `conda spec file <https://github.com/naim-panjwani/LocusFocus/blob/master/spec-file.txt>`_.

