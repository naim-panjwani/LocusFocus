# LocusFocus
## Explore GWAS Results with GTEx's eQTL Association Results 

This web application is to enable the exploration and annotation of GWAS results to assess what gene and tissue are the most probable candidates for the GWAS signal via integration with [GTEx eQTL association results](https://gtexportal.org/home/) into a single interactive plot, and application of the [Simple Sum method](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008007) to visualize the results in a heatmap plot.

LD information is calculated using PLINK on the 1000 Genomes Project to enable viewing of LD given a lead SNP. The default lead SNP is the top SNP.

eQTL association results have been stored locally in a MongoDB database. 

Full documentation is available at https://locusfocus.readthedocs.io.
The application has been published in SickKids servers and may be accessed at https://locusfocus.research.sickkids.ca.
