# Explore GWAS Results with GTEx's eQTL Association Results 
This app is to enable exploration of GWAS results and assessment of what gene and tissue are the most likely responsible candidates for the GWAS signal via integration with [GTEx eQTL association results](https://gtexportal.org/home/) (v6p and v7) into a single interactive plot, and application of the [Simple Sum method](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008007) to visualize the results in a heatmap plot.

LD information is obtained from [LDlink's API](https://ldlink.nci.nih.gov/?tab=apiaccess) to enable viewing of LD given a lead SNP. The default lead SNP is the top SNP.

eQTL association results are obtained from [GTEx's API](https://gtexportal.org/home/api-docs/). 

You may upload your GWAS results to the app, and specify the region of interest (up to 2MB). You may also provide a custom region to apply the Simple Sum calculation to (default is 100kb from the lead SNP).

The gene track shown is the RefSeq gene track on the GRCh37 human genome build (hg19).

