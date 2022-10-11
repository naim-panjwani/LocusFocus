#library('doParallel')
library('data.table')
#library('foreach')
#setwd('/mnt/c/Users/Naim/Desktop/coloc2/original_coloc2/coloc2')
source('functions_coloc_likelihood_summary_integrated2.R')

##perform COLOC2 analysis for one gene and tissue
GWAS_file=read.table(file='sample_data/gwas_file.txt',header=T,stringsAsFactors = F)
# GWAS_file = GWAS_file[, c("CHR","POS","BETA","SE","N","PVAL","MAF","type")]
eQTL_file=read.table(file='sample_data/eqtl_file.txt',header=T,stringsAsFactors = F)
# eQTL_file = eQTL_file[, c("CHR","POS","BETA","SE","N","PVAL","MAF","ProbeID")]
# GWAS_file <- cbind(GWAS_file, A1="A",A2="T")
# GWAS_file <- cbind(GWAS_file, SNPID=with(paste0(gsub("chr","",CHR),"_",POS,"_",A2,"_",A1),data=GWAS_file))
# eQTL_file <- cbind(eQTL_file, A1="A",A2="T")
# eQTL_file <- cbind(eQTL_file, SNPID=with(paste0(gsub("chr","",CHR),"_",POS,"_",A2,"_",A1),data=eQTL_file))
res1=coloc.eqtl.biom(eqtl.df=eQTL_file,biom.df=GWAS_file,p12=1e-6,useBETA=TRUE,outfolder=getwd(),prefix="output", plot=FALSE, save.coloc.output=FALSE, match_snpid=T)
##get colocalization posterior probability:0.9481824
res1$PP.H4.abf

##perform COLOC2 analysis for 40 genes
GWAS_file2=read.table(file='sample_data/gwas_file_for40genes.txt',header=T,stringsAsFactors = F)
eQTL_file2=read.table(file='sample_data/eqtl_file_40genes.txt',header=T,stringsAsFactors = F)
res2=coloc.eqtl.biom(eqtl.df=eQTL_file2,biom.df=GWAS_file2,p12=1e-6,useBETA=FALSE,outfolder=getwd(),prefix="output", plot=FALSE, save.coloc.output=FALSE, match_snpid=T)
res2
