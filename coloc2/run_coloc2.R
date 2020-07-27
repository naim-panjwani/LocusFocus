# Script to obtain the coloc2 PP4 posterior probabilities for a given primary summary statistic dataset (GWAS), 
#  and one or more secondary datasets (eQTLs for tissue/gene pairs)
# Inputs: File name of primary dataset of summary statistics (e.g. GWAS p-values).
#         File name containing secondary dataset(s) of summary statistics
# Ouput: Returns a data.frame with the ProbeID and corresponding PP4 posterior probabilities
# Example: run_coloc2.R primary_dataset_filename secondary_dataset_filename --outfilename output_filename

# Inputs details:
## Primary dataset (GWAS) required columns: CHR, POS, A1 (ALT), A2 (REF), BETA, SE, PVAL, N, MAF, type
## Secondary dataset (eQTL) required columns: CHR, POS, A1 (ALT), A2 (REF), BETA, SE, PVAL, N, MAF, ProbeID
###  ** MAF is used by COLOC2 to filter rare variants with MAF below 0.001 (0.1%) **

options(warn=-1)

#library('doParallel')
library('data.table', quietly = TRUE)
#library('foreach')
#setwd('/mnt/c/Users/Naim/Desktop/coloc2/original_coloc2/coloc2')
library(argparser, quietly = TRUE)
library(here, quietly = TRUE)
source(here('coloc2','functions_coloc_likelihood_summary_integrated2.R'))
source(here("coloc2", "optim_function.R"))


######
# Parse arguments
######

p <- arg_parser("Calculate COLOC2 PP4 Posterior Probabilities")
p <- add_argument(p, "gwas_filename", help = paste0("Filename with GWAS and eQTL p-values - for a set of SNPs, each value tab-separated'\n'"
                                                        ,"with 1st line being the GWAS p-values'\n'"
                                                        ,"and each subsequent line is for eQTL p-values for each tissue/gene combination"))
p <- add_argument(p, "eqtl_filename", help = paste0("The LD matrix filename for the set of SNPs input;'\n'"
                                                         ,"the values per row must be tab-separated; no header"))
p <- add_argument(p, "--outfilename", default = 'Coloc2Pvalues.txt', help = "Output filename")
argv <- parse_args(p)
gwas_filename <- argv$gwas_filename
eqtl_filename <- argv$eqtl_filename
outfilename <- argv$outfilename



##perform COLOC2 analysis for one gene and tissue
# gwas_filename <- 'gwas_file.txt'
# eqtl_filename <- 'eqtl_file.txt'
# GWAS_file=read.table(file=here('coloc2', 'sample_data', gwas_filename),header=T,stringsAsFactors = F)
# # GWAS_file = GWAS_file[, c("PVAL","type","CHR","POS","MAF","N","SNPID")]
# eQTL_file=read.table(file=here('coloc2', 'sample_data', eqtl_filename),header=T,stringsAsFactors = F)
# # eQTL_file = eQTL_file[, c("PVAL","ProbeID","CHR","POS","MAF","N","SNPID")]
# res1=coloc.eqtl.biom(eqtl.df=eQTL_file,biom.df=GWAS_file,p12=1e-6,useBETA=TRUE,outfolder=getwd(),prefix="output", plot=FALSE, save.coloc.output=FALSE, match_snpid=T)
# # get colocalization posterior probability:0.9481824
# res1$PP.H4.abf

##perform COLOC2 analysis for 40 genes
#GWAS_file2=read.table(file=here('coloc2', 'sample_data', 'gwas_file_for40genes.txt'),header=T,stringsAsFactors = F)
#eQTL_file2=read.table(file='sample_data/eqtl_file_40genes.txt',header=T,stringsAsFactors = F)
#res2=coloc.eqtl.biom(eqtl.df=eQTL_file2,biom.df=GWAS_file2,p12=1e-6,useBETA=TRUE,outfolder=getwd(),prefix="output", plot=FALSE, save.coloc.output=FALSE, match_snpid=T)
#res2[,c("ProbeID","PP.H4.abf")]

#gwas_filename <- "static/session_data/coloc2gwas_df-b05aef08-b42c-4f8d-a824-8cb76c8c2f4c.txt"
#eqtl_filename <- "static/session_data/coloc2eqtl_df-b05aef08-b42c-4f8d-a824-8cb76c8c2f4c.txt"

if(!file.exists(gwas_filename)) stop(paste0("COLOC2: Did not find GWAS file: ", gwas_filename))
if(!file.exists(eqtl_filename)) stop(paste0("COLOC2: Did not find eQTL file: ", eqtl_filename))
GWAS_file <- read.table(file=gwas_filename,header=T,stringsAsFactors = F)
colnames(GWAS_file)[which(colnames(GWAS_file) %in% "REF")] <- "A2"
colnames(GWAS_file)[which(colnames(GWAS_file) %in% "ALT")] <- "A1"
if(!("SNPID" %in% colnames(GWAS_file))) GWAS_file <- cbind(GWAS_file, SNPID=with(paste0(gsub("chr","",CHR),"_",POS,"_",A2,"_",A1),data=GWAS_file))
eQTL_file <- read.table(file=eqtl_filename,header=T,stringsAsFactors = F)
colnames(eQTL_file)[which(colnames(eQTL_file) %in% "REF")] <- "A2"
colnames(eQTL_file)[which(colnames(eQTL_file) %in% "ALT")] <- "A1"
if(!("SNPID" %in% colnames(eQTL_file))) eQTL_file <- cbind(eQTL_file, SNPID=with(paste0(gsub("chr","",CHR),"_",POS,"_",A2,"_",A1),data=eQTL_file))

res <- coloc.eqtl.biom(eqtl.df=eQTL_file,biom.df=GWAS_file,p12=1e-6,useBETA=TRUE,outfolder=getwd(),prefix="output", plot=FALSE, save.coloc.output=FALSE, match_snpid=T)
result <- res[,c("ProbeID","PP.H4.abf")]
colnames(result) <- c("ProbeID","PPH4abf")

write.table(result, outfilename, row.names=F, col.names = T, quote = F, sep="\t")
