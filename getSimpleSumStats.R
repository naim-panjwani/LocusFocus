# Adapted from Fan Wang
# Script to obtain the simple sum P-values for a given set of GWAS p-values, and eQTL p-values for each tissue/gene pair
# Inputs: P_values_filename (GWAS p-values - for a set of SNPs - tab-separated, and all in one line)
#         ld_matrix_filename (the LD matrix filename for the set of SNPs input; the values per row must be tab-separated)
# Ouput: Returns the simple sum P-value
# Example: getSimpleSumStat.R P_values_filename ld_matrix_filename

options(warn=-1)

library(argparser, quietly=TRUE)
library(CompQuadForm, quietly=TRUE)
library(clusterGeneration, quietly=TRUE)
library(psych, quietly=TRUE)
library(data.table, quietly=TRUE)

p <- arg_parser("Calculate Simple Sum Statistic")
p <- add_argument(p, "P_values_filename", help = paste0("Filename with GWAS and eQTL p-values - for a set of SNPs, each value tab-separated'\n'"
                                                        ,"with 1st line being the GWAS p-values'\n'"
                                                        ,"and each subsequent line is for eQTL p-values for each tissue/gene combination"))
p <- add_argument(p, "ld_matrix_filename", help = paste0("The LD matrix filename for the set of SNPs input;'\n'"
                                                         ,"the values per row must be tab-separated; no header"))
argv <- parse_args(p)

Pmat <- fread(argv$P_values_filename, header=F, stringsAsFactors=F, na.strings=c("NaN","nan","NA","-1"), sep="\t")
ldmat <- fread(argv$ld_matrix_filename, header=F, stringsAsFactors=F, na.strings=c("NaN","nan","NA","-1"), sep="\t")
#filename = 'C:\\Users\\Naim\\Desktop\\SCS_Data_Analytics\\homework\\24-Final_project\\GWAS-QTL-Explore\\static\\session_data/Pvalues-2ee9162a-aca8-4768-9b27-73a1fb0514bd.txt'
#Pmat <- fread(filename, header=F, stringsAsFactors=F, na.strings=c("NaN","nan","NA","-1"), sep="\t")
#filename = 'C:\\Users\\Naim\\Desktop\\SCS_Data_Analytics\\homework\\24-Final_project\\GWAS-QTL-Explore\\static\\session_data/ldmat-2ee9162a-aca8-4768-9b27-73a1fb0514bd.txt'
#ldmat <- fread(filename, header=F, stringsAsFactors=F, na.strings=c("NaN","nan","NA","-1"), sep="\t")
Pmat <- as.matrix(Pmat)
P_gwas <- Pmat[1,]
P_eqtl <- Pmat[2:nrow(Pmat),]
ldmat <- as.matrix(ldmat)


get_eqtl_evid<-function(P,cut,m){
  if (cut==0){
  covariate=-log10(P)
  }else{
    covariate=rep(1,m)
    for (i in c(1:m)){
      if (P[i]<cut){
        covariate[i]=1
      }else{
        covariate[i]=0
      }
    }
  }
  return(covariate)
}
  
get_simplesumstats<-function(Zsq,eqtl_evid,m){
  if (sum(eqtl_evid)==0){
    SS<--mean(eqtl_evid)
  }else if (sum(eqtl_evid)==m){
    SS<-mean(eqtl_evid)
  }else{
    reg<-lm(Zsq~eqtl_evid)
    SS<-summary(reg)$coefficients[2,1]
  }
  return(SS)
}
get_A<-function(eqtl_evid,m){
  A=matrix(0,m,m)
  a_diag<-rep(0,m)
  for (j in 1:m){
    if (sum(eqtl_evid)==0){
      a_diag[j]<--1/m
    }else if (sum(eqtl_evid)==m){
      a_diag[j]<-1/m
    }else{
      Tbar<-mean(eqtl_evid)
      a_diag[j]=(eqtl_evid[j]-Tbar)/(sum(eqtl_evid^2)-m*(Tbar^2))
    }
  }
  diag(A)=a_diag
  return(A)
}

get_eigenvalue<-function(eqtl_evid,ld.mat,m){
  chol_Sigma=chol(ld.mat)
  result=list()
  matrix_A<-get_A(eqtl_evid,m)
  matrix_mid<-chol_Sigma%*%matrix_A%*%t(chol_Sigma)
  eigenvalues<-eigen(matrix_mid)$values
  return(eigenvalues)
}

get_p<-function(m,eigen_value,teststats){
  l=length(eigen_value)
  pv<-abs(imhof(teststats,eigen_value,h=rep(1,l),delta=rep(0,l))$Qq)
  return(pv)
}

get_simplesumP<-function(P_gwas,P_eqtl,ld.mat,cut,m){
  ##need to match the GWAS SNP with the eQTL SNP and get m
  Z=qnorm(P_gwas/2)
  Zsq=Z^2
  ##get eqtl evidence
  eqtl_evid=get_eqtl_evid(P_eqtl,cut,m)
  #get Simple Sum statistic
  SS_stats<-get_simplesumstats(Zsq,eqtl_evid,m)
  ##get eigenvalues:
  eig_values<-get_eigenvalue(eqtl_evid,ld.mat,m)
  ##get Simple Sum p-values
  pv<-get_p(m,eig_values,SS_stats)
  return(pv)
}

# 
# set.seed(1)
# m=100
# #generate LD matrix 
# covariance=genPositiveDefMat("eigen",dim=m)$Sigma
# var=diag(covariance);var
# Sigma=solve(sqrt(diag(var)))%*%covariance%*%solve(sqrt(diag(var)))
# #generate GWAS summary statistic
# mu_Z=rep(0,m)
# mu_T=rep(0,m)
# Z_resp<-mvrnorm(n=1,mu_Z,Sigma)
# Wald=Z_resp^2
# #generate eQTL summary statistic
# covar_T<-mvrnorm(n=1,mu_T,Sigma)
# #transform eQTL evidence
# #if cut=0, eQTL evidence would be -log10 transform of eQTL p-value;
# #if cut <0 (i.e. cut=0.05),eQTL evidence would be dischotomized eQTL p-value indicator by thresholds of eQTL p<cut.
# cut=0
# P=2*pnorm(abs(covar_T),lower.tail = F)
# #get_simplesumP(Z_resp,P,Sigma,0,m) -- modified by Naim 28-Apr to take GWAS P instead of Z_resp
# ##should equal to 0.474327
# 
# P_gwas = 2*pnorm(-abs(Z_resp))
# get_simplesumP(P_gwas,P,Sigma,0,m) # input P_gwas instead of Z_resp
# 


############
# MAIN
############

Pss <- NULL
num_iterations = 0
if (is.null(nrow(P_eqtl)) & length(P_eqtl)>0) {
  num_iterations=1
} else if(!is.null(nrow(P_eqtl))) {
  num_iterations = nrow(P_eqtl)
} else {
  stop("No eQTL P-values provided")
}
for(i in 1:num_iterations) {
  tempmat <- cbind(P_gwas, P_eqtl[i,])
  NArows = which(is.na(tempmat[,1]) | is.na(tempmat[,2]))
  tempmat = tempmat[-NArows,]
  if(nrow(tempmat) != 0) {
    m <- nrow(tempmat)
    P = get_simplesumP(P_gwas=as.numeric(tempmat[,1]), P_eqtl=as.numeric(tempmat[,2]), ld.mat=ldmat[-NArows, -NArows], cut=0, m=m)
    Pss = c(Pss, P)
  } else {
    Pss = c(Pss, -1)
  }
}

write(Pss, stdout())
