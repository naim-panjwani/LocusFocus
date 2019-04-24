library('CompQuadForm')
library('clusterGeneration')
library('psych')

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

get_simplesumP<-function(Z,P,ld.mat,cut,m){
  ##need to match the GWAS SNP with the eQTL SNP and get m
  Zsq=Z^2
  ##get eqtl evidence
  eqtl_evid=get_eqtl_evid(P,cut,m)
  #get Simple Sum statistic
  SS_stats<-get_simplesumstats(Zsq,eqtl_evid,m)
  ##get eigenvalues:
  eig_values<-get_eigenvalue(eqtl_evid,ld.mat,m)
  ##get Simple Sum p-values
  pv<-get_p(m,eig_values,SS_stats)
  return(pv)
}


set.seed(1)
m=100
#generate LD matrix 
covariance=genPositiveDefMat("eigen",dim=m)$Sigma
var=diag(covariance);var
Sigma=solve(sqrt(diag(var)))%*%covariance%*%solve(sqrt(diag(var)))
#generate GWAS summary statistic
mu_Z=rep(0,m)
mu_T=rep(0,m)
Z_resp<-mvrnorm(n=1,mu_Z,Sigma)
Wald=Z_resp^2
#generate eQTL summary statistic
covar_T<-mvrnorm(n=1,mu_T,Sigma)
#transform eQTL evidence
#if cut=0, eQTL evidence would be -log10 transform of eQTL p-value;
#if cut <0 (i.e. cut=0.05),eQTL evidence would be dischotomized eQTL p-value indicator by thresholds of eQTL p<cut.
cut=0
P=2*pnorm(abs(covar_T),lower.tail = F)
get_simplesumP(Z_resp,P,Sigma,0,m)
##should equal to 0.474327
