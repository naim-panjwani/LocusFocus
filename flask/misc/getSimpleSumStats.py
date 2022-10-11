# Adapted from Fan Wang, coverted to Python by Scott Mastromatteo
# Script to obtain the simple sum P-values for a given set of GWAS p-values, and eQTL p-values for each tissue/gene pair
# Inputs: P_values_filename (GWAS p-values - for a set of SNPs - tab-separated, and all in one line)
#         ld_matrix_filename (the LD matrix filename for the set of SNPs input; the values per row must be tab-separated)
# Ouput: Returns the simple sum P-value, number of SNPs used in each calculation
# Example: getSimpleSumStat.py p_values_filename ld_matrix_filename


import argparse
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm
import statsmodels.formula.api as stats
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

imhof = importr(str("CompQuadForm")).__dict__['imhof']
davies = importr(str("CompQuadForm")).__dict__['davies']

def set_based_test(summary_stats, ld_, num_genes, alpha=0.05):
    z = norm.ppf(summary_stats/2)
    zsq = np.square(z)
    statistic = sum(zsq)
    m = len(zsq)
    eigenvalues = np.linalg.eigvals(ld_).real
    res = imhof(statistic, eigenvalues)
    pv = abs(np.array(res[res.names.index('Qq')]))
    if pv < (alpha / num_genes):
        return True
    else:
        return False


def get_p(ss_stat, eig_values, m, meth='davies'):

    if meth == 'davies':
        res = davies(ss_stat, eig_values)
    elif meth == 'imhof':
        res = imhof(ss_stat, eig_values)
    p = np.array(res[res.names.index('Qq')])

    return abs(p[0])

  
def get_a_diag(eqtl_evid, m):
    s = np.sum(eqtl_evid)

    if s == 0 or s == m:
        a_diag = np.repeat(1.0/m, m)
    else:
        t_bar = np.mean(eqtl_evid)
        denom = np.sum(np.square(eqtl_evid)) - m*(t_bar*t_bar)
      
        a_diag = [(x - t_bar)/denom for x in eqtl_evid]
      
    return a_diag


def get_eigenvalues(eqtl_evid, ld_, m):
        
    chol_sigma = scipy.linalg.cholesky(ld_)
    a_diag = get_a_diag(eqtl_evid, m)
    
    a = np.zeros((m, m)) 
    np.fill_diagonal(a, a_diag)
    mid = np.linalg.multi_dot([chol_sigma, a, np.transpose(chol_sigma)])
    
    return np.linalg.eigvals(mid).real


def get_simple_sum_stats(zsq, eqtl_evid, m):
    s = np.sum(eqtl_evid)
    if s == 0 or s == m:
        return np.mean(eqtl_evid)
    
    df = pd.DataFrame({"zsq": zsq, "eqtl": eqtl_evid})
    lm = stats.ols(formula="zsq ~ eqtl", data=df).fit()
    return lm.params[1]


##if cut = 0, eQTL evidence would be -log10 transform of eQTL p-value;
##if cut < 0 (i.e. cut=0.05), eQTL evidence would be dischotomized eQTL p-value indicator by thresholds of eQTL p<cut.
def get_eqtl_evid(p, cut, m):
    if cut == 0:
        covariate = -np.log10(p)
    else:
        mask = p < cut
        covariate = np.array(p, copy=True) 
        covariate[mask] = 1
        covariate[~mask] = 0
    
    return covariate


def simple_sum_p(gwas_, eqtl_, ld_, cut, m, meth='davies'):
    ##need to match the GWAS SNP with the eQTL SNP and get m
    z = norm.ppf(gwas_/2)
    zsq = np.square(z)

    ##get eqtl evidence
    eqtl_evid = get_eqtl_evid(eqtl_, cut, m)
    
    #get Simple Sum statistic
    ss_stat = get_simple_sum_stats(zsq, eqtl_evid, m)
    
    ##get eigenvalues:
    eig_values = get_eigenvalues(eqtl_evid, ld_, m)
    
    ##get Simple Sum p-values
    return get_p(ss_stat, eig_values, m, meth=meth)


def parse_param():
    parser = argparse.ArgumentParser(description="Calculate Simple Sum Statistic")
    parser.add_argument("p_values_filename", metavar="pvalues", 
                        help="Filename with GWAS and eQTL p-values - for a set of SNPs, " + \
                        "each value tab-separated with 1st line being the GWAS p-values")
    parser.add_argument("ld_matrix_filename", metavar="ld",
                        help="The LD matrix filename for the set of SNPs input; " + \
                        "the values per row must be tab-separated; no header")
    args = parser.parse_args()
        
    p_mat = np.loadtxt(args.p_values_filename)
    ld_mat = np.loadtxt(args.ld_matrix_filename)
    return (p_mat, ld_mat)

def parse_param2():

    path = "p_example.txt" 
    p_mat = np.loadtxt(path)
    path = "ld_example.txt" 
    ld_mat = np.loadtxt(path)
    return (p_mat, ld_mat)


def get_simple_sum_p(p_mat, ld_mat): # main function to use
    # P-values returned can be negative and have the following meanings:
    # -1: there was no eQTL data
    # -2: fails the set_based_test(), so eQTL region is not significant after Bonferroni correction
    # -3: could not compute the Simple Sum p-value; this is likely due to insufficient number of SNPs
    
    #p_mat, ld_mat = parse_param2()
    
    if p_mat.ndim == 2 and np.size(p_mat, 0) > 1:
        n_iter = np.size(p_mat, 0) - 1
    else: #not a 2d matrix or nrow(p_eqtl) < 1 
      return None
    
    p_gwas = np.asarray(p_mat[0,])
    p_eqtl = p_mat[1:,]
   
    pretest, pss, n, comp_used = [], [], [], []
    
    for i in range(n_iter):  
        p_eqtl_i = np.asarray(p_eqtl[i,] if n_iter > 1 else p_eqtl)
        
        #remove nans        
        mask = np.isnan(p_gwas) | np.isnan(p_eqtl_i)
        eqtl_ = p_eqtl_i[~mask]
        gwas_ = p_gwas[~mask]
        ld_ = ld_mat[~mask,]
        ld_ = ld_[:,~mask]

        #count SNPs
        snp_count = sum(~mask)
        n.append(snp_count)
        
        if snp_count < 1:
            pss.append(-1) # no eqtl data
            comp_used.append('NA')
            continue
        
        # do pretest (set_based_test)
        try:
            if set_based_test(eqtl_, ld_, n_iter):
                P = simple_sum_p(gwas_, eqtl_, ld_, cut=0, m=snp_count, meth='davies')
                if P == 0 or P<0:
                    P = simple_sum_p(gwas_, eqtl_, ld_, cut=0, m=snp_count, meth='imhof')
                    comp_used.append('imhof')
                    pss.append(P)
                else:
                    comp_used.append('davies')
                    pss.append(P)
            else:
                pss.append(-2) # not significant eQTL given the # of genes
        except:
            pss.append(-3) # could not compute a SS p-value (SNPs not dense enough?)

                
    return (pss, n, comp_used)





