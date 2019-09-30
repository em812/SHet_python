#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates the multivariate SHet statistic for an individual SNP from 
summary statistics (effect size b and standard error se) obtained with univariate 
mlma analysis in GCTA for all the the tierpsy features of the CeNDR strains.

Created on Fri Sep 13 17:40:07 2019

@author: em812
"""
import numpy as np
from numpy.linalg import pinv
from time import time
import pandas as pd 

def EstimateGamma(N_permut, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = True):
    
    from numpy.random import multivariate_normal as mvrnorm
    from scipy.stats import gamma
    
    # Sample SHet
    Stat = []
    for i in range(N_permut):
        X = mvrnorm(np.zeros(SampleSize.shape[0]), CorrMatrix, check_valid='warn', tol=1e-8)
    
        Stat.append(Trucated_TestScore(X, SampleSize, CorrMatrix, correct=correct,
                                       startCutoff=startCutoff, endCutoff=endCutoff, 
                                       CutoffStep=CutoffStep, isAllpossible=isAllpossible))
    # Fit gamma distribution
    fit_alpha, fit_loc, fit_beta = gamma.fit(Stat)
    k = fit_alpha
    theta = 1/fit_beta
    a = fit_loc
#    a = min(Stat)*3/4
#    ex3 = np.mean(Stat*Stat*Stat)
#    V =	np.var(Stat)
#  
#    for i in np.arange(1,101,1):
#        E = np.mean(Stat)-a
#        k = E**2/V
#        theta = V/E
#        a = (-3*k*(k+1)*theta**2+sqrt(9*k**2*(k+1)**2*theta**4-12*k*theta*(k*(k+1)*(k+2)*theta**3-ex3)))/6/k/theta
  
    para = [k,theta,a]
    return(para)

def Trucated_TestScore(X, SampleSize, CorrMatrix, correct = True, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = False):
    """
    Calculates the SHet statistic.
    params:
        X: Wald summary statistic from each cohort for each trait (for a given SNP)
           Each row represents a cohort and each column represents a trait (number of cohorts x number of traits)
        SampleSize: (size of the cohorts) weights in the metaanalysis (number of cohorts x number of traits)
                    [if only one cohort is assumed, the sample size for each trait is the 
                     same and all the weights are equal]
        CorrMatrix: correlation between phenotypes
        correct: flag defining whether the stats with negative sign will be corrected (multiplied by -1) to add to the evidence
        startCutoff, endCutoff, CutoffStep: define the cutoffs in the wald summary stats values that we want to consider when building SHet
        isAllpossible: if True, all the values in X will be used as cutoffs (only one summary statistic will be dropped at a time)
    return:
        SHet
    """
    
    # Number of traits
    #N = X.shape[1]
    
    # Weights of summary statistics -> sqrt(sample size of cohort) 
    W = 1/np.sqrt(SampleSize)
    
    
    #XX = apply(X, 1, function(x) {
    
    # For each row of X (all traits from a given cohort)
    all_time = time()
    
    # Initialize
    TTT = -1;
    sel_group = None
    
    # Get cutoff levels
    if isAllpossible:
        cutoff = np.sort(np.unique(np.abs(X)));	  ## it will filter out any of them. 
    else:
        cutoff = np.arange(startCutoff, endCutoff+CutoffStep, CutoffStep);		
  
    # Calculate SHet for each cutoff threshold
    for threshold in cutoff:
        ptm = time()
        print("Calculating for threshold ",threshold)
        
        # If all the Wald summary stats are smaller that the threshold, break the for loop 
        # (there aren't any larger stats left to group)
        if np.all(np.abs(X)<threshold):
            break
        
        # Else calculate SHet using the statistics that are >= threshold
        else:
            index = np.abs(X) >= threshold
            x1 = X[index]
            A = CorrMatrix[index, :]   ## update the matrix
            A = A[:, index]
            W1 =  W[index]
      
        
        # Change sign of negative stats to add evidence
        if correct:
            W1[(x1 < 0)] = -W1[(x1 < 0)]    ## update the sign
        
        # Calculate SHet
        try:
            A = pinv(A)
        except:
            pd.DataFrame(A).to_csv('A_non_invertible_threshold={}.csv'.format(threshold))
            print('error')
        x1 = x1.reshape(1,-1) 
        W1 = W1.reshape(1,-1)
        T = np.matmul(np.matmul(W1,A),np.transpose(x1))
        T = (T*T) / np.matmul(np.matmul(W1,A),np.transpose(W1));
        
        # If SHet for the given threshold is larger than all previous SHet values, then keep this value
        if TTT < T[0,0]:
            TTT = T[0,0]	
            sel_group = index
        
        iter_time = time() - ptm
        print("Calculated in ",iter_time)
    
    print('Finished in {} sec.'.format(time()-all_time))
    
    return(TTT,sel_group)


if __name__ == "__main__":
    #%%
    ## Input
    # GCTA summary stats
    res = pd.read_csv("stats_SNP=1:727.csv")
    # Tierpsy features per strain
    features = pd.read_csv("10_tierpsy_features.csv")
    
    #%%
    ## Preprocess stats and features/traits
    res = res.sort_values(by='feature')
    analysed_feat = res["feature"].tolist()
    res = res.set_index('feature')
    res = res[['b','se','p']]
    
    features = features[res.index]
    
    ## Get feature/trait correlation matrix (equivalent to the correlation between Wald statistics of each trait, according to equation 5 in Zhou et al 2015)
    CorrMatrix = features.corr()
    
    # Make sure that the traits are in the same order in res and CorrMatrix
    res = res.loc[CorrMatrix.index,:]
    res = res.values
    CorrMatrix = CorrMatrix.values
    
    ## Calculate Wald statistic per feature/trait
    X = res[:,0]/res[:,1]
    
    # Sample size
    n_sample = features.shape[0]
    # Number of traits
    nvar = res.shape[0]
    
    # Make SampleSize vector
    SampleSize = np.repeat([n_sample],nvar)
    
    ## Calculate SHet statistic
    SHet,sel_group = Trucated_TestScore(X, SampleSize, CorrMatrix, correct = True, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = False)