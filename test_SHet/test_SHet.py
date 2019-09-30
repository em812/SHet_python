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
import pandas as pd 
from SHet import Trucated_TestScore,EstimateGamma
from scipy.stats import gamma

if __name__ == "__main__":
    #%%
    ## Input
    # GCTA summary stats
    res = pd.read_csv("summary_stats_SNP=1:717.csv")
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
    SHet,sel_group = Trucated_TestScore(X, SampleSize, CorrMatrix, correct = True, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = True, return_group=True)
    
    ## Estimate gamma distribution
    N_permut = 100000
    para = EstimateGamma(N_permut, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = True, plot_gamma=True)
    
    ## Calculate pvalue
    pval = 1-gamma.cdf(SHet, a=para[0], loc=para[1], scale=para[2])