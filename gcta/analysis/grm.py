"""
GRM (Genetic Relationship Matrix) 계산 모듈 - 고급 멀티프로세싱 최적화
"""
import numpy as np
import pandas as pd
from scipy.linalg.blas import dsyrk, dgemm
import psutil
import os
from multiprocessing import Pool, cpu_count, shared_memory
from functools import partial
import pickle
import time

# standardization functions
def compute_whitening_coefficients(afreq):
    p1 = 2.0 * afreq * (1.0 - afreq)  # P(hetero)
    p2 = afreq * afreq                # P(homo alt)
    
    # Variance terms
    var1 = p1 * (1.0 - p1)    # Var(I1)
    var2 = p2 * (1.0 - p2)    # Var(I2) 
    cov12 = -p1 * p2          # Cov(I1, I2)
    
    # Cholesky decomposition elements
    l11 = np.sqrt(var1)
    l21 = cov12 / l11
    l22 = np.sqrt(var2 - l21 * l21)
    
    # Whitening matrix (inverse of Cholesky)
    w11 = 1.0 / l11
    w21 = -l21 / (l11 * l22)
    w22 = 1.0 / l22
    
    return p1, p2, w11, w21, w22

def std_additive(genotype, afreq):
    missing_mask = (genotype == -9)
    mean_geno = 2.0 * afreq
    std_geno = np.sqrt(2.0 * afreq * (1.0 - afreq))
    
    geno_std = (genotype - mean_geno) / std_geno
    geno_std[missing_mask] = 0.0
    return geno_std

def std_dominant(genotype, afreq):
    missing_mask = (genotype == -9)
    geno_gramschmidt = np.where(
        genotype == 1, 
        2.0 * afreq, 
        np.where(
            genotype == 2, 
            4.0 * afreq - 2.0, 
            0.0)
    ).astype(np.float64)
    mean_geno = 2.0 * afreq**2
    std_geno = 2.0 * afreq * (1.0 - afreq)
    
    geno_std = (geno_gramschmidt - mean_geno) / std_geno
    geno_std[missing_mask] = 0.0
    return geno_std

def std_indicator(genotype, afreq):
    missing_mask = (genotype == -9)
    p1, p2, w11, w21, w22 = compute_whitening_coefficients(afreq)
    
    hetero_indicator = (genotype == 1).astype(np.float64)
    homo_indicator = (genotype == 2).astype(np.float64)
    
    # centering
    d0 = hetero_indicator - p1
    d1 = homo_indicator - p2

    # whitening
    hetero_std = d0 * w11
    hetero_std[missing_mask] = 0.0        
    homo_std = d0 * w21 + d1 * w22
    homo_std[missing_mask] = 0.0

    return hetero_std, homo_std

# grm function
def calculate_grm_by_std(
    genotype,
    afreq,
    standardization,
    ):
    n_samples, n_snps = genotype.shape
    
    grm_numer = np.zeros((n_samples, n_samples), dtype=np.float64, order='C')
    grm_denom = np.zeros((n_samples, n_samples), dtype=np.int32, order='C')
    missing_mask = (genotype == -9)
    
    if standardization == "add":
        geno_std = std_additive(genotype, afreq)
        grm_numer = np.dot(geno_std, geno_std.T)
        grm_denom = (~missing_mask).astype(float) @ (~missing_mask).astype(float).T
        return [grm_numer/grm_denom]
    
    elif standardization == "dom":
        geno_std = std_dominant(genotype, afreq)
        grm_numer = np.dot(geno_std, geno_std.T)
        grm_denom = (~missing_mask).astype(float) @ (~missing_mask).astype(float).T
        return [grm_numer/grm_denom]
    
    elif standardization == "factor":
        # std genotype
        het_std, hom_std = std_indicator(genotype, afreq)
        
        # make GRMs
        grm_numer = np.dot(het_std, het_std.T)
        grm_denom = (~missing_mask).astype(float) @ (~missing_mask).astype(float).T
        grm_het = grm_numer/grm_denom
        
        grm_numer = np.dot(hom_std, hom_std.T)
        grm_denom = (~missing_mask).astype(float) @ (~missing_mask).astype(float).T
        grm_hom = grm_numer/grm_denom
        
        grm_numer = np.dot(het_std, hom_std.T)
        grm_numer = grm_numer + grm_numer.T
        grm_denom = (~missing_mask).astype(float) @ (~missing_mask).astype(float).T
        grm_denom = grm_denom + grm_denom.T
        grm_hethom = grm_numer / grm_denom
        
        return [grm_het, grm_hom, grm_hethom]
    
    # elif standardization == "hethom":
    #     het_std, hom_std = std_indicator(genotype, afreq)
    #     grm_numer = np.dot(het_std, hom_std.T)
    #     grm_numer = grm_numer + grm_numer.T
        
    #     grm_denom = (~missing_mask).astype(float) @ (~missing_mask).astype(float).T
    #     grm_denom = grm_denom + grm_denom.T
        
    #     return [grm_numer / grm_denom]
    else:
        raise ValueError(f"Invalid standardization: {standardization}")
    

if __name__ == "__main__":
    pass