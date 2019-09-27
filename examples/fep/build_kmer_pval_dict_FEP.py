# -*- coding: utf-8 -*-
"""
Load files containing k-mer and associated p-value information and return
a dataframe combining this info.
"""

import pandas as pd
from pathlib import Path

def get_kmer_subset(total_kmers, indices, one_indexed=True):
    
    if one_indexed:
        indices['indices'] -= 1
    
    indices = indices['indices'].tolist()
    
    subset_kmers = total_kmers[total_kmers.index.isin(indices)]
    
    return subset_kmers
    
    
def main(kmers, kmers_used_index, pvals):
    # Load list of kmers as series
    kmers = pd.read_csv(kmers, names=['kmer'], squeeze=True)
    
    kmers_used_index = pd.read_csv(kmers_used_index, names=['indices'])
    kmers = get_kmer_subset(kmers, kmers_used_index)
    
    # Load p-values as dataframe
    dfpvals = pd.read_csv(pvals, sep='\t')
    
    # Checks on file integrity
    assert kmers.is_unique, "Duplicate k-mer entries encountered in file!"
    assert kmers.shape[0] == dfpvals.shape[0], "K-mer, p-value file lengths do not match!"
    assert kmers.map(len).max() == kmers.map(len).min(), "Not all k-mers same length!"
    
    # Set kmers series as column in p-value DataFrame
    dfpvals = pd.merge(dfpvals, kmers, left_on='LMM_kmers', right_index=True)
    
    # Set "key" column as index and convert to dict
    dfpvals = dfpvals.set_index('kmer')
    kmer_pvalues = dfpvals['p_lrt'].to_dict()
   
    return kmer_pvalues

if __name__ == '__main__':
    # File locations
    folder = Path('C:/Users/Jacob/Downloads/FEP/static_files') #Path('C:/Users/Jacob/Downloads/fusidic_data/static_files')
    kmers = folder / 'ecol241_sir_fep_kmers.txt.' #'fusidic_acid_kmers.txt'
    pvals = folder / 'ecol241_sir_fep_LMM_allkmers_out.txt' #'fusidic_acid_pvals.txt'
    
    kmers_used_index = folder / 'ecol241_sir_fep_kmers_used.txt'
    
    kmer_pvalues = main(kmers, kmers_used_index, pvals)
