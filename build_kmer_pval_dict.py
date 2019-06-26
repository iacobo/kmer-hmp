# -*- coding: utf-8 -*-
"""
Load files containing k-mer and associated p-value information and return
a dataframe combining this info.
"""

import pandas as pd
from pathlib import Path
    
def main(kmers, pvals):
    # Load list of kmers as series
    kmers = pd.read_csv(kmers, names=['kmer'], squeeze=True)
    # Load p-values as dataframe
    dfpvals = pd.read_csv(pvals, names=['p_score'])
    
    # Checks on file integrity
    assert kmers.is_unique, "Duplicate k-mer entries encountered in file!"
    assert kmers.shape[0] == dfpvals.shape[0], "K-mer, p-value file lengths do not match!"
    assert kmers.map(len).max() == kmers.map(len).min(), "Not all k-mers same length!"
    
    # Set kmers series as column in p-value DataFrame
    dfpvals['kmer'] = kmers
    
    # Set "key" column as index and convert to dict
    dfpvals = dfpvals.set_index('kmer')
    kmer_pvalues = dfpvals['p_score'].to_dict()
   
    return kmer_pvalues

if __name__ == '__main__':
    # File locations
    folder = Path('C:/Users/Jacob/Downloads/fusidic_data')
    kmers = folder / 'fusidic_acid_kmers.txt'
    pvals = folder / 'fusidic_acid_pvals.txt'
    
    main(kmers, pvals)
