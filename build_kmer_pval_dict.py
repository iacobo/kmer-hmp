# -*- coding: utf-8 -*-
"""
Script designed to:
    
    1. Load files containing k-mer and associated p-value information and return
       a dataframe combining this info.
"""

import os
import pandas as pd
    
def main(base_path=None):
    # File locations
    total_kmers_file = f'{base_path}cipro_all_kmer_out.kmer.txt'
    total_pvals_file = f'{base_path}saur_992_derval_fusidic_acid_all_kmers_LMM_pvals_only.txt'

    # Load p-value info as pd dataframe
    dfpvals = pd.read_csv(total_pvals_file, names=['p_score'])
    # Load total list of kmers as pd series
    total_kmers = pd.read_csv(total_kmers_file, names=['kmer'], squeeze=True)
    
    # Checks on file integrity
    assert total_kmers.is_unique, "Duplicate k-mer entries encountered in file!"
    assert total_kmers.shape[0] == dfpvals.shape[0], "K-mer, p-value file lengths do not match!"
    assert total_kmers.map(len).max() == total_kmers.map(len).min(), "Not all k-mers same length!"
    
    # Set kmers series as column in p-value DataFrame
    dfpvals['kmer'] = total_kmers
    
    # Set "key" column as index and convert to dict
    dfpvals = dfpvals.set_index('kmer')
    kmer_pvalues = dfpvals['p_score'].to_dict()
   
    return kmer_pvalues

if __name__ == '__main__':
    base_path = 'C:\\Users\\Jacob\\Downloads\\fusidic_data\\'
    os.chdir(base_path)
    main(base_path=base_path)