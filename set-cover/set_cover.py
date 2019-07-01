# -*- coding: utf-8 -*-
"""
Determine smallest subset of genomes required to encapsulate > 80% of distinct
k-mers present in complete genome set.

Implementation: Greedy algorithm. 
                Not guaranteed minimum, but at most factor of log(n) off.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def sample_file(x,y):
    """Return DataFrame of size (1,y) with each element a random binary string
    of length x.
    """
    df = pd.DataFrame([''.join(map(str, np.random.randint(2, size=(x,)))) for _ in range(y)], columns=['bin'])
    return df

def len_first_line(file):
    """Return length of first line of file.
    """
    with open(file) as f:
        length = len(f.readline().strip())
    return length

def plot_curve(kmer_cumcounts, n_kmers, n_genomes, proportion):
    """Plot curve of number of genomes vs k-mer coverage.
    """
    plt.plot(range(len(kmer_cumcounts)), np.array(kmer_cumcounts)/n_kmers) # Remove hardcoding of length
    plt.plot([0,len(kmer_cumcounts)], [proportion, proportion], 'r--')
    plt.xlabel(f'Number of genomes ({n_genomes} total)')
    plt.ylabel(f'Proportion of k-mers covered \n({n_kmers} total)')
    plt.title('Set covering')
    plt.suptitle('Minimum number of genomes required to capture y% of all k-mers')

def main(patterns, pattern_index=None, weighted=True, proportion=0.8):
    # Load Patterns file as DataFrame
    # Row = k-mer (pattern of presence/absence)
    # Col = genome
    n_genomes = len_first_line(patterns)
    df = pd.read_fwf(patterns, header=None, widths = [1]*n_genomes)
    print('DataFrame loaded.')
        
    # Weighted version: pat weighted by number of k-mers sharing that pattern
    if weighted:
        pattern_index = pd.read_csv(pattern_index, names=['i'])
        weights = pattern_index['i'].value_counts().sort_index()
        
        assert len(df) == len(weights), "Index file does not match patterns file!"
        # Times row i by weight i  (1 > i, 0 > 0)
        df = df.multiply(weights, axis=0)
        print("Weights applied.")
    
    # Transpose dataframe to get:
    # Row = genome
    # Col = k-mer (pattern of presence/absence)
    df = df.T
    print('DataFrame transposed.')
    
    if weighted:
        n_kmers = len(pattern_index)
    else:
        # Actually total number of patterns, not k-mers
        n_kmers = len(df.columns)
    
    # Minimum number to catch
    threshold = proportion*n_kmers
    kmer_cumcounts = [0]
    indices = []
    
    # Set cover - greedy algorithm
    while kmer_cumcounts[-1] < threshold:
        sums =  df.sum(axis=1, skipna=True)
        try:
            # Only grab one (first) row (multiple lines may have same sum)
            row = df.loc[[sums.idxmax()]].iloc[[0]]
            index = row.index[0]
            indices.append(index)
            
            # Exclude non-zero columns from current row
            column_mask = np.logical_not(row).values[0]
            # Exclude current row
            row_mask = df.index != index
            # Update DataFrame for next iteration:
            df = df.loc[row_mask,column_mask]
            
            kmer_cumcounts.append(kmer_cumcounts[-1] + sums.max())
            print(f'Number of k-mers covered: {kmer_cumcounts[-1]}')
        except ValueError:
            print("No further novel k-mers present in this iteration.")
            break
    
    print(f'Covering set row numbers: {indices}')
    
    # Plot curve of number of genomes vs k-mer coverage
    plot_curve(kmer_cumcounts, n_kmers, n_genomes, proportion)
    
    return indices

if __name__ == '__main__':
    
    base_path = Path('C:/Users/Jacob/Downloads/fusidic_data/static_files')
    patterns = base_path / 'fusidic_acid_patternKey.txt' 
    pattern_index = base_path / 'fusidic_acid_patternIndex.txt'
    
    main(patterns, pattern_index)
