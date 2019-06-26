# -*- coding: utf-8 -*-
"""
Script to:
    
    1. Determine smallest subset of genomes required to encapsulate >= x% of all
       k-mers present in complete genome set.
       
Effectively a set covering problem.

Implementation: Greedy algorithm. 
                Not guaranteed minimum, but at most factor of log(n) off.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# TODO: Weighted version. Change DF to have n instead of 1 for each pattern.

def sample_file(x,y):
    """Return DataFrame of size (1,y) with each element a random binary string
    of length x.
    """
    df = pd.DataFrame([''.join(map(str, np.random.randint(2, size=(x,)))) for _ in range(y)], columns=['bin'])
    return df

if __name__ == '__main__':
    
    base_path = 'C:\\Users\\Jacob\\Downloads\\fusidic_data\\'
    patterns = f'{base_path}patternindex\\cipro_all_kmer_out.patternKey.txt' 
    pattern_index = f'{base_path}patternindex\\cipro_all_kmer_out.patternIndex.txt'
    
    load_file = True
    whitespace = False
    weighted = True
    
    # Load Patterns file as DataFrame
    if load_file:
        dft = pd.read_csv(patterns, names=['bin']) # Dataframe where columns are genome presence, rows are k-mer patterns
        
        assert max(dft['bin'].apply(len)) == min(dft['bin'].apply(len)), "Not all patterns of same length!"
        
        pattern_count = len(dft)
        genome_count = len(dft['bin'].values[0])
        
        # Only consider n genomes (not enough RAM to deal with all ~1000 at once, ~400 is memory limit
        n = 450
        start = 542 - n
        stop = 542
        dft['bin'] = dft['bin'].str[start:stop]
        print('Dataframe slice taken.')
        
        if whitespace:
            # Remove any trailing spaces
            dft['bin'] = dft['bin'].str.rstrip()
        
        # Split string column into single-char columns
        dft = dft['bin'].apply(lambda x: pd.Series(list(x)))
        print('Binary string converted to columns.')
        
        # Convert '0', '1' to ints
        dft = dft.apply(pd.to_numeric)
        print('Column types converted to numeric.')
        
        # Weighted version
        if weighted:
            patinddf = pd.read_csv(pattern_index, names=['i'])
            weights = patinddf['i'].value_counts().sort_index()
            
            assert len(dft) == len(weights), "Index file does not match patterns file!"
            # Times row i by weight i  (1 > i, 0 > 0)
            dft = dft.multiply(weights, axis=0)
            print("Weights applied.")
        
        #Transpose dataframe to get columns as kmers and rows as genomes
        dft = dft.T
    
    # Initialise cumulative counts at 0
    kmer_cumcounts = [0]
    proportion = 0.9
    
    if weighted:
        total = sum(weights)
    else:
        total = pattern_count
    
    threshold = proportion*total
        
    indices = []
    
    print('DataFrame loaded.')
    
    # Set cover - greedy algorithm
    while kmer_cumcounts[-1] < threshold:
        sums =  dft.sum(axis=1, skipna=True)
        try:
            row = dft.loc[[sums.idxmax()]].iloc[[0]] # ONLY GRAB ONE ROW (multiple lines may have same sum)
            index = row.index[0]
            indices.append(index)
            print('\n', row)
            
            # Update DataFrame for next iteration
            column_mask = np.logical_not(row).values[0] # Exclude columns in present row
            row_mask = dft.index != index # Exclude present row
            dft = dft.loc[row_mask,column_mask]
            
            kmer_cumcounts.append(kmer_cumcounts[-1] + sums.max())
            print(f'Number of k-mers covered: {kmer_cumcounts[-1]}')
        except ValueError:
            print("No further novel k-mers present in this iteration.")
            break
    
    # Give absolute ot relative index values
    indices = np.array(indices) + start
    
    print('-----')
    print(f'Covering set row numbers: {indices}')
    
    # Plot coverage curve
    plt.plot(range(len(kmer_cumcounts)), np.array(kmer_cumcounts)/total) # Remove hardcoding of length
    plt.plot([0,len(kmer_cumcounts)], [proportion, proportion], 'r--')
    plt.xlabel(f'Number of genomes ({genome_count} total)')
    plt.ylabel(f'Proportion of k-mer (patterns) covered \n({total} total)')
    plt.title('Set covering')
    plt.suptitle('Minimum number of genomes required to capture y% of all k-mers')
    plt.text(0, -0.35, f'Note: only considering {n} genomes ({start} - {stop}). \nGlobal performance may be better.')
