# -*- coding: utf-8 -*-
"""
Script to:
    
    1. Determine smallest subset of genomes required to encapsulate >= x% of all
       k-mers present in complete genome set.
       
Effectively a set covering problem.

Idea 1: Greedy algorithm. Not guaranteed minimum, but at most factor of log(n) off.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

base_path = 'C:\\Users\\Jacob\\Downloads\\fusidic_data\\'
patterns = f'{base_path}patternindex\\cipro_all_kmer_out.patternKey.txt' 
indices = f'{base_path}patternindex\\cipro_all_kmer_out.patternIndex.txt'

load_file = True

if load_file:
    dft = pd.read_csv(patterns, names=['bin']) # Dataframe where columns are genome presence, rows are k-mers
    pattern_count = len(dft)
    genome_count = len(dft['bin'].values[0])
    
    # Testing, only consider first n genomes (not enough RAM to deal with all ~1000 at once, ~400 is memory limit
    n = 390
    i = 4
    dft['bin'] = dft['bin'].str[(i-1)*n:i*n]
    
    ## TESTING
    #dft = pd.DataFrame([''.join(map(str, np.random.randint(2, size=(n,))))for x in range(100)], columns=['bin'])
    
    # Remove any trailing spaces
    dft['bin'] = dft['bin'].str.rstrip()
    #Convert single string column to multiple binary columns
    dft = dft['bin'].apply(lambda x: pd.Series(list(x)))
    # Convert strings to ints
    dft = dft.apply(pd.to_numeric)
    
    #Transpose dataframe to get columns as kmers and rows as genomes
    dft = dft.T

kmer_count = [0]
proportion = 0.9
threshold = proportion*len(dft.columns)
indices = []

print(dft)

while kmer_count[-1] < threshold:
    sums =  dft.sum(axis=1, skipna=True)
    try:
        row = dft.loc[[sums.idxmax()]].iloc[[0]] # ONLY GRAB ONE ROW (multiple lines may have same sum)
        index = row.index[0]
        indices.append(index)
        print('\n', row)
        columns = row.loc[:, (row == 0).any()].columns
        
        # Update DataFrame for next iteration
        dft = dft[columns] # Include only good columns
        dft = dft.loc[dft.index != index] # exclude selected row
        
        kmer_count.append(kmer_count[-1] + sums.max())
        print(f'Number of k-mers covered: {kmer_count[-1]}')
    except ValueError:
        print("No further novel k-mers present in this iteration.")
        break
    
print('-----')
print(f'Covering set row numbers: {indices}')

plt.plot(range(len(kmer_count)), np.array(kmer_count)/pattern_count) # Remove hardcoding of length
plt.plot([0,len(kmer_count)], [proportion, proportion], 'r--')
plt.xlabel(f'Number of genomes ({pattern_count} total)')
plt.ylabel('Proportion of k-mer patterns covered')
plt.title('Set covering: minimum number of genomes required to capture x% of all k-mers')
plt.text(0,-0.3,f'Note: only considering {i}th {n} genomes. Global performance may be better.')
