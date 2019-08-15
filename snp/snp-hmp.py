# -*- coding: utf-8 -*-

from pathlib import Path
from scipy.stats import hmean
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

file = Path('C:/Users/Jacob/Downloads/FEP/ecol241sir_fep_LMM_out_biallelic_and_multiallelic.txt')

df = pd.read_csv(file, sep='\t')

print(df.head())

window_width = 2000
stagger = int(window_width/2)
num_tests = len(df)
end_position = int(max(df['Position']))

# Value from Danny
alpha = 0.02787289

# Adjust p-vals so they can be plotted alongside HMPs
df['adjusted_pvals'] = df['pvals']*num_tests
df['colour'] = df['adjusted_pvals'].apply(lambda row: 'blue' if -np.log10(row) < -np.log10(alpha) else 'red')
plt.scatter(df['Position'], -np.log10(df['adjusted_pvals']), c=df['colour'])

for start in range(0, end_position - stagger, stagger):
    end = start + window_width
    window = df[df['Position'].between(start,end)]['pvals']
    hmp = hmean(window)
    num_tests_window = len(window)
    try:
        adjusted_hmp = hmp*(num_tests/num_tests_window)
    except ZeroDivisionError:
        adjusted_hmp = np.nan

    x = [start, end]
    y = [-np.log10(adjusted_hmp), -np.log10(adjusted_hmp)]
    
    if -np.log10(adjusted_hmp) > -np.log10(alpha):
        c = 'y'
    else:
        c = 'cyan'
    
    plt.plot(x,y, c=c)
    
plt.plot([0, end_position], [-np.log10(alpha), -np.log10(alpha)])
plt.title(f'Harmonic mean $p$-value of E.coli ({window_width/1000}kbp windows)')
plt.xlabel('Genome position')
plt.ylabel('$-\log_{10}($Adjusted $p$-value$)$')

fig = plt.gcf()
fig.set_size_inches(40, 20)
fig.savefig(Path(f'C:/Users/Jacob/Downloads/FEP/image_{window_width}.png'), dpi=100)
