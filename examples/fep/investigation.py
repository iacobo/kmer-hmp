import pandas as pd
from pathlib import Path

"""
Script to investigate whether any genome files are missing etc.
Which of the genome files are positive for FEP resistance.
"""

directory = Path('FEP/genome_references')
comid2guid = directory / 'comids.tsv'
comids = directory / 'ecol241.comidPheno.sir_fep.txt'
output = directory / 'guids.txt'

comids_df = pd.read_csv(comids, sep=' ')
comid2guid_df = pd.read_csv(comid2guid, sep='\t', names=['Comid','Guid'])

## Subset of those positive for FEP resistance
fep_comids_df = comids_df[comids_df['sir_fep']==1]

# Only 208 lines in mapping file
intersection = pd.merge(comids_df, comid2guid_df, on='Comid', how='inner')
print(intersection.head())

mask = comids_df['Comid'].isin(comid2guid_df['Comid'])

# 33 lines appear to be missing from map
missing_mask = mask ^ 1

# Write GUIDS to file
with open(output, "w") as file:
    file.write('\n'.join(intersection['Guid']))
