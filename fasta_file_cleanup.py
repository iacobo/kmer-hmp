import textwrap
from pathlib import Path

"""
Code to convert fasta files with overly long lines into 50 char limited new-line
files that Mauve can deal with.
"""

directory = Path('FEP/genomes/draft_genomes')

for file in directory.glob('*.fa'):
        
    outputlines = []

    with open(file, mode='r') as csv_file:
        for line in csv_file:
            lines = textwrap.wrap(line, width=100)
            outputlines.extend(lines)
    
    outputlines = [line + '\n' for line in outputlines]
    
    file_name = file.name.split('.')[0] + '_EDIT.fa'
    output_file = directory / 'edited' / file_name
    
    with open(output_file, mode='w') as csv_file_output:
        csv_file_output.writelines(outputlines)
    
    print('Done')
