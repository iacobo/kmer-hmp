# -*- coding: utf-8 -*-
"""
Script designed to:
    
    1. Reorder contigs for multiple draft files against a single reference file 
       using Mauve on Windows.
    2. Align reordered sequences with reference file using Mauve.
    3. Return MSA as list.
    
See Mauve documentation for further info:
    
    http://darlinglab.org/mauve/user-guide/progressivemauve.html
"""

import os
import glob
import subprocess
import datetime

from Bio import AlignIO
from Bio import Entrez

def order_contigs(reference, draft, output_dir, mauve_dir, java_dir='java'):
    """Call Mauve to reorder the contigs of a draft sequence relative to
    a reference sequence.
    Save results to `output_dir`.
    """
    command = f'{java_dir} -Xmx500m -cp {mauve_dir} org.gel.mauve.contigs.ContigOrderer -output {output_dir} -ref {reference} -draft {draft}'
    subprocess.call(command)
    
def align_seqs(reference, sequences, output='alignment.xmfa'):
    """Align [multiple] sequences to a reference sequence using Mauve.
    """
    command = f'progressiveMauve --output={output} {reference} '
    command += ' '.join(sequences)
    subprocess.call(command)
    
def grab_latest_alignment(suffixes, results_dir):
    """ Returns a list of file locations of the reordered FASTA files for
    multiple alignment.
    Often multiple successive alignments are created by Mauve, this returns
    the most recent.
    
    i.e.
    C:\\Users\\User\\Documents\\alignment_sequence_X\\alignment1
    C:\\Users\\User\\Documents\\alignment_sequence_X\\alignment2
    C:\\Users\\User\\Documents\\alignment_sequence_X\\alignment3   <<< 
    """
    seqs = []
    for suffix in suffixes:
        # WARNING: Max logic will fail if more than 9 files (i.e. '10' < '2', lexicographically)
        try:
            file = max(glob.glob(f'{results_dir}alignment_{suffix}\\alignment*\\{suffix}*.fas'))
            seqs.append(file)
        except FileNotFoundError as fnf_error:
            print(f'Error: {fnf_error}')
        except ValueError as val_error:
            print(f'Error: {val_error}')
    return seqs

def grab_fasta_files(location):
    """Grab locations of all FASTA files in given directory.
    """
    extensions = ['fasta','fa','faa','fna','ffn','frn']
    files = []
    for extension in extensions:
        files.extend(glob.glob(f'{location}*.{extension}'))
    return files

def download_genome(search_term, directory, filetype='fasta', email=None):
    """Download genome from NCBI to given directory.
    Must provide email address to Entrez.
    Return location of downloaded file.
    """
    if not email:
        email = input('Email address to query Entrez with: ')

    Entrez.email = email
    
    handle = Entrez.esearch(db='nucleotide', term=search_term)
    genome_id = Entrez.read(handle)['IdList'][0]
    record = Entrez.efetch(db='nucleotide', id=genome_id, rettype=filetype, retmode='text')
    
    filename = f'{directory}Record_{genome_id}.{filetype}'
    
    print(f'Writing: {filename}')
    with open(filename, 'w') as f:
        f.write(record.read())
    handle.close()
    
    return filename

def sort_multi_alignment_by_reference(multialignment, reference='reference'):
    """Given a multialignment file (multiple alignments), sorts alignments
    
    1) primarily on whether alignment contains a record from the reference genome
    2) secondarily on id of first element of each alignment.
    """
    
    ## TODO: INCLUDE LOGIC TO SORT ON ID OF REFERENCE RECORD
    ## DON'T ASSUME REFERENCE IS FIRST RECORD
    multialignment.sort(key = lambda x: (all([reference not in record.id for record in x]), int(x[0].id.split('-')[-1])))


####################################
## Main program
####################################
    
def main(base_path=None, reference=None, reorder=False, mauve_dir=None):
    
    ############################################
    ## Sequence files
    ############################################
    
    # Download reference genome
    searches = ['staphylococcus[orgn]','MSSA476 complete genome[title]','NC_002953.3[accession]']
    search_term = ' AND '.join(searches)
    output_dir = f'{base_path}reference_genome'
    
    # Download reference file if not specified
    if not reference:
        reference = download_genome(search_term, output_dir)

    # Locations of executables to run in cmd
    # Note: double quotes necessary for whitespace in directory names
    java_dir = '"C:\\Program Files (x86)\\Common Files\\Oracle\\Java\\javapath_target_10534109\\java.exe"'
    mauve_dir = '"C:\\Program Files (x86)\\Mauve 20150226\\Mauve.jar"'
    
    # Seq files of which to reorder contigs
    drafts_dir = f'{base_path}fa_files\\'
    results_dir = f'{base_path}ordered_contigs\\'
    alignment_filename = f'{base_path}alignment\\alignment.xmfa'

    if reorder:
        if mauve_dir:
            # Change directory to execute progressiveMauve (necessary if not in PATH)
            os.chdir(mauve_dir)
        # Container for draft file names
        draft_filenames = []
    
        print('Ordering sequence contigs relative to reference sequence...')
        t0 = datetime.datetime.now() 
        # Loops through files in directory of '.fa' format
        for draft in grab_fasta_files(drafts_dir):
            draft_filename = draft.split('\\')[-1].split('.')[-2]
            draft_filenames.append(draft_filename)
            output_dir = f'{results_dir}alignment_{draft_filename}'
            
            # Orders contigs of each draft sequence relative to reference
            order_contigs(reference, draft, output_dir, mauve_dir, java_dir)
            
        t1 = datetime.datetime.now()
        print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')
        
        print('Aligning sequences...')
        t0 = datetime.datetime.now() 
        # Aligns all draft sequences relative to reference
        seqs = grab_latest_alignment(draft_filenames, results_dir)
        align_seqs(reference, seqs, alignment_filename)
        t1 = datetime.datetime.now()
        print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')
    
    # Parse alignments with Biopython
    print('Loading alignment file...')
    t0 = datetime.datetime.now() 
    alignments = list(AlignIO.parse(alignment_filename, 'mauve'))
    # Sort alignments to match ordering of reference genome
    sort_multi_alignment_by_reference(alignments, reference='reference')
    t1 = datetime.datetime.now()
    print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')
    
    return alignments

if __name__ == '__main__':
    mauve_dir = 'C:\\Program Files (x86)\\Mauve 20150226'
    base_path = 'C:\\Users\\Jacob\\Downloads\\fusidic_data\\'
    reference = f'{base_path}reference_genome\\MSSA476.fasta'
    
    os.chdir(base_path)
    main(base_path=base_path, reference=reference, mauve_dir=mauve_dir)