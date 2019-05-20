# -*- coding: utf-8 -*-
"""
Script designed to:
    
    1. reorder contigs for multiple draft files against a single 
       reference file using Mauve on Windows.
    2. align reordered sequences with reference file using Mauve.
    3. calculate the Harmonic Mean p-value for each sliding window
       across the length of the sequences
       
    
Instructions from Mauve website:

        Reordering contigs from the command-line (batch mode)
        ----------------------------------------------------
        In situations where it is necessary to order contigs in a large number 
        of draft genomes it is often more desirable to automate the process 
        using command-line interfaces and scripts. Mauve Contig Mover supports 
        command-line operation through the Mauve Java JAR file.
        
        Given a reference genome file called “reference.gbk” and a draft genome 
        called “draft.fasta”, one would invoke the reorder program with the 
        following syntax:
        
        java -Xmx500m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output results_dir -ref reference.gbk -draft draft.fasta
        
        The file Mauve.jar is part of the Mauve distribution. On windows 
        systems it can usually be found in C:\Program Files\Mauve X\Mauve.jar 
        where X is the version of Mauve. On Mac OS X it is located inside the 
        Mauve application. For example, if Mauve has been placed in the OS X 
        applications folder, Mauve.jar can be found at 
        /Applications/Mauve.app/Contents/Resources/Java/Mauve.jar. 
        On Linux, Mauve.jar is simply at the top level of the tar.gz archive. 
        In the above example command, it will be necessary to specify the full 
        path to the Mauve.jar file.

Source:
    
    http://darlinglab.org/mauve/user-guide/reordering.html


Overview of objects:
    
              [                      ]  <-- `alignments`
               [XXX,  [XXX,  [XXXXXX,   <-- `alignment`
                XXX,   XXX,   XXXXXX]
                XXX],  XXX,              
                       XXX],
               (        )                  
                   (        )
                       (        )
                           (        )   <-- `window`
                ---    ---    ---       <-- `kmer`
                ---    ---     ---
                       ---      ---
                                 ---
                 .      .      ....     <-- p-values
                 .      .
                        .
                                  
Explanation of logic:
    
    Each k-mer in the sequences has an associated p-value from the GWAS
    performed on the sequences testing for correlation with phenotype.
    
    This program creates a window of defined length, and travelling laterally
    across the aligned sequences, takes the Harmonic Mean of the k-mers in
    said window.
    
    I.e. in a window of length 10, k-mers of k=5, and 4 aligned sequences
    (excluding the reference sequence), assuming no gap characters and assuming
    that the alignment is longer than the window, this contain generate:
        
        (10-5+1)*4 = 6*4 = 24 total k-mers
    
    and thus 24 values of which the Harmonic Mean will be taken.
    
Gap characters:
    
    If there are gap characters in a sequence, e.g.:
        
        ACACAC-GTG
        
    Within a given window, the program will ignore the gap characters, thus the
    k-mers returned for this window (of size 10) for this sequence would be:
        
        ACACA
         CACAC
          ACACG
           CACGT
            ACGTG

Window overlaps disjoint alignments:
    
    If there are multiple alignments ordered sequentially, some sliding windows
    will overlap with the discontinuity between alignments.
    
         [...  ...]
    [ACACGTGT][CGAT]
    
    In these instances k-mers are taken from each alignment intersection only.
    
    I.e. for k=2, window_size=6 in the above example:
        
        TG
         GT
             CG
              GA

Neighbouring alignments with different depths:
    
    Some neighbouring alignments may have different numbers of aligned
    sequences in them. E.g. 3 and 2:
    
     [...  ...]
    [ACTG][GTGTA]
    [CCTG][GTGGA]
    [TCTG]
    
    In these instances, if a window overlaps an alignment boundary, k-mers are
    taken from all aligned sequences. E.g. k=3, window_size=6:
        
     CTG  GTG
     CTG  GTG
     CTG
    
"""

import os
import glob
import subprocess
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import hmean
from bisect import bisect
import random
import matplotlib.pyplot as plt


#####################################
## MAUVE alignment reading functions
#####################################

def order_contigs(reference, draft, output_dir, mauve_dir, java_dir='java'):
    """Call Mauve to reorder the contigs of a draft sequence relative to
    a reference sequence.
    Save results to `output_dir`.
    """
    command = f'{java_dir} -Xmx500m -cp {mauve_dir} org.gel.mauve.contigs.ContigOrderer -output {output_dir} -ref {reference} -draft {draft}'
    subprocess.call(command)
    
def align_seqs(reference, sequences, output='alignment.xmfa'):
    """Align [multiple] sequences to a reference sequence using Mauve.
    Examples given here:
        
    http://darlinglab.org/mauve/user-guide/progressivemauve.html
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
    alignment1
    alignment2
    alignment3 <<<
    """
    seqs = []
    for suffix in suffixes:
        try:
            file = max(glob.glob(f'{results_dir}alignment_{suffix}\\alignment*\\{suffix}*.fas'))
            seqs.append(file)
        except FileNotFoundError as fnf_error:
            print(f'Error: {fnf_error}')
        except ValueError as val_error:
            print(f'Error: {val_error}')
    return seqs

def grab_fasta_files(location):
    """Grab all FASTA files in given location.
    """
    extensions = ['fasta','fa','faa','fna','ffn','frn']
    files = []
    for extension in extensions:
        files.extend(glob.glob(f'{location}*.{extension}'))
    return files

#####################################
## Sequence manipulation functions
#####################################
    
def get_kmers(k, sequence, ignore_gaps=True, upper=True):
    """Return a list of kmers of length `k` from the provided sequence.
    Kmers returned as upper case strings.
    If sequence length < k, return empty list.
    """
    if ignore_gaps:
        sequence = sequence.ungap(gap='-')
    
    length = len(sequence)
    if length >= k:
        if upper:
            return [sequence[i:i+k].upper() for i in range(length - k + 1)]
        else:
            return [sequence[i:i+k] for i in range(length - k + 1)]
    else:
        return []

def align_concat_len(x):
    """Return length of flattened list considering only first sub-elements."""
    return sum(len(sublist[0]) for sublist in x)

def find_2d_idx(cumsum, idx):
    """Convert a flat index into a nested 2d index for a nested list."""
    i1 = bisect(cumsum, idx)
    i2 = (idx - cumsum[i1 - 1]) if i1 > 0 else idx
    return (i1, i2)

def len_cumsum(x):
    """Return the cumsum of the lengths of the first element of each sublist."""
    c = [len(lst[0]) for lst in x]
    for i in range(1, len(c)):
        c[i] += c[i - 1]
    return c

def random_seq(n,alpha='ACGT'):
    """Return random Seq of length n."""
    return ''.join(random.choice(alpha) for _ in range(n))

def remove_reference(alignments, reference):
    """Removes the reference alignment from a multi-alignment."""
    new_alignments = []
    for alignment in alignments:
        # Include only non-reference records
        new_alignment = [record for record in alignment if reference not in record.id]
        # If reference record only record in alignment, include gap sequence
        if len(new_alignment) == 0:
            new_alignment = [SeqRecord(Seq('-'*len(alignment[0])))]
        new_alignments.append(new_alignment)
    return new_alignments

def alignment_window(window_size, stagger, alignments):
    """Given a a multi-alignment file, returns all sub-sequences of each
    record in each alignment for each 'sliding window' across the alignments.
    """
    
    # Grab length of concatenated alignments
    l = align_concat_len(alignments)
    
    # Grab cumulative sum of lengths of each alignment in alignments
    c = len_cumsum(alignments)
    
    # Generate start and end indices for each sliding window
    start_indices = np.array(range(0,l,stagger))
    end_indices = start_indices + window_size
    
    # Generate 2d indices for each window across the alignments
    starts = [find_2d_idx(c, i) for i in start_indices]
    ends = [find_2d_idx(c, i) for i in end_indices]
    
    # Grab all sub-sequence which intersect with each window.
    windows = []
    # Loop through each window's start, end
    for s, e in zip(starts,ends):
        window = []
        start_i, start_j = s
        end_i, end_j = e
        # Loop through each alignment which intersects window
        for i, alignment in enumerate(alignments[start_i:end_i+1]):
            if i == 0:
                start = start_j
            else:
                start = 0
            if i == end_i - start_i:
                end = end_j
            else:
                end = len(alignment[0])
            for record in alignment:
                window.append(record[start:end])
        windows.append(window)
    return windows


####################################
## Main program
####################################

if __name__ == '__main__':
    
    # If debug, code will skip reordering and alignmnet steps
    debug = 1
    
    # kmer k
    k = 31
    
    # Sliding window size for analysing kmers
    window_size = 10000
    stagger = window_size #int(window_size/2)
    
    # Default value for kmers without associated p-val (if nan, values excluded)
    default = np.nan
    
    ############################################
    ## kmer files
    ############################################
    
    # File locations
    base_path = 'C:\\Users\\Jacob\\Downloads\\fusidic_data\\'
    kmers_file = f'{base_path}fuc_400000_kmers.txt'
    pvals_file = f'{base_path}fuc_LMM_results_400000kmers_out.txt'
    total_kmers_file = f'{base_path}cipro_all_kmer_out.kmer.txt'
    
    if not debug:
    
        # Load list of kmers as pd series
        kmers_series = pd.read_csv(kmers_file, header=None, squeeze=True, names=['kmer'])
        
        # Load p-value info as pd dataframe
        dfpvals = pd.read_csv(pvals_file, sep='\t')
        
        # Load total list of kmers as pd dataframe
        total_kmers = pd.read_csv(total_kmers_file, header=None, names=['kmer'])
        # Add default p-value of 1
        total_kmers['p_score'] = 1    
        
        # Checks on file integrity
        assert kmers_series.is_unique, "List of kmers is not unique!"
        assert kmers_series.shape[0] == dfpvals.shape[0], "K-mer file lengths do not match!"
        assert sum(kmers_series.isin(total_kmers['kmer'])) == kmers_series.shape[0], "Not all kmers extract in full kmers list!"
        
        # Grab only necessary columns
        dfpvals['kmer'] = kmers_series
        dfpvals = dfpvals[['kmer','p_score']]
        
        # Set indices for update and dictionary conversion
        dfpvals = dfpvals.set_index('kmer')
        total_kmers = total_kmers.set_index('kmer')
        
        # Update total kmer list with pvals from 400k list
        total_kmers.update(dfpvals)
        
        # Create dictionary of kmers to p-values
        kmer_pvalues = total_kmers['p_score'].to_dict()
    
    
    
    ############################################
    ## Sequence files
    ############################################

    # Locations of executables to run in cmd
    # Note: double quotes necessary for whitespace in directory names
    java_dir = '"C:\\Program Files (x86)\\Common Files\\Oracle\\Java\\javapath_target_10534109\\java.exe"'
    mauve_dir = '"C:\\Program Files (x86)\\Mauve 20150226\\Mauve.jar"'
    
    # Seq files to reorder contigs of
    reference = f'{base_path}reference_genome\\MSSA476.fasta'
    drafts_dir = f'{base_path}fa_files\\'
    results_dir = f'{base_path}ordered_contigs\\'
    alignment_dir = f'{base_path}alignment\\alignment.xmfa'
    
    # Container for draft file names
    draft_filenames = []
    
    # Change directory to execute progressiveMauve (necessary if not in PATH)
    os.chdir('C:\\Program Files (x86)\\Mauve 20150226')
    
    if not debug:
        # Loops through files in directory of '.fa' format
        for draft in grab_fasta_files(drafts_dir):
            draft_filename = draft.split('\\')[-1].rstrip('.fa')
            draft_filenames.append(draft_filename)
            output_dir = f'{results_dir}alignment_{draft_filename}'
            
            # Orders contigs of each draft sequence relative to reference
            order_contigs(reference, draft, output_dir, mauve_dir, java_dir)
    
        # Aligns all draft sequences relative to reference
        seqs = grab_latest_alignment(draft_filenames, results_dir)
        align_seqs(reference, seqs, alignment_dir)
        
    # Parse alignments with Biopython
    alignments = list(AlignIO.parse(alignment_dir, "mauve"))
    alignments = remove_reference(alignments, reference)
    
    # Deal with smaller amount for testing
    if debug:
        alignments = alignments[:5]
    
    # Container for harmonic mean p-values
    hmps = []
    total_kmers_count = 0
    unknown_kmers = 0
    
    ###############################################################
    # Create windows
    windows = alignment_window(window_size, stagger, alignments)
    
    # Calculate HMP for each window
    for window in windows:
        pvalues = []
        for record in window:
            for kmer in get_kmers(k,record.seq):
                total_kmers_count += 1
                pval = kmer_pvalues.get(str(kmer),default)
                # Check reverse complement
                if pval is np.nan:
                    pval = kmer_pvalues.get(str(kmer.reverse_complement()),default)
                if pval is np.nan:
                    unknown_kmers += 1
                else:
                    pvalues.append(pval)
        
        # Calculate harmonic mean
        hmp_window = hmean(pvalues)
        
        # Add hmp for window to sequence container
        hmps.append(hmp_window)
    hmps = np.array(hmps)
    
    # Stats
    percent_unknown_kmers = 100*unknown_kmers/total_kmers_count
    print(f'Proportion of unknown kmers: {percent_unknown_kmers:.2f}%')
    print(f'Number of unknown kmers: {unknown_kmers}')
    print(f'Total number of kmers: {total_kmers_count}')
    
    # Displaying graphs
    alpha = 5*10**(-8)
    plt.scatter(np.arange(len(hmps)), -np.log(hmps))
    plt.plot([0,len(hmps)], [-np.log(alpha), -np.log(alpha)], linestyle='--', color='k')
    plt.xlabel("Window position")
    plt.ylabel("Harmonic Mean p-Value")
    
    plt.close()
    
    for i, hmp in enumerate(hmps[:400]):
        plt.plot([i*window_size, (i+1)*window_size], [-np.log(hmp), -np.log(hmp)], color=window_size)
    
    ###############################################################
#        
