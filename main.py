# -*- coding: utf-8 -*-
"""
Script designed to:
    
    1. reorder contigs for multiple draft files against a single 
       reference file using Mauve on Windows.
    2. align reordered sequences with reference file using Mauve.
    3. calculate the Harmonic Mean p-value for each sliding window
       across the length of the sequences
    4. Plot the HMP's 'Manatan-plot' style (optionall also plot kmers)
       
    
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
                .      .      ....      <-- p-values
                .      .
                       .
                                  
Explanation of logic:
    
    Each k-mer in the sequences has an associated p-value from the GWAS
    performed on the sequences (testing for correlation with phenotype).
    
    This program creates a window of defined length, and travelling laterally
    across the aligned sequences, takes the Harmonic Mean of the p-values of 
    the k-mers in said window.
    
    i.e. in a window of length 10, k-mers of k=5, and 4 aligned sequences
    (excluding the reference sequence), assuming no gap characters and assuming
    that the alignment is longer than the window, this will contain:
        
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
    
    i.e. for k=2, window_size=6 in the above example:
        
        TG
         GT
             CG
              GA

Neighbouring alignments with different depths:
    
    Some neighbouring alignments may have different numbers of aligned
    sequences in them. e.g. 3 and 2:
    
     [...  ...]
    [ACTG][GTGTA]
    [CCTG][GTGGA]
    [TCTG]
    
    In these instances, if a window overlaps an alignment boundary, k-mers are
    taken from all aligned sequences. e.g. k=3, window_size=6:
        
     CTG  GTG
     CTG  GTG
     CTG
    
    k-mes 'across' neighbouring alignments are not taken, as these may not
    correspond to real kmers present in the genome (since the ordering of the
    alignment may not be reflected in reality).
    
"""

import os
import glob
import subprocess
import random
import datetime

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import hmean


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
    C:\\Users\\User\\Documents\\alignment_sequence_X\\alignment1
    C:\\Users\\User\\Documents\\alignment_sequence_X\\alignment2
    C:\\Users\\User\\Documents\\alignment_sequence_X\\alignment3   <<< 
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


#############################################################
## ALTERNATIVE MAIN PROGRAM LOGIC TESTING
#############################################################

## Alt logic designed to creat DF of p-values per position for reuse of values 
# for different window plots, and plotting individual kmers Manhattan-plot style

def alignment2df(alignments,k,dictionary,default=np.nan):
    """Create DataFrame from multi-alignment. Dissects alignment to produce
    k-mers and given a dictionary of k-mer:p-value creates a row in the dataframe
    for each position in each sequence with an associated p-value (with respect 
    to the k-mer whose leftmost basepair starts at that position).
    
    Gap characters are ignored for generating k-mers.
    
    No p-values are associated to positions which are a gap character.
    """
    columns = ['kmer','p_val','position','absolute_pos','alignment','sequence']
    
    temp = []
    
    for i, alignment in enumerate(alignments):
        surplus = sum(len(alignment[0]) for alignment in alignments[:i])
        for record in alignment:
            seq_id = record.id.split("/")[0].split("\\")[-1].split(".")[0]
            seq_id = record_dict[seq_id]
            # Each position in sequence gets associated p-value
            for position, char in enumerate(record):
                # Position determined by leftmost value of k-mer
                # If no k-mer starts here (i.e. this pos is a gap char), default (NAN) given
                #if char = '-':
                #    temp.append(['', default, position, position+surplus, i, seq_id])
                if char != '-':
                    kmer = record[position:].seq.ungap(gap='-')[:k].upper()
                    pval = dictionary.get(kmer,default)
                    # Check reverse complement
                    if pval is default:
                        kmer = kmer.reverse_complement()
                        pval = dictionary.get(kmer,default)
                    if pval is not default: # not elif! checking condition twice
                        temp.append([str(kmer), pval, position, position+surplus, i, seq_id])
    
    # Convert list to dataframe
    df = pd.DataFrame(data=temp, columns=columns)
    df = df.set_index(['absolute_pos','alignment','sequence','position'])
    df = df.sort_index()
    
    return df

def get_hmps(df, window_size):
    """
    Get hmps for each sliding window.
    """
    hmps = []
    stagger = int(window_size/2)
    
    # Grab length of concatenated alignments (largest absolute position, + 1 for 0 indexing)
    #l = align_concat_len(alignments)
    l = max(df.index)[0]+1
    
    # Generate start and end indices for each sliding window
    start_indices = np.array(range(0,l,stagger))
    end_indices = start_indices + window_size
    
    # Takes slice of absolute positions (index 0)
    idx = pd.IndexSlice
    hmps = np.array([hmean(df.loc[idx[start:end,:,:,:],:]['p_val'].dropna()) for start, end in zip(start_indices,end_indices)])
    
    return hmps

#############################################################
## Visualisation functions
#############################################################

def remove_duplicate_legends():
    """Remove duplicate labels on a matplotlib legend.
    """
    handles, labels = plt.gca().get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)
    
    plt.legend(newHandles, newLabels)
    

def plot_hmps(df, hmps, window_size):
    """
    Plot p-values of sliding windows vs window position across sequence.
    Optionall also plot 'Manhattan Plot' of individual kmer p-values vs
    kmer position in multi-sequence.
    """
    stagger = int(window_size/2)
        
    # Sliding windows plot
    for i, hmp in enumerate(hmps):
        plt.plot([i*stagger, i*stagger+window_size], [-np.log(hmp), -np.log(hmp)], color=colors[window_size], label=f'{window_size}bp')
    
    remove_duplicate_legends()
    plt.xlabel("Genome position")
    plt.ylabel("Adjusted Harmonic Mean p-value (-log)")


####################################
## Main program
####################################

if __name__ == '__main__':
    
    # If debug, code will skip reordering and alignmnet steps
    debug = 1
    
    # kmer k
    k = 31
    
    # Sliding window size for analysing kmers
    window_sizes = (100,1000,10000,100000)
    # INSERT LOGIC TO DETERMINE STAGGER SIZE
    colors = {10:'purple',100:'blue',1000:'red',10000:'yellow',100000:'green',1000000:'orange'}
    
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
    total_pvals_file = f'{base_path}saur_992_derval_fusidic_acid_all_kmers_LMM_pvals_only.txt.'
    
    if not debug:
        print("Loading k-mers and p-values...")
        t0 = datetime.datetime.now() 
        # Load p-value info as pd dataframe
        dfpvals = pd.read_csv(total_pvals_file, header=None, names=['p_score'])
        
        # Load total list of kmers as pd series
        total_kmers = pd.read_csv(total_kmers_file, header=None, squeeze=True, names=['kmer'])
        
        # Checks on file integrity
        assert total_kmers.is_unique, "List of kmers is not unique!"
        assert total_kmers.shape[0] == dfpvals.shape[0], "K-mer, p-value file lengths do not match!"
        
        # Combine add kmers series as column to p-value DataFrame
        dfpvals['kmer'] = total_kmers
        
        # Set index for dictionary conversion
        dfpvals = dfpvals.set_index('kmer')
        
        # Create dictionary of kmers to p-values
        kmer_pvalues = dfpvals['p_score'].to_dict()
        
        print("K-mers and p-values loaded.")
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    
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
        print("Ordering sequence contigs relative to reference sequence...")
        t0 = datetime.datetime.now() 
        # Loops through files in directory of '.fa' format
        for draft in grab_fasta_files(drafts_dir):
            draft_filename = draft.split('\\')[-1].rstrip('.fa')
            draft_filenames.append(draft_filename)
            output_dir = f'{results_dir}alignment_{draft_filename}'
            
            # Orders contigs of each draft sequence relative to reference
            order_contigs(reference, draft, output_dir, mauve_dir, java_dir)
        print("All sequences ordered.")
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
        
        print("Aligning sequences...")
        # Aligns all draft sequences relative to reference
        seqs = grab_latest_alignment(draft_filenames, results_dir)
        align_seqs(reference, seqs, alignment_dir)
        print("Sequences aligned.")
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
        
    # Parse alignments with Biopython
    print("Loading alignment file...")
    t0 = datetime.datetime.now() 
    alignments = list(AlignIO.parse(alignment_dir, "mauve"))
    alignments = remove_reference(alignments, reference)
    print("Alignment file loaded (and reference sequence removed).")
    t1 = datetime.datetime.now()
    print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    
    
    ############################################
    ## P-Value/Position DataFrame
    ############################################
    
    # Dictionary for sequence ids
    # NEED TO NOT USE GLOBALS
    record_ids = sorted(list(set(record.id.split("/")[0].split("\\")[-1].split(".")[0] for alignment in alignments for record in alignment)))
    record_dict = {record_id:i for i, record_id in enumerate(record_ids)}
    
    # Create DF
    print("Converting alignment file to DataFrame with p-value/position as row...")
    t0 = datetime.datetime.now() 
    df = alignment2df(alignments[:50],k,kmer_pvalues)
    print("Dataframe successfully created.")
    t1 = datetime.datetime.now()
    print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    
    
    ############################################
    ## HMPs, Windows, Visualisation
    ############################################
    
    total_sequence_length = l = max(df.index)[0]+1
    upper_exp = int(math.log10(total_sequence_length))
    
    window_sizes = (10**e for e in range(4,upper_exp+1))
    
    alpha = 5*10**(-8)
    
    for window_size in window_sizes:
        print(f"Generating visualisation for window size {window_size}bp...")
        t0 = datetime.datetime.now() 
        stagger = int(window_size/2)
        
        # Calculate HMP's for each window size
        hmps = get_hmps(df,window_size)
        plot_hmps(df,hmps,window_size)
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
        
    # Plot kmers Manhattan-plot style (scatter graph)
    print(f"Generating Manhattan-plot...")
    t0 = datetime.datetime.now() 
    plt.scatter(df.index.to_frame()['absolute_pos'], -np.log(df['p_val']))
    t1 = datetime.datetime.now()
    print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    plt.show()
    
    # Save image
    plt.savefig('C:\\Users\\Jacob\\Downloads\\fig.svg')
        
        ###############################################################
        # Stats
        ###############################################################
        
        # Stats
#        percent_unknown_kmers = 100*unknown_kmers/total_kmers_count
#        print(f'Proportion of unknown kmers: {percent_unknown_kmers:.2f}%')
#        print(f'Number of unknown kmers: {unknown_kmers}')
#        print(f'Total number of kmers: {total_kmers_count}')
