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
    
              [                    ]  <-- `alignments`
               [XXX,  [XXX,  [XXXX,   <-- `alignment`
                XXX],  XXX,   XXXX]            
                       XXX],
               (        )                  
                   (        )
                       (        )
                           (        )   <-- `window`
                ---    ---    ---       <-- `kmer`
                ---    ---     ---
                       ---    ---
                               ---
                .      .      ..        <-- p-values
                .      .      ..
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

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import plotly.plotly as py
import plotly.graph_objs as go

from scipy.stats import hmean
from scipy.integrate import dblquad

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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

def subset_proportion(subset, completeset):
    """Given a subset, returns the relative size of the subset compared to
    completeset.
    Raise error if subset not subset of completeset.
    """
    subset = set(subset)
    completeset = set(completeset)
    
    len_intersection = len(subset.intersection(completeset))
    assert len_intersection == len(subset), "First element is not a subset of second element."
    
    return len_intersection/len(completeset)
    
    
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
    """Removes the reference alignment from a multi-alignment. If an alignment
    consists of only the reference sequence, it replaces it with a gapped
    sequence of the same length.
    """
    new_alignments = []
    for alignment in alignments:
        # Include only non-reference records
        new_alignment = [record for record in alignment if reference not in record.id]
        # If reference record only record in alignment, include gap sequence
        if len(new_alignment) == 0:
            new_alignment = [SeqRecord(Seq('-'*len(alignment[0])))]
        new_alignments.append(new_alignment)
    return new_alignments

def sort_multi_alignment_by_reference(multialignment, reference='reference'):
    """Given a multialignment file (multiple alignments), sorts
    
    1) primarily on whether alignment contains a record fromt he reference genome
    2) secondarily on id of first element of each alignment.
    """
    multialignment.sort(key = lambda x: (all([reference not in record.id for record in x]), int(x[0].id.split('-')[-1])))

#############################################################
## HMP FUNCS
#############################################################
    
def landau_integrand(t,x,mu,sigma):
    return (1/(np.pi*sigma))*np.exp(-t*(x-mu)/sigma - (2/np.pi)*t*np.log(t))*np.sin(2*t)

def combinedp(mu,sigma,hmp):
    return dblquad(landau_integrand, 0, np.inf, 1/hmp, np.inf, args=(mu,sigma))

############################################################
## TEXT FORMATTING ETC
############################################################

def sigfigs(n):
    """Return string formatted integer with K or M for thousands, millions
    sig figs.
    
    E.g. 1000 > 1K
    5000000 > 5M
    """
    if n > 1000000:
        return f'{n/1000000}M'
    elif n > 1000:
        return f'{n/1000}K'

#############################################################
## CONVERT ALIGNMENTS OBJECT TO DATAFRAME
#############################################################

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
                # If no k-mer starts here (i.e. this pos is a gap char), no value given
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

def alignment2array(alignments,k,dictionary,default=np.nan):
    """
    """
    length = sum(len(alignment[0]) for alignment in alignments)
    width = max(len(alignment) for alignment in alignments)
    array = np.full((length, width), default)
    surplus = 0
    
    for alignment in alignments:
        for j, record in enumerate(alignment):
            # Each position in sequence gets associated p-value
            for position, char in enumerate(record):
                # Position determined by leftmost value of k-mer
                # If no k-mer starts here (i.e. this pos is a gap char), 
                # default (NAN) given
                if char != '-':
                    kmer = record[position:].seq.ungap(gap='-')[:k].upper()
                    pval = dictionary.get(kmer,default)
                    # Check reverse complement
                    if pval is default:
                        kmer = kmer.reverse_complement()
                        pval = dictionary.get(kmer,default)
                    array[position+surplus][j] = pval
        surplus += len(record)
    return array

def get_hmps(df, window_size, weighted=True):
    """
    Get hmps for each sliding window.
    """
    hmps = []
    stagger = int(window_size/2)
    
    # Grab length of concatenated alignments (largest absolute position, + 1 for 0 indexing)
    l = max(df.index)[0]+1
    num_tests = len(df)
    
    # Generate start and end indices for each sliding window
    start_indices = np.array(range(0,l,stagger))
    end_indices = start_indices + window_size
    
    # Takes slice of absolute positions (index 0)
    idx = pd.IndexSlice
    if weighted:
        for start, end in zip(start_indices,end_indices):
            df_window = df.loc[idx[start:end,:,:,:],:]['p_val']
            num_tests_window = len(df_window)
            hmp = hmean(df_window)
            # Adjusted by factor of weight*(-1)
            try:
                adjusted_hmp = hmp*(num_tests/num_tests_window)
            except ZeroDivisionError:
                adjusted_hmp = 0
            hmps.append(adjusted_hmp)
    else:
        hmps = np.array([hmean(df.loc[idx[start:end,:,:,:],:]['p_val']) for start, end in zip(start_indices,end_indices)])
    return hmps

def get_hmps_array(array, window_size):
    """
    Get hmps for each sliding window.
    """
    hmps = []
    stagger = int(window_size/2)
    
    # Grab length of concatenated alignments (largest absolute position, + 1 for 0 indexing)
    #l = align_concat_len(alignments)
    l =  array.shape[0]
    
    # Generate start and end indices for each sliding window
    start_indices = np.array(range(0,l,stagger))
    end_indices = start_indices + window_size
    
    # Takes slice of absolute positions (index 0)
    hmps = np.array([hmean(array[start:end+1][~np.isnan(array[start:end+1])]) for start, end in zip(start_indices,end_indices)])
    
    return hmps

#############################################################
## Visualisation functions
#############################################################

def remove_duplicate_legends(ax=None):
    """Remove duplicate labels on a matplotlib legend.
    """
    ax = ax or plt.gca()
    
    handles, labels = ax.get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)
    
    ax.legend(newHandles, newLabels)

def plot_multiple_hmps(df, window_sizes, ax=None):
    """
    """
    for window_size in window_sizes:
        print(f"Generating visualisation for window size {window_size}bp...")
        t0 = datetime.datetime.now() 
        
        # Calculate HMP's for each window size
        hmps = get_hmps(df,window_size)
        plot_hmps(df,hmps,window_size,ax)
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    
    # Plot kmers Manhattan-plot style (scatter graph)
    print(f"Generating Manhattan-plot...")
    t0 = datetime.datetime.now()
    plot_manhattan(df)
    t1 = datetime.datetime.now()
    print(f"Total time: {(t1-t0).total_seconds()}.\n\n")

    plt.show()
    

def plot_hmps(df, hmps, window_size, ax=None):
    """Plot p-values of sliding windows vs window position across sequence.
    Optionall also plot 'Manhattan Plot' of individual kmer p-values vs
    kmer position in multi-sequence.
    """
    ax = ax or plt.gca()
    stagger = int(window_size/2)
        
    # Sliding windows plot
    for i, hmp in enumerate(hmps):
        ax.plot([i*stagger, i*stagger+window_size], [-np.log10(hmp), -np.log10(hmp)], color=colors[window_size], label=f'{window_size}bp')
    
    remove_duplicate_legends()
    ax.set_xlabel("Genome position")
    ax.set_ylabel("-log10(Adjusted p-value)")


def plot_manhattan(df, alpha=0.05, thresh=0.001, ax=None, millions=True):
    """Plot Manhattan-plot of k-mer p-values.
    """
    ax = ax or plt.gca()
    
    num_tests = len(df)
    
    # Grab bottom thresh% values by p-value (for plotting purposes)
    thresh_amount = int(thresh*num_tests)
    df_temp = df.nsmallest(thresh_amount, 'p_val')
    
    # Alternate plot colour for each alignment
    colors = df_temp.index.get_level_values('alignment') % 2
    colors = cm.Paired(colors)
    
    # Adjust p-values by weight (1/len(df))**(-1) = len(df)
    # To enable comparison with alpha
    adjusted_pvals = df_temp['p_val']*num_tests
    
    ax.scatter(df_temp.index.to_frame()['absolute_pos'], -np.log10(adjusted_pvals), edgecolors=colors, facecolors='none')
    # Plot alpha value line
    if alpha:
        ax.plot([min(df.index)[0], max(df.index)[0]], [-np.log10(alpha), -np.log10(alpha)], 'r', ls='--')
    
    # Plot x-ticks as XM instead of X000000 etc
    if millions:
        xlabels = ['{:,.2f}M'.format(x) for x in ax.get_xticks()/1000000]
        ax.set_xticklabels(xlabels)


def plot_manhattan_plotly(df,window_sizes, alpha=0.05, thresh=0.01, ax=None, millions=True):
    """Plot manhattan plot to plotly interactive graph.
    """
    num_tests = len(df)
    # Grab bottom thresh% values by p-value (for plotting purposes)
    thresh_amount = int(thresh*num_tests)
    df_temp = df.nsmallest(thresh_amount, 'p_val')
    
    ############### EDIT!!!!!!!!!!!!!!!!!!
    reference_only = False
    if reference_only:
        df_temp = df_temp[df_temp.index.get_level_values('sequence') == record_dict['MSSA476']]
    
    # Alternate colours for alignment blocks
    alignment_colors = df_temp.index.get_level_values('alignment') % 2
    
    ### INCLUDE SECTION TO HIGHLIGHT PEAKS IN ORANGE
    
    # Adjust p-values by weight (1/len(df))**(-1) = len(df)
    # To enable comparison with alpha
    adjusted_pvals = df_temp['p_val']*num_tests
    
    # Scatter graph for kmer p-values
    kmers = go.Scattergl(
        x = df_temp.index.get_level_values('absolute_pos'),
        y = -np.log10(adjusted_pvals),
        # Add kmer and sequence name text info
        text = df_temp['kmer'] + '<br>' + df_temp.index.get_level_values('sequence').map(record_dict_reverse),
        mode = 'markers',
        opacity = 0.5,
        marker = dict(color = alignment_colors, 
                      colorscale = [[0.0,'grey'],[1.0,'skyblue']]),
        name = '31-mer'
    )
    
    data = [kmers]
    
    # Adjusted alpha (Bonferroni) for comparison
    adjusted_alpha = alpha/len(df)
    alpha_trace = go.Scatter(
                x = [min(df.index)[0], max(df.index)[0]+1],
                y = [-np.log10(adjusted_alpha), -np.log10(adjusted_alpha)],
                mode = 'lines',
                line = dict(color = '#d62728',
                            dash = 'dash'),
                name='alpha')
    
    data.append(alpha_trace)
    
    # HMP windows over graph
    for window_size in window_sizes:
        hmps = get_hmps(df,window_size)
        stagger = int(window_size/2)
        # Shift necessary if not starting at basepair 0
        shift = min(df.index)[0]
        showlegend=True
        for i, hmp in enumerate(hmps):
            # Don't show duplicate legends
            if i > 0:
                showlegend = False
            windows = go.Scattergl(
                    x = [i*stagger+shift, i*stagger+window_size+shift],
                    y = [-np.log10(hmp), -np.log10(hmp)],
                    mode='lines',
                    marker = dict(color=colors[window_size]),
                    name = f'{sigfigs(window_size)} bp',
                    showlegend=showlegend)
            data.append(windows)
    
    layout = dict(title='31-mer p-values for multi-alignment of Staph a.',
                  xaxis=dict(title='Genome position'),
                  yaxis=dict(title='-log10(adjusted p-val)'))
    
    fig = dict(data=data, layout=layout)
    
    py.plot(fig, layout=layout, filename='kmer-test', auto_open=True)
    

####################################
## Main program
####################################

if __name__ == '__main__':
    
    # If debug, code will skip reordering and alignment steps
    debug = 1
    
    # kmer k
    k = 31
    
    # Sliding window size for analysing kmers
    window_sizes = (100000)
    # INSERT LOGIC TO DETERMINE STAGGER SIZE
    colors = {10:'springgreen',100:'steelblue',1000:'thistle',10000:'teal',100000:'seagreen',1000000:'turquoise'}
    
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
    
    # If dictionary in memory...
    try:
        kmer_pvalues
        print("Dictionary pre-loaded.")
    except NameError:
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
        
        # CLEAR UO MEMORY ???
        dfpvals = None
        
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
    alignment_filename = f'{base_path}alignment\\alignment.xmfa'
    
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
        align_seqs(reference, seqs, alignment_filename)
        print("Sequences aligned.")
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    
    try:
        alignments
        print("Alignments file pre-loaded.")
    except NameError:
        # Parse alignments with Biopython
        print("Loading alignment file...")
        t0 = datetime.datetime.now() 
        alignments = list(AlignIO.parse(alignment_filename, "mauve"))
        # Sort alignments to match ordering of reference genome
        sort_multi_alignment_by_reference(alignments, reference='reference')
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    
    
    ############################################
    ## P-Value/Position DataFrame
    ############################################
    
    # Dictionary for sequence ids
    try:
        record_dict
        df
        print("DataFrame pre-constructed.")
    except NameError:
        record_ids = sorted(list(set(record.id.split("/")[0].split("\\")[-1].split(".")[0] for alignment in alignments for record in alignment)))
        record_dict = {record_id:i for i, record_id in enumerate(record_ids)}
        record_dict_reverse = {i:record_id for i, record_id in enumerate(record_ids)}
        
        # Create DF
        print("Converting alignment file to DataFrame with p-value/position as row...")
        t0 = datetime.datetime.now() 
        df = alignment2df(alignments,k,kmer_pvalues)
        print("Dataframe successfully created.")
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds()}.\n\n")
    
    #plt.savefig(f'{base_path}images\\fig_{str(t1).replace(":",".")}.svg')
    
    ############################################
    ## HMPs, Windows, Visualisation
    ############################################
    
    fig, ax = plt.subplots()
    
    # Calculating window sizes - powers of 10 less than total sequence length
    total_sequence_length = max(df.index)[0]
    upper_exp = int(math.log10(total_sequence_length))+1
    lower_exp = 5
    
    window_sizes = (10**e for e in range(lower_exp,upper_exp))
    
    # Plot HMP windows and Manhattan of kmers on plot.ly
    plot_manhattan_plotly(df,window_sizes)
    
    ###############################################################
    # Stats
    ###############################################################

    percent_kmers_captured = 100*subset_proportion(df['kmer'], total_kmers)
    print(f'Proportion of kmers present in sequences tested: {percent_kmers_captured:.2f}%')
    
