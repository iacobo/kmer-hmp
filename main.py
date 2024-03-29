# -*- coding: utf-8 -*-
"""
    1. Create a DataFrame storing the k-mer and associated p-value info
       for each position of each sequence in the alignment.
    2. Calculate the Harmonic Mean p-value for each sliding window
       across the length of the sequences (for multiple window sizes)
    3. Plot the k-mer and HMP p-values 'Manhattan-plot' style
"""

import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hmean

import plotly.plotly as py
import plotly.graph_objs as go

import build_msa_mauve
import build_kmer_pval_dict


def subset_proportion(subset, completeset):
    """Given `subset`, return the relative size of the `subset` compared to
    `completeset`.
    Raise error if `subset` not a subset of `completeset`.
    """
    subset = set(subset)
    completeset = set(completeset)
    len_intersection = len(subset.intersection(completeset))
    
    assert len_intersection == len(subset), "First element is not a subset of second element."
    
    return len_intersection/len(completeset)

def alignment2df(alignments,k,dictionary,record_dict,default=np.nan):
    """Create DataFrame from multi-alignment.
    
    Given MSA and dictionary of form {k-mer:p-value} create df of form:
        
                                                      kmer      pval
    absolute_pos alignment sequence position                                           
    0            0         0        0        TAATCGGACCTGG  0.876878
                           1        0        TAATCGGACCCGG  0.840357
    1            0         0        1        AATAGGACCTGGT  0.792824
                           2        1        AATCGGACCTGGT  0.897792
                           3        1        AATCGGACCTGGT  0.876878
    """
    
    columns = ['kmer','pval','position','absolute_pos','alignment','sequence']
    
    temp = []
    
    for i, alignment in enumerate(alignments):
        surplus = sum(len(alignment[0]) for alignment in alignments[:i])
        for record in alignment:
            # id of format 'C:\Users\User\location\genomes\reference_genome\Seq_Id/xxx-yyy'
            # Grabs 'Seq_Id' part
            seq_id = Path(record.id).parent.stem
            seq_id = record_dict[seq_id]
            # Each position in sequence gets associated p-value
            for position, char in enumerate(record):
                # Fetch k-mer starting at this position
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

def get_hmps(df, window_size, weighted=True):
    """
    Create array of Harmonic Mean p-values for sliding windows across
    the dataframe.
    
    If `weighted`, multiplies each value by a factor of:
        
        total values in dataframe / number of values in window
    """
    hmps = []
    stagger = int(window_size/2)
    
    # Grab length of concatenated alignments (largest absolute position)
    first = min(df.index)[0]
    last = max(df.index)[0]
    num_tests = len(df)
    
    # Generate start and end indices for each sliding window
    start_indices = np.array(range(first,last+1,stagger))
    end_indices = start_indices + window_size
    
    # Takes slice of absolute positions (index 0)
    idx = pd.IndexSlice
    if weighted:
        for start, end in zip(start_indices,end_indices):
            df_window = df.loc[idx[start:end,:,:,:],:]['pval']
            num_tests_window = len(df_window)
            hmp = hmean(df_window)
            # Adjusted by factor of weight**(-1)
            try:
                adjusted_hmp = hmp*(num_tests/num_tests_window)
            except ZeroDivisionError:
                adjusted_hmp = np.nan
            hmps.append(adjusted_hmp)
    else:
        hmps = np.array([hmean(df.loc[idx[start:end,:,:,:],:]['pval']) for start, end in zip(start_indices,end_indices)])
    return hmps


def plot_manhattan_plotly(df, window_sizes, record_dict_reverse, alpha=0.05, thresh=0.01):
    """Plot manhattan plot to plotly interactive graph.
    Alternate alignments coloured blue/grey.
    """
    num_tests = len(df)
    # Grab bottom thresh% values by p-value (for plotting purposes)
    thresh_amount = int(thresh*num_tests)
    df_temp = df.nsmallest(thresh_amount, 'pval')
    
    # Alternate colours for alignment blocks
    alignment_colors = df_temp.index.get_level_values('alignment') % 2
    
    ###TODO: INCLUDE SECTION TO HIGHLIGHT PEAKS IN ORANGE
    
    # Adjust p-values by weight (1/len(df))**(-1) = len(df)
    # To enable comparison with alpha
    adjusted_pvals = df_temp['pval']*num_tests
    
    # Adjusted alpha (Bonferroni) for comparison
    adjusted_alpha = alpha/len(df)
    y = -np.log10(adjusted_alpha)
    alpha_trace = go.Scatter(x = [min(df.index)[0], max(df.index)[0]+1],
                             y = [y, y],
                             name = 'alpha',
                             mode = 'lines',
                             line = dict(color = '#d62728',
                                         dash = 'dash'))
    
    data = [alpha_trace]
    
    # Scatter graph for kmer p-values
    kmers = go.Scattergl(x = df_temp.index.get_level_values('absolute_pos'),
                         y = -np.log10(adjusted_pvals),
                         # Add kmer and sequence name text info
                         name = f'{k}-mer',
                         text = df_temp['kmer'] + '<br>' + df_temp.index.get_level_values('sequence').map(record_dict_reverse),
                         mode = 'markers',
                         opacity = 0.5,
                         marker = dict(color = alignment_colors,
                                       colorscale = [[0.0,'grey'],[1.0,'skyblue']]))
    
    data.append(kmers)
    
    # Colour-map for different window_sizes
    # +1 to ensure start color != end color
    colors = [f'hsl({h},50%,50%)' for h in np.linspace(0, 360, len(window_sizes)+1)]
    
    # HMP windows over graph
    for j, window_size in enumerate(window_sizes):
        stagger = int(window_size/2)
        hmps = get_hmps(df,window_size)
        name = f'{sigfigs(window_size)} bp'
        # Shift necessary if df doesn't start at basepair 0
        shift = min(df.index)[0]
        showlegend = True
        for i, hmp in enumerate(hmps):
            windows = go.Scatter(x = [i*stagger+shift, i*stagger+window_size+shift],
                                 y = [-np.log10(hmp), -np.log10(hmp)],
                                 mode = 'lines',
                                 line = dict(color=colors[j]),
                                 name = name,
                                 legendgroup = name,
                                 showlegend = showlegend)
            data.append(windows)
            # Hide duplicate legends
            showlegend = False
    
    layout = dict(title=f'{k}-mer p-values for multi-alignment of Staph a.',
                  xaxis=dict(title='Genome position'),
                  yaxis=dict(title='-log10(adjusted p-val)'))
    
    fig = dict(data=data, layout=layout)
    py.plot(fig, layout=layout, filename='kmer-test', auto_open=True)

def sigfigs(n):
    """Return string formatted integer with K or M for thousands, millions
    sig figs e.g:

    1000    > 1K
    5000000 > 5M
    """
    if n >= 1000000: 
        if (n/1000000).is_integer():
            return f'{int(n/1000000)}M'
        else:
            return f'{n/1000000:.1f}M'
    elif n >= 1000: 
        if (n/1000).is_integer():
            return f'{int(n/1000)}K'
        else:
            return f'{n/1000:.1f}K'
    else: 
        return f'{n}'


def main(k, alignments, kmer_pvalues):
    
    # Dictionary for sequence ids
    record_ids = sorted(list(set(Path(record.id).parent.stem for alignment in alignments for record in alignment)))
    record_dict = {record_id:i for i, record_id in enumerate(record_ids)}
    record_dict_reverse = {i:record_id for i, record_id in enumerate(record_ids)}
    
    # Create DF
    print('Converting alignment file to DataFrame with p-value/position as row...')
    t0 = datetime.datetime.now() 
    df = alignment2df(alignments,k,kmer_pvalues,record_dict)
    print('Dataframe successfully created.')
    t1 = datetime.datetime.now()
    print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')
    
    # Calculate window sizes - powers of 10 less than total sequence length
    total_sequence_length = max(df.index)[0]
    upper_exp = int(np.log10(total_sequence_length))+1
    lower_exp = 5
    
    window_sizes = [10**e for e in range(lower_exp,upper_exp)]
    
    # Plot HMP windows and Manhattan of kmers on plot.ly
    print('Calculating HMP and plotting data to plot.ly...')
    t0 = datetime.datetime.now() 
    plot_manhattan_plotly(df,window_sizes,record_dict_reverse)
    t1 = datetime.datetime.now()
    print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')

    percent_kmers_captured = 100*subset_proportion(df['kmer'], kmer_pvalues.keys())
    print(f'Proportion of k-mers present in sequences tested: {percent_kmers_captured:.2f}%')
    
    return df


if __name__ == '__main__':
    mauve_dir = Path('C:/Program Files (x86)/Mauve 20150226')
    base_path = Path('C:/Users/Jacob/Downloads/fusidic_data')
    reference = base_path / 'genomes/reference_genome/Record_49484912.fasta'
    drafts_dir = base_path / 'genomes/draft_genomes'
    kmers = base_path / 'static_files/fusidic_acid_kmers.txt'
    pvals = base_path / 'static_files/fusidic_acid_pvals.txt'
    
    print('Loading k-mer/p-values dictionary...')
    t0 = datetime.datetime.now() 
    kmer_pvalues = build_kmer_pval_dict.main(kmers, pvals)
    t1 = datetime.datetime.now()
    print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')
    
    k = len(list(kmer_pvalues.keys())[0])
    alignments = build_msa_mauve.main(base_path=base_path / 'genomes', reference=reference, drafts_dir=drafts_dir, mauve_dir=mauve_dir)
    df = main(k, alignments, kmer_pvalues)
