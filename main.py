# -*- coding: utf-8 -*-
"""
    1. Create a DataFrame storing the k-mer and associated p-value info
       for each position of each sequence in the alignment.
    2. Calculate the Harmonic Mean p-value for each sliding window
       across the length of the sequences (for multiple window sizes)
    3. Plot the k-mer and HMP p-values 'Manhattan-plot' style
"""

import os
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hmean

import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

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

def get_hmps(df, window_size, weighted=True):
    """
    Create array of Harmonic Mean p-values for each sliding window across
    the dataframe of given window size.
    
    If `weighted`, multiplies each value by a factor of:
        
        total values in dataframe / number of values in window
    """
    hmps = []
    stagger = int(window_size/2)
    
    # Grab length of concatenated alignments (largest absolute position, + 1 for 0 indexing)
    length = max(df.index)[0]+1
    first = min(df.index)[0]
    num_tests = len(df)
    
    # Generate start and end indices for each sliding window
    start_indices = np.array(range(first,length,stagger))
    end_indices = start_indices + window_size
    
    # Takes slice of absolute positions (index 0)
    idx = pd.IndexSlice
    if weighted:
        for start, end in zip(start_indices,end_indices):
            df_window = df.loc[idx[start:end,:,:,:],:]['p_val']
            num_tests_window = len(df_window)
            hmp = hmean(df_window)
            # Adjusted by factor of weight**(-1)
            try:
                adjusted_hmp = hmp*(num_tests/num_tests_window)
            except ZeroDivisionError:
                adjusted_hmp = 0
            hmps.append(adjusted_hmp)
    else:
        hmps = np.array([hmean(df.loc[idx[start:end,:,:,:],:]['p_val']) for start, end in zip(start_indices,end_indices)])
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
    """ Call plot_hmps for multiple values of window_size.
    Automatically generate distinct colours for different plots.
    """
    
    for i, window_size in enumerate(window_sizes):
        print(f"Generating visualisation for window size {window_size}bp...")
        t0 = datetime.datetime.now() 
        
        # Calculate HMP's for each window size
        hmps = get_hmps(df,window_size)
        color = plt.get_cmap('jet')(i/len(window_sizes))
        plot_hmps(hmps,window_size,color,ax)
        t1 = datetime.datetime.now()
        print(f"Total time: {(t1-t0).total_seconds():.2f}s.\n")
    
    # Plot k-mers Manhattan-plot style (scatter graph)
    print(f"Generating Manhattan-plot...")
    t0 = datetime.datetime.now()
    plot_manhattan(df, ax=ax)
    t1 = datetime.datetime.now()
    print(f"Total time: {(t1-t0).total_seconds():.2f}s.\n")
    

def plot_hmps(hmps, window_size, color, ax=None):
    """Plot p-values of sliding windows vs window position across sequence.
    Optionall also plot 'Manhattan Plot' of individual kmer p-values vs
    kmer position in multi-sequence.
    """
    ax = ax or plt.gca()
    stagger = int(window_size/2)
    
    # Sliding windows plot
    for i, hmp in enumerate(hmps):
        ax.plot([i*stagger, i*stagger+window_size], 
                [-np.log10(hmp), -np.log10(hmp)], 
                c = color, 
                label = f'{window_size}bp')
    
    remove_duplicate_legends()
    ax.set_xlabel('Genome position')
    ax.set_ylabel(r'$-\log_{10}($Adjusted $p$-value$)$')


def plot_manhattan(df, alpha=0.05, thresh=0.001, ax=None, millions=True):
    """Plot Manhattan-plot of k-mer p-values.
    """
    ax = ax or plt.gca()
    
    num_tests = len(df)
    
    # Grab bottom thresh% values by p-value (for plotting purposes)
    thresh_amount = int(thresh*num_tests)
    df_small = df.nsmallest(thresh_amount, 'p_val')
    
    # Alternate plot colour for each alignment
    colors = df_small.index.get_level_values('alignment') % 2
    colors = cm.Paired(colors)
    
    # Adjust p-values by weight (1/len(df))**(-1) = len(df)
    # To enable comparison with alpha
    adjusted_pvals = df_small['p_val']*num_tests
    
    ax.scatter(df_small.index.to_frame()['absolute_pos'], -np.log10(adjusted_pvals), edgecolors=colors, facecolors='none')
    # Plot alpha value line
    if alpha:
        ax.plot([min(df.index)[0], max(df.index)[0]], [-np.log10(alpha), -np.log10(alpha)], 'r', ls='--')
    
    # Plot x-ticks as XM instead of X000000 etc
    if millions:
        xlabels = ['{:,.2f}M'.format(x) for x in ax.get_xticks()/1000000]
        ax.set_xticklabels(xlabels)


def plot_manhattan_plotly(df, window_sizes, record_dict_reverse, alpha=0.05, thresh=0.01):
    """Plot manhattan plot to plotly interactive graph.
    Alternate alignments coloured blue/grey.
    """
    num_tests = len(df)
    # Grab bottom thresh% values by p-value (for plotting purposes)
    thresh_amount = int(thresh*num_tests)
    df_temp = df.nsmallest(thresh_amount, 'p_val')
    
    # Alternate colours for alignment blocks
    alignment_colors = df_temp.index.get_level_values('alignment') % 2
    
    ###TODO: INCLUDE SECTION TO HIGHLIGHT PEAKS IN ORANGE
    
    # Adjust p-values by weight (1/len(df))**(-1) = len(df)
    # To enable comparison with alpha
    adjusted_pvals = df_temp['p_val']*num_tests
    
    # Adjusted alpha (Bonferroni) for comparison
    adjusted_alpha = alpha/len(df)
    y = -np.log10(adjusted_alpha)
    alpha_trace = go.Scatter(
                x = [min(df.index)[0], max(df.index)[0]+1],
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
                         name = '31-mer',
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
        hmps = get_hmps(df,window_size)
        stagger = int(window_size/2)
        name = f'{sigfigs(window_size)} bp'
        # Shift necessary if not starting at basepair 0
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
    
    layout = dict(title='31-mer p-values for multi-alignment of Staph a.',
                  xaxis=dict(title='Genome position'),
                  yaxis=dict(title='-log10(adjusted p-val)'))
    
    fig = dict(data=data, layout=layout)
    
    py.plot(fig, layout=layout, filename='kmer-test', auto_open=True)

def sigfigs(n):
    """Return string formatted integer with K or M for thousands, millions
    sig figs.
    
    E.g. 1000 > 1K
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


####################################
## Main program
####################################

def main(k, alignments, kmer_pvalues, df=None):
    
    ############################################
    ## P-Value/Position DataFrame
    ############################################
    
    # Dictionary for sequence ids
    record_ids = sorted(list(set(record.id.split('/')[0].split('\\')[-1].split('.')[0] for alignment in alignments for record in alignment)))
    record_dict = {record_id:i for i, record_id in enumerate(record_ids)}
    record_dict_reverse = {i:record_id for i, record_id in enumerate(record_ids)}
    
    if df is None:
        # Create DF
        print('Converting alignment file to DataFrame with p-value/position as row...')
        t0 = datetime.datetime.now() 
        df = alignment2df(alignments,k,kmer_pvalues,record_dict)
        print('Dataframe successfully created.')
        t1 = datetime.datetime.now()
        print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')
    
    ############################################
    ## HMPs, Windows, Visualisation
    ############################################
    
    # Calculating window sizes - powers of 10 less than total sequence length
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
    
    ###############################################################
    ## Stats
    ###############################################################

    percent_kmers_captured = 100*subset_proportion(df['kmer'], kmer_pvalues.keys())
    print(f'Proportion of k-mers present in sequences tested: {percent_kmers_captured:.2f}%')
    
    return df


if __name__ == '__main__':
    mauve_dir = Path('C:/Program Files (x86)/Mauve 20150226')
    base_path = Path('C:/Users/Jacob/Downloads/fusidic_data')
    reference = base_path / 'genomes/reference_genome/Record_49484912.fasta'
    kmers = base_path / 'static_files/fusidic_acid_kmers.txt'
    pvals = base_path / 'static_files/fusidic_acid_pvals.txt'
    
    os.chdir(base_path)
    
    print('Loading k-mer/p-values dictionary...')
    t0 = datetime.datetime.now() 
    kmer_pvalues = build_kmer_pval_dict.main(kmers, pvals)
    t1 = datetime.datetime.now()
    print(f'Total time: {(t1-t0).total_seconds():.2f}s.\n')
    
    alignments = build_msa_mauve.main(base_path=base_path / 'genomes', reference=reference, mauve_dir=mauve_dir)
    df = main(31, alignments, kmer_pvalues)
