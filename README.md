# kmer-hmp

## Overview of project

#### Goal:

- Develop method for applying Harmonic Mean P-Value to a multiple sequence alignment (MSA)

#### Steps:

1.  Choose twelve S. aureus genomes from a [previous analysis](http://sro.sussex.ac.uk/id/eprint/63252/1/Earle%20SG%202016.pdf) (main example of paper) conducted into fusidic acid resistance. Select files for:
    - contigs from each of the twelve genomes
    - list of k-mers from the analysis of all ~1000 genomes
    - list of associated p-values for these k-mers
2. Reorder the contigs on a genome-by-genome basis with [Mauve](http://darlinglab.org/mauve/mauve.html) with trespect to the [MSSA476 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/BX571857.1).
3. Use Mauve to align the twelve genomes into an MSA.
4. K-merize (into 31-mers) the MSA and identify the p-value corresponding to each kmer in the MSA.
5. Compute HMPs for overlapping sliding windows at different scales, e.g. 10bp, 100bp, 1kb, 10kb, 100kb, 1Mb. Suggestion: stagger each sliding window by about 50% of its length.


## Data

### /fusidic_data

This folder contains the files:

- `fusidic_acid_staph992_12_genome_subset_for_jacob_with_paths.txt`
containing sample names, guids (new sample names), contig files, md5sums for the contig files, and phenotypes for fusidic acid resistance
- `cipro_all_kmer_out.kmer.txt`
list of k-mer sequences
- `saur_992_derval_fusidic_acid_all_kmers_LMM_pvals_only.txt`
list of p-values from the LMM analysis for the above k-mers (identified on a row-number basis)
- the 12 contig files
- the reference file downloaded from Entrez


## Scripts

### `build_msa_mauve.py`

1. Reorder contigs for multiple draft files against a single reference file using Mauve on Windows.
2. Align reordered sequences with reference file using Mauve.

### `build_kmer_pval_dict.py`

3. Create dictionary of `{k-mers: p-values}` from static files

### `main.py`

4. Create a DataFrame storing the k-mer and associated p-value info for each position of each sequence in the alignment.
5. Calculate the Harmonic Mean p-value for each sliding window across the length of the sequences (for multiple window sizes)
6. Plot the k-mer and HMP p-values 'Manhattan-plot' style

#### Overview of "sliding window" logic

Each k-mer in the sequences has an associated p-value from the GWAS performed on the sequences (testing for correlation with phenotype).

This program creates a window of defined length, and travelling laterally across the aligned sequences, takes the Harmonic Mean of the p-values of the k-mers in said window.

Example:
- 2 alignments (first contianing 3 sequences, second containing 2)
- window size = 4
- k = 3

![Sliding windows animation](https://github.com/ja-ox/kmer-hmp/blob/master/images/sliding_windows_animation.gif)

The 3-mers included in each window in the animation are:

1. ACG, ACG, ACG, CGT, CGA, CGG
2. CGT, CGA, CGG, GTG, GAG, GGC
3. GTG, GAG, GGC, TGC, AGC, GCC
4. TGC, AGC, GCC, GCA, GCA, CCA
5. GCA, GCA, CCA
6.
7. TGT, TGT
8. TGT, TGT, GTA, GTA

Thus, in e.g. window 1 the HMP of the 6 p-values corresponding to the 6 k-mers it includes will be calculated and assigned position 1.
    
#### Handling of gap characters 
    
If there are gap characters in a sequence, e.g.:

    ACA-CACGTG

the program will ignore the gap characters, thus the k-mers returned for this window (of size 10) for this sequence would be:

    ACACA
     CACAC
      ACACG
      
       CACGT
        ACGTG

Note that no value is associated to position 4.

#### What happens when a window overlaps multiple alignments?

If there are multiple alignments ordered sequentially, some sliding windows will overlap with the discontinuity between alignments.

         [...  ...]
    [ACACGTGT][CGAT]

In these instances k-mers are taken from each alignment intersection only.

i.e. for `k=2`, `window_size=6` in the above example:

          TG   CG
           GT   GA

K-mers 'across' neighbouring alignments are not considered, as these may not correspond to real k-mers present in the genome (since the ordering of the alignment may not be reflected in reality). E.g. in this instance, the (possibly fictional) inter-alignment k-mer `TC` is not recognised.

#### What happens when the overlapped alignments have different depths (numbers of sequences in them)?
    
Some neighbouring alignments may have different numbers of aligned sequences in them. e.g. 3 and 2:

     [...  ...]
    [ACTG][GTGTA]
    [CCTG][GTGGA]
    [TCTG]

In these instances, if a window overlaps an alignment boundary, k-mers are taken from all aligned sequences. e.g. k=3, window_size=6:

     CTG  GTG
     CTG  GTG
     CTG


## Mauve instructions

>#### [Reordering contigs from the command-line (batch mode)](http://darlinglab.org/mauve/user-guide/reordering.html) 
>In situations where it is necessary to order contigs in a large number of draft genomes it is often more desirable to automate the process using command-line interfaces and scripts. Mauve Contig Mover supports command-line operation through the Mauve Java JAR file.
>
>Given a reference genome file called `reference.gbk` and a draft genome called `draft.fasta`, one would invoke the reorder program with the following syntax:
>
>      java -Xmx500m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output results_dir -ref reference.gbk -draft draft.fasta
>
>The file Mauve.jar is part of the Mauve distribution. On windows systems it can usually be found in `C:\Program Files\Mauve X\Mauve.jar` where `X` is the version of Mauve. On Mac OS X it is located inside the Mauve application. For example, if Mauve has been placed in the OS X applications folder, Mauve.jar can be found at `/Applications/Mauve.app/Contents/Resources/Java/Mauve.jar`. 
>On Linux, `Mauve.jar` is simply at the top level of the `tar.gz` archive. In the above example command, it will be necessary to specify the full path to the `Mauve.jar` file.
