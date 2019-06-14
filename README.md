# kmer-hmp
Repo for code and data related to project developing methodology for applying Harmonic Mean P-Values to sliding windows of multiple aligned genomes.

Steps:

1.  Choose twelve S. aureus genomes from a [previous analysis](http://sro.sussex.ac.uk/id/eprint/63252/1/Earle%20SG%202016.pdf) (main example of paper) conducted into fusidic acid resistance. Select contigs from the twelve genomes, and a list of kmers and associated p-values from the analysis of all ~1000 genomes.
2. Using the [MSSA476 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/BX571857.1), reorder the contigs on a genome-by-genome basis with [Mauve](http://darlinglab.org/mauve/mauve.html).
3. Use Mauve to align the twelve genomes.
4. Kmerize (into 31mers) the multiple sequence alignment (MSA) and identify the p-value corresponding to each kmer in the MSA.
5. Compute HMPs for overlapping sliding windows at different scales, e.g. 10bp, 100bp, 1kb, 10kb, 100kb, 1Mb. Suggestion: stagger each sliding window by about 50% of its length.

---


## /fusidic_data

This folder contains the files:

- `fusidic_acid_staph992_12_genome_subset_for_jacob_with_paths.txt`
containing sample names, guids (new sample names), contig files, md5sums for the contig files, and phenotypes for fusidic acid resistance
- `fuc_LMM_results_400000kmers_out.txt`
containing p-values from the LMM analysis for the 400,000 kmers 
- `fuc_400000_kmers.txt`
the kmer sequences for the p-values in the above file (i.e. the kmer in line `n` of this file has associated p-value from line `n+1` (due to the header line) in the previous file)
- plus the 12 contig files.

---

## main.py

Script designed to:
    
1. Reorder contigs for multiple draft files against a single 
   reference file using Mauve on Windows.
2. Align reordered sequences with reference file using Mauve.
3. Create a DataFrame storing the k-mer and associated p-value info
   for each position of each sequence in the alignment.
3. Calculate the Harmonic Mean p-value for each sliding window
   across the length of the sequences (for multiple window sizes)
4. Plot the k-mer and HMP p-values 'Manatan-plot' style


### Overview of objects in main.py
    
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

### Explanation of logic in main.py
    
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
    
### Handling of gap characters 
    
If there are gap characters in a sequence, e.g.:

    ACACAC-GTG

Within a given window, the program will ignore the gap characters, thus the
k-mers returned for this window (of size 10) for this sequence would be:

    ACACA
     CACAC
      ACACG
       CACGT
        ACGTG

### What happens when a window overlaps multiple alignments?

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

### What happens when the overlapped alignments have different depths (numbers of sequences in them)?
    
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

---

## Mauve instructions

>#### Reordering contigs from the command-line (batch mode)
>In situations where it is necessary to order contigs in a large number 
of draft genomes it is often more desirable to automate the process 
using command-line interfaces and scripts. Mauve Contig Mover supports 
command-line operation through the Mauve Java JAR file.
>
>Given a reference genome file called “reference.gbk” and a draft genome 
called “draft.fasta”, one would invoke the reorder program with the 
following syntax:
>
>    java -Xmx500m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output results_dir -ref reference.gbk -draft draft.fasta
>
>The file Mauve.jar is part of the Mauve distribution. On windows 
systems it can usually be found in C:\Program Files\Mauve X\Mauve.jar 
where X is the version of Mauve. On Mac OS X it is located inside the 
Mauve application. For example, if Mauve has been placed in the OS X 
applications folder, Mauve.jar can be found at 
/Applications/Mauve.app/Contents/Resources/Java/Mauve.jar. 
>On Linux, Mauve.jar is simply at the top level of the tar.gz archive. 
In the above example command, it will be necessary to specify the full 
path to the Mauve.jar file.

-  http://darlinglab.org/mauve/user-guide/reordering.html
