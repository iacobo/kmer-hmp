# kmer-hmp
Repo for code and data related to project developing methodology for applying Harmonic Mean P-Values to sliding windows of multiple aligned genomes.

Steps:

1.  Choose twelve S. aureus genomes from a [previous analysis](http://sro.sussex.ac.uk/id/eprint/63252/1/Earle%20SG%202016.pdf) (main example of paper) conducted into fusidic acid resistance. Select contigs from the twelve genomes, and a list of kmers and associated p-values from the analysis of all ~1000 genomes.
2. Using the [MSSA476 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/BX571857.1), reorder the contigs on a genome-by-genome basis with [Mauve](http://darlinglab.org/mauve/mauve.html).
3. Use Mauve to align the twelve genomes.
4. Kmerize (into 31mers) the multiple sequence alignment (MSA) and identify the p-value corresponding to each kmer in the MSA.
5. Compute HMPs for overlapping sliding windows at different scales, e.g. 10bp, 100bp, 1kb, 10kb, 100kb, 1Mb. Suggestion: stagger each sliding window by about 50% of its length.


### /fusidic_data

This folder contains the files:

- `fusidic_acid_staph992_12_genome_subset_for_jacob_with_paths.txt`
containing sample names, guids (new sample names), contig files, md5sums for the contig files, and phenotypes for fusidic acid resistance
- `fuc_LMM_results_400000kmers_out.txt`
containing p-values from the LMM analysis for the 400,000 kmers 
- `fuc_400000_kmers.txt`
the kmer sequences for the p-values in the above file (i.e. the kmer in line `n` of this file has associated p-value from line `n+1` (due to the header line) in the previous file)
- plus the 12 contig files.
