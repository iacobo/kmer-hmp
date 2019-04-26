# kmer-hmp
Repo for code and data related to project developing methodology for applying Harmonic Mean P-Values to sliding windows of multiple aligned genomes.

Steps:

1.  Choose twelve S. aureus genomes from a previous analysis conducted into fucidic acid resistance. Select contigs from the twelve genomes, and a list of kmers and associated p-values from the analysis of all ~1000 genomes.
2. Using the MSSA476 reference genome, reorder the contigs on a genome-by-genome basis with Mauve.
3. Use Mauve to align the twelve genomes.
4. Kmerize (into 31mers) the multiple sequence alignment (MSA) and identify the p-value corresponding to each kmer in the MSA.
5. Compute HMPs for overlapping sliding windows at different scales, e.g. 10bp, 100bp, 1kb, 10kb, 100kb, 1Mb. Suggestion: stagger each sliding window by about 50% of its length.
