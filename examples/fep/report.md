## Project goal:

- Find an example of a k-mer GWAS where 'by visual inspection' there are some regions of probable significance but they fall just short of the Bonferroni corrected significance threshold.
- Calculate Harmonic Mean P-values of sliding windows of various sizes across Multi-Sequence-Alignment of genomes from study (or possible of multi-aligned genes, with sequences split and collated on a gene-by-gene basis).
- Determine whether windows of x size on regions of possible 'regional' significance are actually significant.

## 1. Data-set

- Method tested on FUC resistant S. Aureus from study ...... Pipeline code written.
- Chosen FEP resistant E. Coli from study ...... for actual analysis since it appears to match criteria above. I.e. consider the SNP Manhattan plot below:

- Blue dots = SNP adjusted p-values below threshold
- Red dots = SNP adjusted p-values above threshold
- Cyan bars = HMP values below threshold (2kbp)
- Gold bars = HMP values above threshold (2kbp)
- Red circles = regions where HMP exceeds threshold but no individual p-value in region does

## 2. K-mer analysis

Running into issues when testing code on k-mer data from study. Data appears to have significant 'banding' obscuring peaks of significance:

Unsure what is causing such large disparity between SNP and k-mer based plots. Further investigation needed. Have tested using `p_score` and `p_lrt` values for p-values in data.

Possible reason: p-values are not matched to correct k-mer in initial dictionary.
