## Project goal:

- Find an example of a k-mer GWAS where 'by visual inspection' there are some regions of probable significance but they fall just short of the Bonferroni corrected significance threshold.
- Calculate Harmonic Mean P-values of sliding windows of various sizes across Multi-Sequence-Alignment of genomes from study (or possible of multi-aligned genes, with sequences split and collated on a gene-by-gene basis).
- Determine whether windows of x size on regions of possible 'regional' significance are actually significant.

## 1. Data-set

- Method tested on FUC resistant S. Aureus from study ...... Pipeline code written.
- Chosen FEP resistant E. Coli from study ...... for actual analysis since it appears to match criteria above. I.e. consider the SNP Manhattan plot below:

![SNP manhattan plot](/examples/fep/image_2000_2.png)

- Blue dots = SNP adjusted p-values below threshold
- Red dots = SNP adjusted p-values above threshold
- Cyan bars = HMP values below threshold (2kbp)
- Gold bars = HMP values above threshold (2kbp)
- Red circles = regions where HMP exceeds threshold but no individual p-value in region does

## 2. K-mer analysis

Using a set of 4-12 genomes from set of ~200, and the subset of 400,000 k-mers for which the p-values had been calculated in the GWAS, ran previous code from FUC project to produce the following plot:

![k-mer manhattan plot](/examples/fep/Figure_2.png)

**Todo:** create HMP's, plot alpha, test thresholds.
