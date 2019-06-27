# Set Cover Problem

Algorithm to select smallest subset of genomes from GWAS which contain > x% of total unique k-mers appearing in complete set of genomes.

Implementation: greedy algorithm for solving [Set Cover problem](https://en.wikipedia.org/wiki/Set_cover_problem#Greedy_algorithm):

> Inapproximability results show that the greedy algorithm is essentially the best-possible polynomial time approximation algorithm for set cover up to lower order terms (see Inapproximability results below), under plausible complexity assumptions.   

Since we are dealing with k-mer **pattern** files, this becomes a **weighted** set cover problem, where the weight of each pattern is the number of k-mers which share that pattern.

- 01 genomes: 19.2%
- 02 genomes: 29.2%
- 05 genomes: 43.0%
- 10 genomes: 54.7%
- 20 genomes: 66.9%
- 50 genomes: 79.9%

![Graph of coverage vs number of genomes](/images/set_cover_graph.png)
