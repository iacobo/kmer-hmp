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

    Covering set row numbers: [139, 364, 86, 641, 696, 168, 980, 372, 578, 846, 676, 82, 197, 162, 370, 109, 754, 533, 525, 163, 765, 84, 502, 934, 875, 698, 680, 79, 701, 229, 238, 849, 883, 751, 650, 165, 188, 198, 960, 695, 433, 561, 87, 337, 486, 504, 806, 600, 267, 874, 83]
