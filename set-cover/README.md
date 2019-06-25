#Set Cover Problem

Algorithm to select smallest subset of genomes from GWAS which contain > x% of total unique k-mers appearing in complete set of genomes.

Implementation: greedy algorithm for solving [Set Cover problem](https://en.wikipedia.org/wiki/Set_cover_problem#Greedy_algorithm):

> Inapproximability results show that the greedy algorithm is essentially the best-possible polynomial time approximation algorithm for set cover up to lower order terms (see Inapproximability results below), under plausible complexity assumptions.   

Current best (limited by RAM to only considering ~450/992 genomes in any one iteration):

- 01 genomes: 56.7%
- 04 genomes: 80.2%
- 11 genomes: 90.3%
- 30 genomes: 95.1%

Considering slice of length 450: `genomes[406:856]`

![Graph of coverage vs number of genomes](/images/set_cover_graph.png)
