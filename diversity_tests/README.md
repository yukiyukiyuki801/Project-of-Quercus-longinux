# Pairwise diversity tests

The analysis uses unpaired, two-sided Wilcoxon rank-sum tests from the R
package `stats` to compare the three varieties for He, Ho, FIS, and nucleotide
diversity. The three pairwise P values within each metric are adjusted with
the Holm method. The core test is unchanged from the original analysis; the
revision corrects the erroneous label "signed-rank" to "rank-sum" and adds
the requested multiple-testing correction.

```bash
module load R/4.4.1-foss-2022b
Rscript diversity_tests/scripts/pairwise_wilcoxon.R \
  diversity_tests/data/diversity_by_population.tsv \
  diversity_tests/work/pairwise_wilcoxon_rank_sum.tsv
```

The archived result is in `outputs/pairwise_wilcoxon_rank_sum.tsv`.
