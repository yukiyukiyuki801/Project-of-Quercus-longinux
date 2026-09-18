# Pairwise diversity tests

This workflow performs two-sided Mann-Whitney/Wilcoxon rank-sum tests between
the three varieties for He, Ho, FIS, and nucleotide diversity, then applies
Holm correction within each metric.

```bash
python3 diversity_tests/scripts/verify_diversity_tests.py \
  --input diversity_tests/data/diversity_by_population.tsv \
  --output diversity_tests/work/pairwise_wilcoxon_rank_sum.tsv
```

The archived result is in `outputs/pairwise_wilcoxon_rank_sum.tsv`.
