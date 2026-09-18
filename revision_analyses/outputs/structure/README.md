# STRUCTURE replication and optimal-K summary

- Input: 180 individuals and 1,291 LD-pruned RAD SNPs.
- Runs: K=1-7, 10 independent replicates per K, 25,000 burn-in iterations and 100,000 MCMC iterations.
- Maximum Evanno Delta K: **K=3** (Delta K = **174.613**).
- Replicate cluster labels were aligned to the highest-likelihood run at each K before averaging.

## Mean ancestry by named variety at K=3

- kuoi (n=20): C1=0.011, C2=0.047, C3=0.941
- lativiolaciifolia (n=57): C1=0.317, C2=0.667, C3=0.016
- longinux (n=103): C1=0.343, C2=0.613, C3=0.044

## Interpretation

Evanno Delta K identifies the strongest upper level of hierarchical structure and is undefined for K=1, so it should be interpreted together with the mean log likelihoods, replicate variability, ancestry plots, geography, and biological question. The alignment here uses minimum-distance assignment of cluster columns, a reproducible analogue of the label-switching correction performed by CLUMPP. CLUMPP aligns replicate Q matrices; it does not calculate Delta K.

The available VCF is a filtered, LD-pruned subset. These runs reproduce broad ancestry structure for the revision but should not be described as an independent genome-wide discovery analysis.
