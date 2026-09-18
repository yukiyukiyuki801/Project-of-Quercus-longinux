# ABBA-BABA / D-statistic sensitivity analysis

## Result

- Topology tested: `((lativiolaciifolia,longinux),kuoi)`, rooted with *Quercus glauca* DH01.
- Patterson's D = **0.0147**; block-jackknife Z = **1.012**; two-sided P = **0.3114**.
- f4-ratio = **0.0573**.
- Weighted site-pattern counts: ABBA = 23.98; BABA = 23.29.
- Q. glauca genotypes passing DP >= 3: **1060 / 1291** (82.1%).

The test does not detect significant asymmetric allele sharing among the three varieties.

## Outgroup-depth sensitivity

| Minimum DH01 depth | Usable loci | D | Z | P |
|---:|---:|---:|---:|---:|
| 3 | 1060 | 0.0147 | 1.012 | 0.3114 |
| 5 | 983 | 0.0194 | 1.316 | 0.1882 |
| 10 | 854 | 0.0244 | 1.480 | 0.1388 |

The inference remains nonsignificant under stricter minimum-depth thresholds.

## Interpretation limits

This is an exploratory, dataset-wide sensitivity test using only the 1,291 LD-pruned, variant-only RAD SNPs retained in the published ingroup VCF. It is not a genome-wide or windowed introgression scan. The two non-kuoi varieties are not reciprocally monophyletic in SVDquartets, so asymmetry—if present—could reflect population structure, incomplete lineage sorting, reference/mapping bias, or introgression. A nonsignificant result is not proof of no historical gene flow. One RAD-seq outgroup individual and sparse informative site patterns further limit power. This result should therefore be presented as a supplementary robustness check, not as decisive evidence for or against reinforcement or species status.

## Reproducibility

- Outgroup reads: public ddRAD-seq sample *Q. glauca* DH01 (SRS19813215; run SRR27151155; SbfI/MspI).
- Mapping reference: *Q. robur* PM1N, the same assembly coordinates used by the ingroup VCF.
- Mapping/calling: BWA-MEM followed by samtools/bcftools, mapping and base quality >= 20; outgroup DP >= 3.
- Statistic: Dsuite Dtrios v0.5 r58, expected variety topology, 20 block-jackknife blocks.
