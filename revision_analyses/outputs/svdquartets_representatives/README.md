# SVDquartets result summary

- Input: 27 individuals (one lowest-missing representative per ingroup population plus DH01), treated as 27 tips, and 1,291 LD-pruned RAD SNPs.
- Method: PAUP* SVDquartets, all 17,550 possible quartets, and 100 standard bootstrap replicates.
- Var. *kuoi* populations (GS, LZ, SK) are monophyletic: **True**.
- Bootstrap support for the complete var. *kuoi* population clade: **62%**.
- The figure displays bootstrap support for every internal split. The split adjacent to the *Q. glauca* root has **11%** support; the post-hoc root itself has no separate bootstrap value.
- The combined var. *longinux* + var. *lativiolaciifolia* group is monophyletic after rooting on *Q. glauca*: **False**.
- The internal pairing within var. *kuoi* is weak/moderate and should not be emphasized.
- Var. *longinux* and var. *lativiolaciifolia* are not reciprocally monophyletic (see `svdquartets_variety_monophyly.tsv`).

## Comparison with existing trees

The SVDquartets tree recovers the three var. *kuoi* populations as a moderately supported lineage, whereas populations assigned to var. *longinux* and var. *lativiolaciifolia* are interspersed. The rooted SVDquartets topology does **not** recover all non-kuoi populations as a clade, so it agrees with SNAPPER/TreeMix on the grouping of the three *kuoi* populations but not necessarily on their basal placement. SVDquartets estimates an unrooted topology; the displayed tree is rooted on the included *Q. glauca* population tip (`Q_glauca`).

## Limitations

The analysis uses sparse, ascertainment-filtered RAD SNPs. The one-representative sampling matches the SNAPPER design but discards within-population information. Bootstrap values therefore quantify repeatability under this marker set and sampling scheme, not certainty about the complete genome history. The tree should be supplied as a supplementary robustness analysis rather than replacing SNAPPER.
