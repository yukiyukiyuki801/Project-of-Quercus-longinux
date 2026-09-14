# Reproducible analyses for the *Quercus longinux* project

Code and available inputs for the manuscript **Postglacial expansion and
sieving of adaptive variation contribute to speciation in a subtropical oak
complex**.

This curated repository supports:

1. Gradient Forest projections using the adaptive candidate overlap or all
   1,291 LD-pruned SNPs, plus a comparison of their predictions.
2. Population directionality (ψ) and range-origin inference with the
   `rangeExpansion` package.
3. Exploratory Tajima's D from the filtered, LD-pruned, variant-only RAD VCF.

The unfiltered VCF is intentionally not included. Incomplete working fragments
for other analyses were also excluded rather than presented as reproducible
code.

## Contents

```text
180_ql.ld_prune.lfmm.csv          180 × 1,291 genotype matrix
ql_180_env_variables.csv          sample coordinates and predictors
ql_rda_fst_overlap_outliers.csv   41 candidate loci
current_bio19/                    current climate rasters
lgm_bio19/                        Last Glacial Maximum rasters
mask_layer/                       habitat and background masks
data/
  ql_180_LD_prune_0.2.recode.vcf  filtered LD-pruned VCF
  sample_metadata.tsv             sample, population, variety, and coordinates
scripts/
  gf_workflow.R
  run_gradient_forest.R
  compare_gradient_forest.R
  run_range_expansion.R
  prepare_tajima_sample_lists.py
  run_tajimas_d.sh
  summarize_tajimas_d.R
run_gf_model_projection.R         compatibility entry point
```

Generated files are written below `results/` and ignored by Git.

## Software

- R 4.2 or a compatible version
- R packages listed in `requirements-r.txt`
- Python 3.9 or newer
- VCFtools 0.1.17 for Tajima's D

`rgdal` has been retired from CRAN but is retained because it was used for the
manuscript analysis. A frozen R environment or container is recommended for
long-term reproduction.

Install `rangeExpansion` at the upstream revision used to construct this
workflow:

```r
remotes::install_github("BenjaminPeter/rangeExpansion@d485202")
```

The Quercus range-expansion wrapper was adapted from the corresponding workflow
in the author's
[`Phalaenopsis-parallel-adaptation`](https://github.com/yukiyukiyuki801/Phalaenopsis-parallel-adaptation)
repository and from the upstream `rangeExpansion` example.

## Gradient Forest

The models use BIO01, BIO07, BIO12, BIO17, and BIO18; 500 trees by default;
`corr.threshold = 0.50`; and the manuscript's maximum-level calculation.

The candidate table has 41 loci. Twenty overlap the 1,291-SNP LD-pruned matrix,
so the matched adaptive model contains 20 loci. The other 21 candidates were
removed by LD pruning.

```bash
Rscript scripts/run_gradient_forest.R adaptive
Rscript scripts/run_gradient_forest.R all
Rscript scripts/compare_gradient_forest.R
```

The all-SNP model is computationally expensive (approximately 64 GB and eight
hours in the original cluster run). For a quick test:

```bash
GF_NTREE=5 Rscript scripts/run_gradient_forest.R adaptive
```

The original entry point remains available and invokes the cleaned adaptive
workflow:

```bash
Rscript run_gf_model_projection.R
```

## Range expansion

The range-expansion script converts the LD-pruned genotype matrix to the input
format required by `rangeExpansion`, groups individuals into the 26 named
populations, calculates pairwise ψ, and runs the TDOA origin analysis across
all populations.

```bash
Rscript scripts/run_range_expansion.R
```

The inferred origin is sensitive to SNP ascertainment, allele polarization,
population sampling, and the old package implementation. The provided matrix
has no outgroup for ancestral-state polarization; results therefore reproduce
the manuscript workflow as an exploratory directionality analysis and should
be interpreted with those limitations.

## Exploratory Tajima's D

The included VCF contains 1,291 filtered and LD-pruned variant sites for 180
individuals. It lacks invariant callable sites. Its Tajima's D values are not
absolute genome-wide estimates and must not be used as strong demographic
evidence; they are only exploratory relative comparisons.

The genomic population assignments contain 20 kuoi, 103 longinux, and 57
lativiolaciifolia individuals. The manuscript's 20/88/72 counts refer to the
morphological dataset and should not be substituted for these VCF sample
assignments.

```bash
python3 scripts/prepare_tajima_sample_lists.py
bash scripts/run_tajimas_d.sh
Rscript scripts/summarize_tajimas_d.R
```

The default window is 10 kb. It can be changed with, for example,
`WINDOW=20000 bash scripts/run_tajimas_d.sh`.

## Scope

This is a reproducibility package, not an archive of every exploratory command.
STRUCTURE, SNAPPER, TreeMix, ecological niche modeling, and fastsimcoal2 are not
included because complete manuscript-specific code and inputs were unavailable.
The raw RAD reads remain in the public archive cited by the manuscript, and the
unfiltered VCF is not distributed here.

Exact sampling coordinates are present in `ql_180_env_variables.csv` and
`data/sample_metadata.tsv`, consistent with the study's sampling information.
