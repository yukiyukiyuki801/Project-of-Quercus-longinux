# Reproducible analyses for the *Quercus longinux* project

Code and available inputs for the manuscript **Postglacial expansion and
sieving of adaptive variation contribute to speciation in a subtropical oak
complex**.

This repository supports:

1. Genetic PCA and visualization.
2. Gradient Forest projections using adaptive or genome-wide SNPs, plus a
   comparison of their predictions.
3. Population directionality (ψ) and range-origin inference with the
   `rangeExpansion` package.
4. Tajima's D from the filtered LD-pruned RAD VCF.

The unfiltered VCF is too large to include. Please contact the authors if it is
needed.

## Contents

```text
180_ql.ld_prune.lfmm.csv          genotype matrix
ql_180_env_variables.csv          sample coordinates and predictors
ql_rda_fst_overlap_outliers.csv   candidate loci
current_bio19/                    current climate rasters
lgm_bio19/                        Last Glacial Maximum rasters
mask_layer/                       habitat and background masks
data/
  ql_180_LD_prune_0.2.recode.vcf  filtered LD-pruned VCF
  sample_metadata.tsv             sample, population, variety, and coordinates
scripts/
  run_genetic_pca.R
  gf_workflow.R
  run_gradient_forest.R
  compare_gradient_forest.R
  run_range_expansion.R
  prepare_tajima_sample_lists.py
  run_tajimas_d.sh
  summarize_tajimas_d.R
run_gf_model_projection.R         compatibility entry point
```

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

## Genetic PCA

Run PCA on the LD-pruned genotype matrix and plot the first two axes by variety:

```bash
Rscript scripts/run_genetic_pca.R
```

## Gradient Forest

The models use BIO01, BIO07, BIO12, BIO17, and BIO18; 500 trees by default;
`corr.threshold = 0.50`; and the manuscript's maximum-level calculation.

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


## Exploratory Tajima's D

```bash
python3 scripts/prepare_tajima_sample_lists.py
bash scripts/run_tajimas_d.sh
Rscript scripts/summarize_tajimas_d.R
```

The default window is 10 kb. It can be changed with, for example,
`WINDOW=20000 bash scripts/run_tajimas_d.sh`.
