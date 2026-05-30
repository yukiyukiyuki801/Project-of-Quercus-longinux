# Gradient forest current/LGM projection

This directory contains the code and input files needed to fit the gradient forest model and generate `gf_current.pdf` for the current and LGM genomic composition projection in:

**Postglacial expansion and sieving of adaptive variation contribute to speciation in a subtropical oak complex**

The color palette is assigned from a shared PCA of current and LGM GF predictions:

```r
R = PC1
G = PC2
B = PC3
```

Because current and LGM predictions are colored in the same PCA space, similar colors indicate similar GF-predicted genomic composition.

## Run

```bash
Rscript run_gf_model_projection.R
```

Required R packages: `raster`, `rgdal`, `vegan`, and `gradientForest`.

The script writes only:

```text
gf_pca_combine.pdf
gf_current.pdf
```
