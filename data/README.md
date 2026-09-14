# Data included in this directory

## `ql_180_LD_prune_0.2.recode.vcf`

Filtered, LD-pruned, variant-only RAD-seq VCF containing 1,291 SNPs and 180
individuals. This is not the unfiltered VCF. Invariant callable sites are absent,
so statistics such as Tajima's D must be interpreted only as exploratory
relative comparisons.

## `sample_metadata.tsv`

Tab-delimited metadata for the same 180 VCF samples. Columns are:

- `sample_id`: exact VCF sample identifier
- `population`: population code
- `variety`: population-level taxonomic assignment
- `analysis_group`: `kuoi` or combined `rest`
- `longitude`, `latitude`: sampling coordinates in decimal degrees

Every VCF sample occurs exactly once in this metadata table.
