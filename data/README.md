# Data included in this directory

## `ql_180_LD_prune_0.2.recode.vcf`

Filtered, LD-pruned RAD-seq VCF used by the included scripts. The unfiltered VCF
is too large to include; please contact the authors if it is needed.

## `sample_metadata.tsv`

Tab-delimited metadata for the VCF samples. Columns are:

- `sample_id`: exact VCF sample identifier
- `population`: population code
- `variety`: population-level taxonomic assignment
- `analysis_group`: `kuoi` or combined `rest`
- `longitude`, `latitude`: sampling coordinates in decimal degrees

Every VCF sample occurs exactly once in this metadata table.
