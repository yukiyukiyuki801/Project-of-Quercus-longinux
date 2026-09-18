# SVDquartets

This workflow uses PAUP* 4.0a168 to analyze the merged ingroup and *Quercus
glauca* VCF produced by the Dsuite workflow. It runs the all-individual
population analysis and a one-representative-per-population analysis. Install
the Python packages listed in `requirements-python.txt`, then run:

```bash
export PAUP_BIN=/path/to/paup4a168
bash svdquartets/scripts/run_svdquartets.sh
```

To supply a merged VCF at another location:

```bash
bash svdquartets/scripts/run_svdquartets.sh /path/to/merged.vcf.gz
```

Archived trees, support tables, and figures are in `outputs/all_individuals/`
and `outputs/representatives/`.
