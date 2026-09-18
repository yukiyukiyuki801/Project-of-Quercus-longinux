# Gene annotation and enrichment

This workflow maps shared FST-RDA candidate SNPs to *Quercus robur* PM1N gene
models, searches reviewed *Arabidopsis thaliana* proteins with BLASTP, obtains
UniProt annotations, and queries g:Profiler for GO and KEGG enrichment.

Install the Python dependencies from `requirements-python.txt`, then run:

```bash
export QROBUR_ASSEMBLY_REPORT=/path/to/GCA_900291515.1_assembly_report.txt
export QROBUR_GFF=/path/to/Qrob_PM1N_genes.gff.gz
export QROBUR_PROTEINS=/path/to/Qrob_PM1N_CDS_aa.fa.gz
export ARABIDOPSIS_REVIEWED_FASTA=/path/to/arabidopsis_reviewed_uniprot.fasta
bash gene_annotation/scripts/run_gene_annotation.sh
```

The workflow also requires BLAST+. Archived annotation and enrichment files
are in `outputs/`.
