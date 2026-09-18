# STRUCTURE and Evanno Delta K

Prepare the STRUCTURE input from the LD-pruned VCF:

```bash
mkdir -p structure/work
python3 structure/scripts/03_prepare_structure.py \
  --vcf data/ql_180_LD_prune_0.2.recode.vcf \
  --metadata data/sample_metadata.tsv \
  --output structure/work/structure_input.tsv \
  --population-key structure/work/population_key.tsv
```

Run 10 replicates for each K from 1 to 7 with STRUCTURE 2.3.4:

```bash
export STRUCTURE_BIN=/path/to/structure
sbatch --export=ALL,STRUCTURE_BIN structure/scripts/run_structure_array.slurm
```

Align cluster labels, calculate Evanno Delta K, average Q matrices, and create
the summary plots:

```bash
python3 structure/scripts/09_summarize_structure.py \
  --raw-dir structure/work/raw \
  --metadata data/sample_metadata.tsv \
  --output-dir structure/work/summary
```

Archived tables and figures are in `outputs/`.
