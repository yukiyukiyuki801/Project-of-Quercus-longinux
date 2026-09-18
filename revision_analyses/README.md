# Reviewer-requested and revision analyses

This directory contains code and compact outputs added during manuscript
revision. The published input VCF and sample metadata are in `../data/`.
Large or readily downloadable intermediates—raw reads, BAM files, reference
genomes, BLAST databases, and raw STRUCTURE runs—are intentionally excluded.

## Included analyses

| Analysis | Run/preparation code | Summary code | Archived output |
|---|---|---|---|
| Pairwise diversity tests | `scripts/verify_diversity_tests.py` | same script | `outputs/rank_tests/` |
| STRUCTURE, K = 1–7 | `03_prepare_structure.py`, `run_structure_array.slurm` | `09_summarize_structure.py` | `outputs/structure/` |
| SVDquartets | `05_prepare_svdquartets.py`, `17_prepare_svdquartets_representatives.py`, `run_svdquartets.sh` | `10_summarize_svdquartets.py` | `outputs/svdquartets/`, `outputs/svdquartets_representatives/` |
| Patterson's D / Dsuite | `07_prepare_outgroup_sites.py`, `map_qglauca_outgroup.slurm`, `11_merge_qglauca_outgroup.py`, `run_dsuite_depths.sh` | `13_summarize_abba_baba.py` | `outputs/dsuite/` |
| Outlier functional annotation | `01_annotate_outliers.py`, `run_gene_annotation.sh` and helper scripts | `02_finalize_gene_annotation.py`, `run_gprofiler.py` | `outputs/gene_annotation/` |

Create a Python environment with:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r revision_analyses/requirements.txt
```

All commands below are run from the repository root.

## Independent-group diversity tests

Population estimates for He, Ho, FIS, and nucleotide diversity are supplied
in `data/diversity_by_population.tsv`. The script performs the three pairwise,
two-sided Mann–Whitney/Wilcoxon rank-sum tests within each metric and applies
Holm correction to those three P values.

```bash
python3 revision_analyses/scripts/verify_diversity_tests.py \
  --input revision_analyses/data/diversity_by_population.tsv \
  --output revision_analyses/work/rank_tests/pairwise_wilcoxon_rank_sum.tsv
```

## STRUCTURE and Evanno Delta K

Prepare the two-row-per-diploid-individual input:

```bash
mkdir -p revision_analyses/work/structure
python3 revision_analyses/scripts/03_prepare_structure.py \
  --vcf data/ql_180_LD_prune_0.2.recode.vcf \
  --metadata data/sample_metadata.tsv \
  --output revision_analyses/work/structure/structure_input.tsv \
  --population-key revision_analyses/work/structure/population_key.tsv
```

Run 70 jobs (10 replicates for each K = 1–7) with STRUCTURE 2.3.4:

```bash
export STRUCTURE_BIN=/path/to/structure
sbatch --export=ALL,STRUCTURE_BIN revision_analyses/scripts/run_structure_array.slurm
```

After all jobs finish, align replicate cluster labels, calculate Evanno Delta
K from replicate likelihoods, average Q matrices, and plot K = 2–7:

```bash
python3 revision_analyses/scripts/09_summarize_structure.py \
  --raw-dir revision_analyses/work/structure/raw \
  --metadata data/sample_metadata.tsv \
  --output-dir revision_analyses/work/structure/summary
```

The archived run has maximum Delta K at K = 3 (174.613).

## Q. glauca outgroup and Dsuite

The outgroup is the public ddRAD sample DH01, SRA run SRR27151155. Download the
Q. robur PM1N FASTA and its assembly report before running. The included BED
and coordinate map restrict calling to the same 1,291 loci as the ingroup VCF.
The mapping job requires SRA Toolkit, BWA 0.7.17, SAMtools, and BCFtools 1.21.

```bash
export QROBUR_FASTA=/path/to/Qrob_PM1N.fa
sbatch --export=ALL,QROBUR_FASTA revision_analyses/scripts/map_qglauca_outgroup.slurm
```

Run Dsuite 0.5 r58 at outgroup minimum depths 3, 5, and 10. The script creates
the variety assignment file and topology, runs `Dtrios` with 20 jackknife
blocks, and writes a guarded summary.

```bash
export DSUITE_BIN=/path/to/Dsuite
bash revision_analyses/scripts/run_dsuite_depths.sh
```

The primary result is D = 0.0147, Z = 1.012, P = 0.311; all depth sensitivity
tests are nonsignificant. This is a sparse, dataset-wide RAD sensitivity test,
not a windowed genome scan or proof that historical gene flow was absent.

## SVDquartets

SVDquartets uses the merged ingroup plus DH01 VCF produced by the Dsuite
workflow. PAUP* 4.0a168 evaluates 100,000 randomly sampled lineage quartets
with 100 standard bootstrap replicates and seed 20260917. The summarizer roots
the unrooted consensus post hoc on the Q. glauca population partition and
plots every internal support value.

```bash
export PAUP_BIN=/path/to/paup4a168
bash revision_analyses/scripts/run_svdquartets.sh
```

The three kuoi populations form a clade with 100% bootstrap support. The
outgroup-adjacent split has 27% support, so the basal ingroup relationship is
not resolved by this analysis. The same wrapper also selects the
lowest-missing-data individual from each population, evaluates all 17,550
possible quartets, and summarizes this sampling sensitivity analysis; it
recovers the kuoi clade with 62% bootstrap support.

## Functional annotation and enrichment

Download the Q. robur PM1N assembly report, GFF3, and protein FASTA and the
reviewed Arabidopsis thaliana UniProt protein FASTA. The workflow translates
VCF accessions to PM1N sequence names, classifies overlap/proximity to gene
models, extracts proteins, searches reviewed Arabidopsis proteins with BLASTP,
retrieves UniProt annotations, and queries g:Profiler for GO and KEGG terms.

```bash
export QROBUR_ASSEMBLY_REPORT=/path/to/GCA_900291515.1_assembly_report.txt
export QROBUR_GFF=/path/to/Qrob_PM1N_genes.gff.gz
export QROBUR_PROTEINS=/path/to/Qrob_PM1N_CDS_aa.fa.gz
export ARABIDOPSIS_REVIEWED_FASTA=/path/to/arabidopsis_reviewed_uniprot.fasta
bash revision_analyses/scripts/run_gene_annotation.sh
```

The archived table contains all 41 shared FST–RDA candidates. Twenty-four are
coding and five are intronic/transcribed. Twenty-one of 26 assigned proteins
have reviewed Arabidopsis homologues. The only corrected enrichment result is
the small ABC-transporter KEGG term (two proteins; adjusted P = 0.0241). These
are proximity- and homology-based hypotheses, not causal annotations.

## Archived-output integrity

`outputs/SHA256SUMS` records checksums for all compact result files deposited
with the revision. Regenerate it from this directory with:

```bash
find outputs -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > outputs/SHA256SUMS
```
