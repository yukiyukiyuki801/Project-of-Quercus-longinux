#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
analysis_dir="${repo_dir}/revision_analyses"
assembly_report="${QROBUR_ASSEMBLY_REPORT:?Set QROBUR_ASSEMBLY_REPORT to the PM1N assembly report}"
gff="${QROBUR_GFF:?Set QROBUR_GFF to the PM1N GFF3, optionally gzip-compressed}"
proteins="${QROBUR_PROTEINS:?Set QROBUR_PROTEINS to the PM1N protein FASTA, optionally gzip-compressed}"
arabidopsis="${ARABIDOPSIS_REVIEWED_FASTA:?Set ARABIDOPSIS_REVIEWED_FASTA to reviewed Arabidopsis UniProt FASTA}"
work_dir="${analysis_dir}/work/gene_annotation"
mkdir -p "${work_dir}"

python3 "${analysis_dir}/scripts/01_annotate_outliers.py" \
  --outliers "${repo_dir}/ql_rda_fst_overlap_outliers.csv" \
  --assembly-report "${assembly_report}" \
  --gff "${gff}" \
  --proteins "${proteins}" \
  --output "${work_dir}/outlier_nearest_genes.tsv" \
  --protein-output "${work_dir}/outlier_nearest_genes.proteins.fa"

makeblastdb -in "${arabidopsis}" -dbtype prot -out "${work_dir}/arabidopsis_reviewed"
blastp -query "${work_dir}/outlier_nearest_genes.proteins.fa" \
  -db "${work_dir}/arabidopsis_reviewed" -evalue 1e-5 -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle' \
  -out "${work_dir}/outlier_protein_blastp_arabidopsis.tsv"

python3 "${analysis_dir}/scripts/select_blast_top_hits.py" \
  --blast "${work_dir}/outlier_protein_blastp_arabidopsis.tsv" \
  --output "${work_dir}/outlier_protein_top_hits.tsv"
python3 "${analysis_dir}/scripts/fetch_uniprot_annotations.py" \
  --hits "${work_dir}/outlier_protein_top_hits.tsv" \
  --output "${work_dir}/arabidopsis_top_hit_uniprot_annotations.tsv"
python3 "${analysis_dir}/scripts/run_gprofiler.py" \
  --hits "${work_dir}/outlier_protein_top_hits.tsv" \
  --output "${work_dir}/gprofiler_result.json"
python3 "${analysis_dir}/scripts/02_finalize_gene_annotation.py" \
  --proximity "${work_dir}/outlier_nearest_genes.tsv" \
  --hits "${work_dir}/outlier_protein_top_hits.tsv" \
  --uniprot "${work_dir}/arabidopsis_top_hit_uniprot_annotations.tsv" \
  --output "${work_dir}/outlier_functional_annotation.tsv" \
  --summary "${work_dir}/README.md"
