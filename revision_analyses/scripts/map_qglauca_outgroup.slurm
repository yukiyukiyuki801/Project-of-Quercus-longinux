#!/usr/bin/env bash
#SBATCH --job-name=ql_qglauca
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/qglauca_%j.out
#SBATCH --error=logs/qglauca_%j.err

set -euo pipefail

repo_dir="${SLURM_SUBMIT_DIR}"
analysis_dir="${repo_dir}/revision_analyses"
ref="${QROBUR_FASTA:?Set QROBUR_FASTA to the uncompressed PM1N FASTA}"
outdir="${analysis_dir}/work/dsuite/outgroup"
mkdir -p "${repo_dir}/logs" "${outdir}/fastq" "${outdir}/tmp"

fasterq-dump SRR27151155 --split-files --threads "${SLURM_CPUS_PER_TASK}" \
  --temp "${outdir}/tmp" --outdir "${outdir}/fastq"
bwa index "${ref}"
samtools faidx "${ref}"
bwa mem -t "${SLURM_CPUS_PER_TASK}" -R '@RG\tID:DH01\tSM:DH01\tPL:ILLUMINA' \
  "${ref}" "${outdir}/fastq/SRR27151155_1.fastq" "${outdir}/fastq/SRR27151155_2.fastq" |
  samtools sort -@ 4 -o "${outdir}/DH01.Qrobur_PM1N.bam" -
samtools index "${outdir}/DH01.Qrobur_PM1N.bam"

bcftools mpileup -f "${ref}" -R "${analysis_dir}/config/dsuite/known_sites.pm1n.bed" \
  -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -Ou "${outdir}/DH01.Qrobur_PM1N.bam" |
  bcftools call -Am -Ou |
  bcftools norm -f "${ref}" -Oz -o "${outdir}/DH01.known_sites.vcf.gz"
bcftools index -f "${outdir}/DH01.known_sites.vcf.gz"
