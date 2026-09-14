#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd -- "${script_dir}/.." && pwd)"
vcf="${VCF:-${project_root}/data/ql_180_LD_prune_0.2.recode.vcf}"
sample_dir="${project_root}/results/tajimas_d/sample_lists"
output_dir="${project_root}/results/tajimas_d/windows"
window="${WINDOW:-10000}"

if ! command -v vcftools >/dev/null 2>&1; then
  echo "ERROR: vcftools is not available on PATH" >&2
  exit 1
fi
if [[ ! -s "${vcf}" ]]; then
  echo "ERROR: filtered LD-pruned VCF not found: ${vcf}" >&2
  exit 1
fi
if [[ ! "${window}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: WINDOW must be a positive integer" >&2
  exit 1
fi

python3 "${script_dir}/prepare_tajima_sample_lists.py"
mkdir -p "${output_dir}"
shopt -s nullglob
sample_lists=("${sample_dir}"/*.samples.txt)
if (( ${#sample_lists[@]} == 0 )); then
  echo "ERROR: no sample lists found under ${sample_dir}" >&2
  exit 1
fi

for keep_file in "${sample_lists[@]}"; do
  label="$(basename "${keep_file}" .samples.txt)"
  prefix="${output_dir}/${label}_${window}bp_SNP_ONLY_EXPLORATORY"
  vcftools --vcf "${vcf}" --keep "${keep_file}" --TajimaD "${window}" --out "${prefix}"
done

printf '%s\n' \
  "These Tajima's D values use a filtered, LD-pruned, variant-only RAD VCF." \
  "Invariant callable sites were unavailable. Treat values only as exploratory" \
  "relative comparisons, not absolute genome-wide demographic estimates." \
  > "${output_dir}/README.txt"
echo "Wrote exploratory Tajima's D windows to ${output_dir}"
