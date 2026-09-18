#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
analysis_dir="${repo_dir}/dsuite"
dsuite_bin="${DSUITE_BIN:?Set DSUITE_BIN to the Dsuite executable}"
outgroup_vcf="${1:-${analysis_dir}/work/outgroup/DH01.known_sites.vcf.gz}"
work_root="${analysis_dir}/work"
mkdir -p "${work_root}"

for depth in 3 5 10; do
  if [[ "${depth}" == 3 ]]; then
    run_dir="${work_root}"
  else
    run_dir="${work_root}/sensitivity_dp${depth}"
  fi
  mkdir -p "${run_dir}"
  python3 "${analysis_dir}/scripts/11_merge_qglauca_outgroup.py" \
    --ingroup-vcf "${repo_dir}/data/ql_180_LD_prune_0.2.recode.vcf" \
    --outgroup-vcf "${outgroup_vcf}" \
    --coordinate-map "${analysis_dir}/config/known_sites.coordinate_map.tsv" \
    --metadata "${repo_dir}/data/sample_metadata.tsv" \
    --output-vcf "${run_dir}/ql_180_plus_qglauca.vcf" \
    --sets "${run_dir}/sets.tsv" \
    --tree "${run_dir}/tree.nwk" \
    --qc "${run_dir}/qc.json" \
    --minimum-depth "${depth}"
  bgzip -f "${run_dir}/ql_180_plus_qglauca.vcf"
  tabix -f -p vcf "${run_dir}/ql_180_plus_qglauca.vcf.gz"
  "${dsuite_bin}" Dtrios -k 20 -t "${run_dir}/tree.nwk" -o "${run_dir}/trio" \
    "${run_dir}/ql_180_plus_qglauca.vcf.gz" "${run_dir}/sets.tsv"
done

python3 "${analysis_dir}/scripts/13_summarize_abba_baba.py" \
  --tree-result "${work_root}/trio_tree.txt" \
  --qc "${work_root}/qc.json" \
  --output "${work_root}/D_STATISTIC_SUMMARY.md"
