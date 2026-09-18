#!/usr/bin/env bash
#SBATCH --job-name=ql_structure
#SBATCH --array=0-69%35
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/structure_%A_%a.out
#SBATCH --error=logs/structure_%A_%a.err

set -euo pipefail

repo_dir="${SLURM_SUBMIT_DIR}"
analysis_dir="${repo_dir}/structure"
structure_bin="${STRUCTURE_BIN:?Set STRUCTURE_BIN to the STRUCTURE 2.3.4 executable}"
task="${SLURM_ARRAY_TASK_ID}"
k=$((task / 10 + 1))
rep=$((task % 10 + 1))
seed=$((2026091700 + task))

mkdir -p "${repo_dir}/logs" "${analysis_dir}/work/raw/K${k}"
"${structure_bin}" \
  -m "${analysis_dir}/config/mainparams" \
  -e "${analysis_dir}/config/extraparams" \
  -i "${analysis_dir}/work/structure_input.tsv" \
  -o "${analysis_dir}/work/raw/K${k}/K${k}_rep${rep}" \
  -K "${k}" -D "${seed}"
