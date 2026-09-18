#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
analysis_dir="${repo_dir}/svdquartets"
paup_bin="${PAUP_BIN:?Set PAUP_BIN to the PAUP* 4.0a168 executable}"
merged_vcf="${1:-${repo_dir}/dsuite/work/ql_180_plus_qglauca.vcf.gz}"
work_dir="${analysis_dir}/work/all_individuals"
mkdir -p "${work_dir}"

python3 "${analysis_dir}/scripts/05_prepare_svdquartets.py" \
  --vcf "${merged_vcf}" \
  --metadata "${repo_dir}/data/sample_metadata.tsv" \
  --outgroup-sample DH01_Qglauca \
  --outgroup-population Q_glauca \
  --output "${work_dir}/matrix.nex"

cat > "${work_dir}/run.nex" <<EOF
#NEXUS
begin paup;
  execute ${work_dir}/matrix.nex;
  log file=${work_dir}/svdquartets.log replace;
  svdquartets taxpartition=populations evalquartets=random nquartets=100000 showScores=no bootstrap nreps=100 seed=20260917 treeFile=${work_dir}/bootstrap.tre;
  savetrees file=${work_dir}/consensus.tre format=newick brlens=no replace;
  log stop;
  quit;
end;
EOF

"${paup_bin}" -n -u "${work_dir}/run.nex"
python3 "${analysis_dir}/scripts/10_summarize_svdquartets.py" \
  --tree "${work_dir}/consensus.tre" \
  --metadata "${repo_dir}/data/sample_metadata.tsv" \
  --outgroup Q_glauca \
  --sampling all_individuals \
  --output-dir "${work_dir}/summary"

# Sampling sensitivity matching the one-representative-per-population design.
representative_dir="${analysis_dir}/work/representatives"
mkdir -p "${representative_dir}"
python3 "${analysis_dir}/scripts/17_prepare_svdquartets_representatives.py" \
  --vcf "${merged_vcf}" \
  --metadata "${repo_dir}/data/sample_metadata.tsv" \
  --outgroup-sample DH01_Qglauca \
  --output "${representative_dir}/matrix.nex" \
  --mapping "${representative_dir}/selected_samples.tsv"
cat > "${representative_dir}/run.nex" <<EOF
#NEXUS
begin paup;
  execute ${representative_dir}/matrix.nex;
  log file=${representative_dir}/svdquartets.log replace;
  svdquartets evalquartets=all showScores=no bootstrap nreps=100 seed=20260917 treeFile=${representative_dir}/bootstrap.tre;
  savetrees file=${representative_dir}/consensus.tre format=newick brlens=no replace;
  log stop;
  quit;
end;
EOF
"${paup_bin}" -n -u "${representative_dir}/run.nex"
python3 "${analysis_dir}/scripts/10_summarize_svdquartets.py" \
  --tree "${representative_dir}/consensus.tre" \
  --metadata "${repo_dir}/data/sample_metadata.tsv" \
  --outgroup Q_glauca \
  --sampling representatives \
  --output-dir "${representative_dir}/summary"
