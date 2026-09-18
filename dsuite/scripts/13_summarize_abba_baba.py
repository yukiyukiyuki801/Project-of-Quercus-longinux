#!/usr/bin/env python3
"""Create a concise numeric summary of the variety-level Dsuite result."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree-result", type=Path, required=True)
    ap.add_argument("--qc", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    qc = json.loads(args.qc.read_text(encoding="utf-8"))
    with args.tree_result.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit(f"Expected one topology-constrained trio, found {len(rows)}")
    r = rows[0]
    d = float(r["Dstatistic"])
    z = float(r["Z-score"])
    p = float(r["p-value"])
    f4 = float(r["f4-ratio"])
    abba = float(r["ABBA"])
    baba = float(r["BABA"])
    sensitivity_lines = ["| Minimum DH01 depth | Usable loci | D | Z | P |", "|---:|---:|---:|---:|---:|"]
    sensitivity_lines.append(
        f"| {qc['minimum_depth']} | {qc.get('usable', 0)} | {d:.4f} | {z:.3f} | {p:.4g} |"
    )
    sensitivity_dirs = sorted(
        args.tree_result.parent.glob("sensitivity_dp*"),
        key=lambda path: int(path.name.removeprefix("sensitivity_dp")),
    )
    for directory in sensitivity_dirs:
        result_path = directory / "trio_tree.txt"
        qc_path = directory / "qc.json"
        if not result_path.exists() or not qc_path.exists():
            continue
        with result_path.open(encoding="utf-8") as handle:
            sr = next(csv.DictReader(handle, delimiter="\t"))
        sq = json.loads(qc_path.read_text(encoding="utf-8"))
        sensitivity_lines.append(
            f"| {sq['minimum_depth']} | {sq.get('usable', 0)} | {float(sr['Dstatistic']):.4f} | "
            f"{float(sr['Z-score']):.3f} | {float(sr['p-value']):.4g} |"
        )

    text = f"""# ABBA-BABA / D-statistic sensitivity analysis

## Result

- Topology tested: `(({r['P1']},{r['P2']}),{r['P3']})`, rooted with *Quercus glauca* DH01.
- Patterson's D = **{d:.4f}**; block-jackknife Z = **{z:.3f}**; two-sided P = **{p:.4g}**.
- f4-ratio = **{f4:.4f}**.
- Weighted site-pattern counts: ABBA = {abba:.2f}; BABA = {baba:.2f}.
- Q. glauca genotypes passing DP >= {qc['minimum_depth']}: **{qc.get('usable', 0)} / {qc['ingroup_sites']}** ({100*qc['usable_fraction']:.1f}%).

## Outgroup-depth sensitivity

{chr(10).join(sensitivity_lines)}

## Reproducibility

- Outgroup reads: public ddRAD-seq sample *Q. glauca* DH01 (SRS19813215; run SRR27151155; SbfI/MspI).
- Mapping reference: *Q. robur* PM1N, the same assembly coordinates used by the ingroup VCF.
- Mapping/calling: BWA-MEM followed by samtools/bcftools, mapping and base quality >= 20; outgroup DP >= 3.
- Statistic: Dsuite Dtrios v0.5 r58, expected variety topology, 20 block-jackknife blocks.
"""
    args.output.write_text(text, encoding="utf-8")
    print(f"D={d:.4f}; Z={z:.3f}; P={p:.4g}")


if __name__ == "__main__":
    main()
