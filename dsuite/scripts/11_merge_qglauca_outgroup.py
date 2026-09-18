#!/usr/bin/env python3
"""Append the Q. glauca genotype at matched PM1N loci to the ingroup VCF.

The ingroup VCF uses INSDC accession names, whereas the PM1N FASTA and the
Q. glauca call use scaffold/chromosome names.  The coordinate-map file makes
that one-to-one naming translation without changing positions.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_sample(fmt: str, sample: str) -> dict[str, str]:
    return dict(zip(fmt.split(":"), sample.split(":")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingroup-vcf", type=Path, required=True)
    ap.add_argument("--outgroup-vcf", type=Path, required=True)
    ap.add_argument("--coordinate-map", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--output-vcf", type=Path, required=True)
    ap.add_argument("--sets", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--qc", type=Path, required=True)
    ap.add_argument("--minimum-depth", type=int, default=3)
    args = ap.parse_args()

    accession_to_pm1n = {}
    with args.coordinate_map.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession_to_pm1n[(row["accession"], int(row["position"]))] = (
                row["pm1n_sequence"], int(row["position"])
            )

    calls = {}
    with open_text(args.outgroup_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            data = parse_sample(f[8], f[9])
            calls[(f[0], int(f[1]))] = {
                "alleles": [f[3], *f[4].split(",")],
                "gt": data.get("GT", "./."),
                "dp": data.get("DP", "."),
                "ad": data.get("AD", "."),
            }

    metadata = {}
    with args.metadata.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata[row["sample_id"]] = row["variety"]

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    samples = []
    with args.ingroup_vcf.open(encoding="utf-8") as src, args.output_vcf.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("##"):
                dst.write(line)
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                samples = fields[9:]
                dst.write(line.rstrip() + "\tDH01_Qglauca\n")
                continue

            f = line.rstrip().split("\t")
            counts["ingroup_sites"] += 1
            pm1n_key = accession_to_pm1n.get((f[0], int(f[1])))
            call = calls.get(pm1n_key) if pm1n_key else None
            target_alleles = [f[3], f[4]]
            out_gt, out_dp, out_ad = "./.", ".", ".,."
            reason = "call_missing"
            if not pm1n_key:
                reason = "coordinate_missing"
            elif call:
                if call["alleles"][0] == target_alleles[1]:
                    counts["pm1n_reference_equals_ingroup_alt"] += 1
                elif call["alleles"][0] != target_alleles[0]:
                    counts["pm1n_reference_is_third_allele"] += 1
                try:
                    dp = int(call["dp"])
                except (TypeError, ValueError):
                    dp = -1
                gt_tokens = call["gt"].replace("|", "/").split("/")
                called_bases = []
                try:
                    called_bases = [call["alleles"][int(x)] for x in gt_tokens if x != "."]
                except (ValueError, IndexError):
                    called_bases = []
                if dp < args.minimum_depth:
                    reason = "depth_below_threshold"
                elif len(called_bases) != 2:
                    reason = "genotype_missing"
                elif any(base not in target_alleles for base in called_bases):
                    reason = "non_target_allele"
                else:
                    translated = [str(target_alleles.index(base)) for base in called_bases]
                    out_gt = "/".join(translated)
                    out_dp = str(dp)
                    if call["ad"] not in {None, "."}:
                        source_ad = call["ad"].split(",")
                        ad_by_base = dict(zip(call["alleles"], source_ad))
                        out_ad = ",".join(ad_by_base.get(base, "0") for base in target_alleles)
                    reason = "usable"
            counts[reason] += 1
            # The ingroup has GT:PL:DP:AD:GP:GQ at every retained site.
            dst.write(line.rstrip() + f"\t{out_gt}:.:{out_dp}:{out_ad}:.:.\n")

    missing_metadata = sorted(set(samples) - set(metadata))
    if missing_metadata:
        raise SystemExit(f"Samples absent from metadata: {missing_metadata}")
    args.sets.parent.mkdir(parents=True, exist_ok=True)
    with args.sets.open("w", encoding="utf-8") as handle:
        for sample in samples:
            handle.write(f"{sample}\t{metadata[sample]}\n")
        handle.write("DH01_Qglauca\tOutgroup\n")
    args.tree.write_text("((longinux,lativiolaciifolia),kuoi);\n", encoding="utf-8")

    report = {
        "minimum_depth": args.minimum_depth,
        "outgroup_vcf_records": len(calls),
        **counts,
        "usable_fraction": counts["usable"] / counts["ingroup_sites"] if counts["ingroup_sites"] else 0,
    }
    args.qc.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
