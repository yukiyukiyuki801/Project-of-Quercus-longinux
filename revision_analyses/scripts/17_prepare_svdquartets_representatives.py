#!/usr/bin/env python3
"""Make a one-lowest-missing-representative-per-population SVDQ matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


IUPAC = {
    frozenset(("A", "G")): "R", frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S", frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K", frozenset(("A", "C")): "M",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--outgroup-sample", default="DH01_Qglauca")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--mapping", type=Path, required=True)
    args = ap.parse_args()

    with args.metadata.open(encoding="utf-8") as handle:
        meta = {r["sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")}

    opener = gzip.open if args.vcf.suffix == ".gz" else open
    samples = []
    sequences = {}
    with opener(args.vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            f = line.rstrip().split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                sequences = {s: [] for s in samples}
                continue
            ref, alt = f[3].upper(), f[4].upper()
            gt_index = f[8].split(":").index("GT")
            for sample, value in zip(samples, f[9:]):
                gt = value.split(":")[gt_index].replace("|", "/").split("/")
                if len(gt) != 2 or "." in gt:
                    state = "?"
                else:
                    alleles = [ref if x == "0" else alt if x == "1" else "?" for x in gt]
                    if "?" in alleles:
                        state = "?"
                    elif alleles[0] == alleles[1]:
                        state = alleles[0]
                    else:
                        state = IUPAC[frozenset(alleles)]
                sequences[sample].append(state)

    by_population = defaultdict(list)
    for sample in samples:
        if sample != args.outgroup_sample:
            by_population[meta[sample]["population"]].append(sample)
    representatives = {
        pop: min(members, key=lambda s: (sequences[s].count("?"), s))
        for pop, members in by_population.items()
    }

    rows = []
    for pop in sorted(representatives):
        sample = representatives[pop]
        rows.append({
            "tip": pop, "sample_id": sample,
            "variety": meta[sample]["variety"],
            "missing_sites": sequences[sample].count("?"),
        })
    rows.append({
        "tip": "Q_glauca", "sample_id": args.outgroup_sample,
        "variety": "outgroup", "missing_sites": sequences[args.outgroup_sample].count("?"),
    })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as out:
        out.write("#NEXUS\n\nbegin data;\n")
        out.write(f"  dimensions ntax={len(rows)} nchar={len(next(iter(sequences.values())))};\n")
        out.write("  format datatype=dna missing=? gap=- interleave=no;\n  matrix\n")
        width = max(len(r["tip"]) for r in rows)
        for row in rows:
            out.write(f"  {row['tip'].ljust(width)}  {''.join(sequences[row['sample_id']])}\n")
        out.write("  ;\nend;\n")

    with args.mapping.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} tips; selected one lowest-missing individual for each of {len(representatives)} populations")


if __name__ == "__main__":
    main()
