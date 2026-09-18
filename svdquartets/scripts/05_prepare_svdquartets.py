#!/usr/bin/env python3
"""Convert the biallelic VCF to a PAUP* NEXUS matrix with population partition."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--outgroup-sample")
    parser.add_argument("--outgroup-population", default="Q_glauca")
    args = parser.parse_args()

    meta_rows = list(csv.DictReader(args.metadata.open(encoding="utf-8"), delimiter="\t"))
    meta = {x["sample_id"]: x for x in meta_rows}
    populations: dict[str, list[str]] = defaultdict(list)
    samples: list[str] = []
    sequences: dict[str, list[str]] = {}
    nchar = 0

    opener = gzip.open if args.vcf.suffix == ".gz" else open
    with opener(args.vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip().split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                sequences = {sample: [] for sample in samples}
                for sample in samples:
                    if sample == args.outgroup_sample:
                        populations[args.outgroup_population].append(sample)
                    elif sample not in meta:
                        raise SystemExit(f"Missing metadata for {sample}")
                    else:
                        populations[meta[sample]["population"]].append(sample)
                continue
            ref, alt = fields[3].upper(), fields[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                raise SystemExit(f"Non-SNP/non-biallelic record at {fields[0]}:{fields[1]}")
            gt_index = fields[8].split(":").index("GT")
            for sample, value in zip(samples, fields[9:]):
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
            nchar += 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as out:
        out.write("#NEXUS\n\n")
        out.write("begin data;\n")
        out.write(f"  dimensions ntax={len(samples)} nchar={nchar};\n")
        out.write("  format datatype=dna missing=? gap=- interleave=no;\n")
        out.write("  matrix\n")
        width = max(map(len, samples))
        for sample in samples:
            out.write(f"  {sample.ljust(width)}  {''.join(sequences[sample])}\n")
        out.write("  ;\nend;\n\n")
        out.write("begin sets;\n")
        out.write("  taxpartition populations =\n")
        parts = []
        for population in sorted(populations):
            parts.append(f"    {population}: " + " ".join(populations[population]))
        out.write(",\n".join(parts) + ";\n")
        out.write("end;\n")
    print(f"Wrote {len(samples)} samples, {nchar} SNPs, and {len(populations)} populations")


if __name__ == "__main__":
    main()
