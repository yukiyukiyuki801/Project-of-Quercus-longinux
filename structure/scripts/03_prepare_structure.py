#!/usr/bin/env python3
"""Convert the LD-pruned VCF to diploid STRUCTURE format."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--population-key", type=Path, required=True)
    args = parser.parse_args()

    metadata_rows = list(csv.DictReader(args.metadata.open(encoding="utf-8"), delimiter="\t"))
    metadata = {x["sample_id"]: x for x in metadata_rows}
    populations = sorted({x["population"] for x in metadata_rows})
    pop_number = {name: i + 1 for i, name in enumerate(populations)}

    samples: list[str] = []
    markers: list[str] = []
    genotypes: dict[str, list[tuple[str, str]]] = {}
    with args.vcf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip().split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                missing = [x for x in samples if x not in metadata]
                if missing:
                    raise SystemExit(f"Samples absent from metadata: {','.join(missing)}")
                genotypes = {x: [] for x in samples}
                continue
            marker = f"{fields[0]}_{fields[1]}"
            markers.append(marker)
            fmt = fields[8].split(":")
            gt_index = fmt.index("GT")
            for sample, value in zip(samples, fields[9:]):
                gt = value.split(":")[gt_index].replace("|", "/").split("/")
                if len(gt) != 2 or "." in gt:
                    alleles = ("-9", "-9")
                else:
                    # STRUCTURE accepts arbitrary positive integer allele labels.
                    alleles = tuple(str(int(x) + 1) for x in gt)
                genotypes[sample].append(alleles)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as out:
        out.write("\t\t" + "\t".join(markers) + "\n")
        for sample in samples:
            pop = pop_number[metadata[sample]["population"]]
            for chromosome in (0, 1):
                allele_row = [x[chromosome] for x in genotypes[sample]]
                out.write(f"{sample}\t{pop}\t" + "\t".join(allele_row) + "\n")

    with args.population_key.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["population_number", "population"])
        writer.writerows((number, name) for name, number in sorted(pop_number.items(), key=lambda x: x[1]))
    print(f"Wrote {len(samples)} individuals and {len(markers)} loci to {args.output}")


if __name__ == "__main__":
    main()
