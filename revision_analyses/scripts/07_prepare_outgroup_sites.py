#!/usr/bin/env python3
"""Translate known VCF loci from accessions to PM1N FASTA sequence names."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    args = parser.parse_args()

    accession_to_sequence = {}
    with args.assembly_report.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            accession_to_sequence[fields[4]] = fields[0]

    rows = []
    with args.vcf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            accession, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            sequence = accession_to_sequence.get(accession)
            if not sequence:
                raise SystemExit(f"No PM1N sequence-name mapping for {accession}")
            rows.append((sequence, pos - 1, pos, accession, ref, alt))

    args.bed.parent.mkdir(parents=True, exist_ok=True)
    with args.bed.open("w", encoding="utf-8") as handle:
        for sequence, start, end, *_ in sorted(rows):
            handle.write(f"{sequence}\t{start}\t{end}\n")
    with args.map.open("w", encoding="utf-8") as handle:
        handle.write("pm1n_sequence\tposition\taccession\tref\talt\n")
        for sequence, _, end, accession, ref, alt in rows:
            handle.write(f"{sequence}\t{end}\t{accession}\t{ref}\t{alt}\n")
    print(f"Wrote {len(rows)} known sites")


if __name__ == "__main__":
    main()
