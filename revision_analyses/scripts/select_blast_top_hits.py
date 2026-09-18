#!/usr/bin/env python3
"""Select the highest-bit-score reviewed Arabidopsis BLASTP hit per query."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fields = [
        "qseqid", "sseqid", "percent_identity", "alignment_length",
        "query_length", "subject_length", "query_coverage", "evalue",
        "bitscore", "title",
    ]
    hits = []
    with args.blast.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, fieldnames=fields, delimiter="\t"):
            row["_bitscore"] = float(row["bitscore"])
            hits.append(row)

    best = {}
    for row in hits:
        if row["qseqid"] not in best or row["_bitscore"] > best[row["qseqid"]]["_bitscore"]:
            best[row["qseqid"]] = row

    output_fields = [
        "qseqid", "arabidopsis_accession", "gene_symbol", "percent_identity",
        "alignment_length", "query_length", "subject_length", "query_coverage",
        "evalue", "bitscore", "hit_description",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for query in sorted(best):
            row = best[query]
            accession = row["sseqid"].split("|")[1] if "|" in row["sseqid"] else row["sseqid"]
            gene_match = re.search(r"\bGN=([^ ]+)", row["title"])
            description = row["title"].split(" OS=")[0]
            description = re.sub(r"^\S+\s+", "", description)
            writer.writerow(
                {
                    "qseqid": query,
                    "arabidopsis_accession": accession,
                    "gene_symbol": gene_match.group(1) if gene_match else "",
                    "percent_identity": row["percent_identity"],
                    "alignment_length": row["alignment_length"],
                    "query_length": row["query_length"],
                    "subject_length": row["subject_length"],
                    "query_coverage": row["query_coverage"],
                    "evalue": row["evalue"],
                    "bitscore": row["bitscore"],
                    "hit_description": description,
                }
            )


if __name__ == "__main__":
    main()
