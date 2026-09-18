#!/usr/bin/env python3
"""Download reviewed-hit annotations from the UniProt REST API."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import requests


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.hits.open(newline="", encoding="utf-8") as handle:
        accessions = sorted({row["arabidopsis_accession"] for row in csv.DictReader(handle, delimiter="\t") if row["arabidopsis_accession"]})
    query = " OR ".join(f"accession:{item}" for item in accessions)
    response = requests.get(
        "https://rest.uniprot.org/uniprotkb/search",
        params={
            "query": f"({query})",
            "format": "tsv",
            "fields": "accession,id,gene_names,protein_name,go_p,go_f,go_c",
            "size": 500,
        },
        timeout=120,
    )
    response.raise_for_status()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(response.text, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
