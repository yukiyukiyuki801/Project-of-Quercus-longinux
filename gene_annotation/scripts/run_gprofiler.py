#!/usr/bin/env python3
"""Run g:Profiler enrichment for unique Arabidopsis top-hit gene symbols."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import requests


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.hits.open(newline="", encoding="utf-8") as handle:
        genes = sorted({row["gene_symbol"] for row in csv.DictReader(handle, delimiter="\t") if row["gene_symbol"]})
    payload = {
        "organism": "athaliana",
        "query": genes,
        "sources": ["GO:MF", "GO:CC", "GO:BP", "KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "domain_scope": "annotated",
    }
    response = requests.post(
        "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
        json=payload,
        timeout=120,
    )
    response.raise_for_status()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(response.json(), indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
