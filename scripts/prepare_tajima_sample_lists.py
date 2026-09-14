#!/usr/bin/env python3
"""Create VCFtools population, variety, and analysis-group sample lists."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "data" / "sample_metadata.tsv"
OUTPUT_DIR = ROOT / "results" / "tajimas_d" / "sample_lists"


def safe_group(value: str) -> str:
    allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-"
    if not value or any(char not in allowed for char in value):
        raise ValueError(f"Unsafe or empty group label: {value!r}")
    return value


def main() -> None:
    groups: dict[str, list[str]] = defaultdict(list)
    seen: set[str] = set()
    with METADATA.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "population", "variety", "analysis_group"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{METADATA} must contain: {', '.join(sorted(required))}")
        for row in reader:
            sample = row["sample_id"].strip()
            if not sample or sample in seen:
                raise ValueError(f"Missing or duplicate sample ID: {sample!r}")
            seen.add(sample)
            population = safe_group(row["population"].strip())
            variety = safe_group(row["variety"].strip())
            analysis_group = safe_group(row["analysis_group"].strip())
            groups[f"population_{population}"].append(sample)
            groups[f"variety_{variety}"].append(sample)
            groups[f"group_{analysis_group}"].append(sample)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for label, samples in sorted(groups.items()):
        (OUTPUT_DIR / f"{label}.samples.txt").write_text(
            "\n".join(samples) + "\n", encoding="utf-8"
        )
    print(f"Wrote {len(groups)} sample lists for {len(seen)} samples to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
