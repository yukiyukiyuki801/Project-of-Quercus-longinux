#!/usr/bin/env python3
"""Pairwise independent-group rank-sum tests with per-metric Holm correction."""

from __future__ import annotations

import argparse
import csv
from itertools import combinations
from pathlib import Path

from scipy.stats import mannwhitneyu


def holm(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [0.0] * len(values)
    running = 0.0
    total = len(values)
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (total - rank) * values[index]))
        adjusted[index] = running
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    groups = ("kuoi", "lativiolaciifolia", "longinux")
    output = []
    for metric in ("He", "Ho", "FIS", "pi"):
        tests = []
        for first, second in combinations(groups, 2):
            x = [float(row[metric]) for row in rows if row["variety"] == first]
            y = [float(row[metric]) for row in rows if row["variety"] == second]
            test = mannwhitneyu(x, y, alternative="two-sided", method="auto")
            tests.append(
                {
                    "metric": metric,
                    "group_1": first,
                    "group_2": second,
                    "n_1": len(x),
                    "n_2": len(y),
                    "W": test.statistic,
                    "p_raw": test.pvalue,
                }
            )
        for row, adjusted in zip(tests, holm([row["p_raw"] for row in tests])):
            row["p_holm"] = adjusted
            output.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["metric", "group_1", "group_2", "n_1", "n_2", "W", "p_raw", "p_holm"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in output:
            writer.writerow(
                {
                    **row,
                    "W": f"{row['W']:.6g}",
                    "p_raw": f"{row['p_raw']:.8g}",
                    "p_holm": f"{row['p_holm']:.8g}",
                }
            )
    print(args.output)


if __name__ == "__main__":
    main()
