#!/usr/bin/env python3
"""Summarize replicate STRUCTURE runs, calculate Evanno Delta K, and plot Q."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment


LIKELIHOOD = re.compile(r"Estimated Ln Prob of Data\s*=\s*([-+0-9.eE]+)")
QROW = re.compile(r"^\s*\d+\s+(\S+)\s+\(\s*\d+\)\s+\d+\s*:\s*(.+?)\s*$")
FILE_RE = re.compile(r"K(\d+)_rep(\d+)_f$")


def parse_run(path: Path) -> tuple[float, list[str], np.ndarray]:
    likelihood = None
    samples: list[str] = []
    values: list[list[float]] = []
    in_q = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = LIKELIHOOD.search(line)
        if match:
            likelihood = float(match.group(1))
        if line.startswith("Inferred ancestry of individuals:"):
            in_q = True
            continue
        if in_q:
            match = QROW.match(line)
            if match:
                samples.append(match.group(1))
                values.append([float(x) for x in match.group(2).split()])
            elif values and not line.strip():
                break
    if likelihood is None or not values:
        raise ValueError(f"Incomplete STRUCTURE output: {path}")
    return likelihood, samples, np.asarray(values)


def align_to_reference(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
    if reference.shape != query.shape:
        raise ValueError("Q matrices have different dimensions")
    cost = np.zeros((reference.shape[1], query.shape[1]))
    for i in range(reference.shape[1]):
        for j in range(query.shape[1]):
            cost[i, j] = np.sum((reference[:, i] - query[:, j]) ** 2)
    ref_cols, query_cols = linear_sum_assignment(cost)
    order = [query_cols[np.where(ref_cols == i)[0][0]] for i in range(reference.shape[1])]
    return query[:, order]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    runs = []
    for path in sorted(args.raw_dir.glob("K*/K*_rep*_f")):
        match = FILE_RE.search(path.name)
        if not match:
            continue
        k, rep = map(int, match.groups())
        likelihood, samples, q = parse_run(path)
        runs.append({"K": k, "replicate": rep, "likelihood": likelihood, "samples": samples, "q": q, "path": str(path)})
    counts = defaultdict(int)
    for run in runs:
        counts[run["K"]] += 1
    expected = set(range(1, 8))
    if set(counts) != expected or any(counts[k] < 2 for k in expected):
        raise SystemExit(f"Need >=2 completed runs for every K=1..7; found {dict(counts)}")

    with (args.output_dir / "structure_run_likelihoods.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["K", "replicate", "likelihood", "path"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for run in sorted(runs, key=lambda x: (x["K"], x["replicate"])):
            writer.writerow({key: run[key] for key in writer.fieldnames})

    means, sds = {}, {}
    for k in sorted(counts):
        vals = [x["likelihood"] for x in runs if x["K"] == k]
        means[k] = statistics.mean(vals)
        sds[k] = statistics.stdev(vals)
    rows = []
    for k in sorted(counts):
        lprime = means[k] - means[k - 1] if k > min(counts) else None
        ldouble = abs(means[k + 1] - 2 * means[k] + means[k - 1]) if min(counts) < k < max(counts) else None
        delta = ldouble / sds[k] if ldouble is not None and sds[k] > 0 else None
        rows.append({"K": k, "replicates": counts[k], "mean_LnP_K": means[k], "sd_LnP_K": sds[k],
                     "Ln_prime_K": lprime, "abs_Ln_double_prime_K": ldouble, "Delta_K": delta})
    with (args.output_dir / "evanno_delta_k.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    metadata_rows = list(csv.DictReader(args.metadata.open(encoding="utf-8"), delimiter="\t"))
    metadata = {x["sample_id"]: x for x in metadata_rows}
    consensus = {}
    for k in sorted(counts):
        k_runs = [x for x in runs if x["K"] == k]
        best = max(k_runs, key=lambda x: x["likelihood"])
        aligned = []
        for run in k_runs:
            if run["samples"] != best["samples"]:
                raise ValueError(f"Sample order mismatch at K={k}")
            aligned.append(align_to_reference(best["q"], run["q"]))
        mean_q = np.mean(aligned, axis=0)
        mean_q /= mean_q.sum(axis=1, keepdims=True)
        consensus[k] = (best["samples"], mean_q)

    with (args.output_dir / "structure_consensus_q.tsv").open("w", newline="", encoding="utf-8") as handle:
        max_k = max(counts)
        fields = ["K", "sample_id", "population", "variety"] + [f"cluster_{i}" for i in range(1, max_k + 1)]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for k, (samples, q) in consensus.items():
            for sample, values in zip(samples, q):
                row = {"K": k, "sample_id": sample, "population": metadata[sample]["population"], "variety": metadata[sample]["variety"]}
                row.update({f"cluster_{i + 1}": value for i, value in enumerate(values)})
                writer.writerow(row)

    variety_rows = []
    for k, (samples, q) in consensus.items():
        for variety in sorted({metadata[s]["variety"] for s in samples}):
            indices = [i for i, sample in enumerate(samples) if metadata[sample]["variety"] == variety]
            mean_values = np.mean(q[indices, :], axis=0)
            row = {"K": k, "variety": variety, "n": len(indices)}
            row.update({f"cluster_{i + 1}": value for i, value in enumerate(mean_values)})
            variety_rows.append(row)
    with (args.output_dir / "structure_variety_mean_q.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["K", "variety", "n"] + [f"cluster_{i}" for i in range(1, max(counts) + 1)]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(variety_rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].errorbar(sorted(means), [means[k] for k in sorted(means)], yerr=[sds[k] for k in sorted(means)], marker="o", capsize=3)
    axes[0].set(xlabel="K", ylabel="Mean Ln P(K) ± SD", xticks=sorted(means))
    delta_rows = [x for x in rows if x["Delta_K"] is not None]
    axes[1].plot([x["K"] for x in delta_rows], [x["Delta_K"] for x in delta_rows], marker="o")
    axes[1].set(xlabel="K", ylabel="Evanno ΔK", xticks=[x["K"] for x in delta_rows])
    fig.tight_layout()
    fig.savefig(args.output_dir / "structure_likelihood_deltaK.pdf")
    fig.savefig(args.output_dir / "structure_likelihood_deltaK.png", dpi=300)
    plt.close(fig)

    plot_ks = list(range(2, 8))
    fig, axes = plt.subplots(len(plot_ks), 1, figsize=(13, 8), sharex=True)
    colors = plt.get_cmap("tab10")(np.arange(7))
    for ax, k in zip(axes, plot_ks):
        samples, q = consensus[k]
        bottom = np.zeros(len(samples))
        for cluster in range(k):
            ax.bar(np.arange(len(samples)), q[:, cluster], bottom=bottom, width=1.0, color=colors[cluster], linewidth=0)
            bottom += q[:, cluster]
        ax.set_ylim(0, 1)
        ax.set_ylabel(f"K={k}", rotation=0, labelpad=20, va="center")
        previous = None
        centers = []
        start = 0
        for i, sample in enumerate(samples + [None]):
            pop = metadata[sample]["population"] if sample else None
            if previous is not None and pop != previous:
                ax.axvline(i - 0.5, color="white", linewidth=0.35)
                centers.append(((start + i - 1) / 2, previous))
                start = i
            previous = pop
        if k == plot_ks[-1]:
            ax.set_xticks([x[0] for x in centers], [x[1] for x in centers], rotation=90, fontsize=7)
        else:
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.set_yticks([0, 1])
    fig.supxlabel("Population (individuals in VCF order)")
    fig.supylabel("Ancestry proportion")
    fig.tight_layout()
    fig.savefig(args.output_dir / "structure_consensus_K2-K7.pdf")
    fig.savefig(args.output_dir / "structure_consensus_K2-K7.png", dpi=300)
    plt.close(fig)

    best_delta = max(delta_rows, key=lambda x: x["Delta_K"])
    best_k = best_delta["K"]
    best_means = [x for x in variety_rows if x["K"] == best_k]
    composition_lines = []
    for row in best_means:
        components = ", ".join(
            f"C{i}={float(row[f'cluster_{i}']):.3f}" for i in range(1, best_k + 1)
        )
        composition_lines.append(f"- {row['variety']} (n={row['n']}): {components}")
    readme = f"""# STRUCTURE replication and optimal-K summary

- Input: 180 individuals and 1,291 LD-pruned RAD SNPs.
- Runs: K=1-7, 10 independent replicates per K, 25,000 burn-in iterations and 100,000 MCMC iterations.
- Maximum Evanno Delta K: **K={best_k}** (Delta K = **{best_delta['Delta_K']:.3f}**).
- Replicate cluster labels were aligned to the highest-likelihood run at each K before averaging.

## Mean ancestry by named variety at K={best_k}

{chr(10).join(composition_lines)}

## Interpretation

Evanno Delta K identifies the strongest upper level of hierarchical structure and is undefined for K=1, so it should be interpreted together with the mean log likelihoods, replicate variability, ancestry plots, geography, and biological question. The alignment here uses minimum-distance assignment of cluster columns, a reproducible analogue of the label-switching correction performed by CLUMPP. CLUMPP aligns replicate Q matrices; it does not calculate Delta K.

The available VCF is a filtered, LD-pruned subset. These runs reproduce broad ancestry structure for the revision but should not be described as an independent genome-wide discovery analysis.
"""
    (args.output_dir / "README.md").write_text(readme, encoding="utf-8")
    print(f"Parsed {len(runs)} runs; maximum Delta K is K={best_delta['K']} ({best_delta['Delta_K']:.3f})")


if __name__ == "__main__":
    main()
