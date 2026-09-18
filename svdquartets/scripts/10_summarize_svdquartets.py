#!/usr/bin/env python3
"""Root, plot, and summarize the SVDquartets population consensus tree."""

from __future__ import annotations

import argparse
import copy
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo


COLORS = {"kuoi": "#3A7D44", "longinux": "#2F6FA3", "lativiolaciifolia": "#8A5FBF"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--outgroup")
    parser.add_argument("--sampling", choices=("all_individuals", "representatives"), default="all_individuals")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    metadata_rows = list(csv.DictReader(args.metadata.open(encoding="utf-8"), delimiter="\t"))
    pop_variety = {}
    for row in metadata_rows:
        pop_variety.setdefault(row["population"], row["variety"])

    tree = Phylo.read(args.tree, "newick")
    # PAUP*'s "raw Newick" consensus export writes bootstrap percentages in
    # the branch-length slot (for example, ):100).  They are support values,
    # not evolutionary branch lengths.  Preserve terminal values only as an
    # export artefact and transfer internal values to Biopython confidence.
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.branch_length is not None:
            clade.confidence = clade.branch_length
        clade.branch_length = None
    kuoi = sorted(pop for pop, variety in pop_variety.items() if variety == "kuoi")
    kuoi_terminals = [next(tree.find_clades(name=name)) for name in kuoi]
    kuoi_clade = tree.common_ancestor(kuoi_terminals)
    kuoi_descendants = sorted(x.name for x in kuoi_clade.get_terminals())
    kuoi_support = kuoi_clade.confidence
    kuoi_monophyletic = kuoi_descendants == kuoi
    display_tree = copy.deepcopy(tree)
    if args.outgroup:
        display_tree.root_with_outgroup(next(display_tree.find_clades(name=args.outgroup)))
    elif kuoi_monophyletic:
        display_kuoi = display_tree.common_ancestor(
            [next(display_tree.find_clades(name=name)) for name in kuoi]
        )
        display_tree.root_with_outgroup(display_kuoi)

    for terminal in display_tree.get_terminals():
        terminal.color = COLORS.get(pop_variety.get(terminal.name, ""), "black")

    def branch_label(clade):
        # Show every internal split value, including weak support adjacent to
        # the outgroup, so rooting uncertainty is visible in the figure.
        if clade.is_terminal() or clade.confidence is None:
            return None
        return str(int(round(clade.confidence)))

    fig = plt.figure(figsize=(9, 10))
    ax = fig.add_subplot(111)
    Phylo.draw(display_tree, axes=ax, do_show=False, show_confidence=False, branch_labels=branch_label)
    ax.set_title("SVDquartets population tree (100 bootstrap replicates; all node supports)")
    ax.set_xlabel("Topology only (branch lengths not estimated)")
    outgroup_support_note = ""
    if args.outgroup:
        ingroup_child = next(
            clade for clade in display_tree.root.clades
            if args.outgroup not in {tip.name for tip in clade.get_terminals()}
        )
        if ingroup_child.confidence is not None:
            outgroup_support_note = (
                "- The figure displays bootstrap support for every internal split. "
                f"The split adjacent to the *Q. glauca* root has **{ingroup_child.confidence:.0f}%** support; "
                "the post-hoc root itself has no separate bootstrap value.\n"
            )
            ax.text(
                0.99,
                0.985,
                f"Outgroup-adjacent split: {ingroup_child.confidence:.0f}% bootstrap",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
            )
    legend_handles = [
        Patch(color=COLORS["kuoi"], label="var. kuoi"),
        Patch(color=COLORS["longinux"], label="var. longinux"),
        Patch(color=COLORS["lativiolaciifolia"], label="var. lativiolaciifolia"),
    ]
    if args.outgroup:
        legend_handles.append(Patch(color="black", label="Q. glauca (outgroup)"))
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output_dir / "svdquartets_population_tree.pdf")
    fig.savefig(args.output_dir / "svdquartets_population_tree.png", dpi=300)
    plt.close(fig)

    # With an explicit outgroup, evaluate rooted variety clades on the rooted
    # display copy so MRCA membership is biologically interpretable.
    evaluation_tree = display_tree if args.outgroup else tree
    variety_rows = []
    group_specs = [
        ("kuoi", lambda value: value == "kuoi"),
        ("longinux", lambda value: value == "longinux"),
        ("lativiolaciifolia", lambda value: value == "lativiolaciifolia"),
        ("non_kuoi_combined", lambda value: value != "kuoi"),
    ]
    for variety, selector in group_specs:
        pops = sorted(pop for pop, value in pop_variety.items() if selector(value))
        terminals = [next(evaluation_tree.find_clades(name=name)) for name in pops]
        clade = evaluation_tree.common_ancestor(terminals)
        descendants = sorted(x.name for x in clade.get_terminals())
        variety_rows.append(
            {
                "variety": variety,
                "populations": ",".join(pops),
                "monophyletic": descendants == pops,
                "crown_support": clade.confidence if descendants == pops and clade.confidence is not None else "NA",
                "mrca_descendants": ",".join(descendants),
            }
        )
    with (args.output_dir / "svdquartets_variety_monophyly.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(variety_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(variety_rows)

    support_text = f"{kuoi_support:.0f}%" if kuoi_support is not None else "not available"
    non_kuoi_monophyletic = next(x["monophyletic"] for x in variety_rows if x["variety"] == "non_kuoi_combined")
    if args.sampling == "representatives":
        input_text = "27 individuals (one lowest-missing representative per ingroup population plus DH01), treated as 27 tips"
        method_text = "PAUP* SVDquartets, all 17,550 possible quartets, and 100 standard bootstrap replicates"
    else:
        input_text = "181 individuals, 27 population/species tips" if args.outgroup else "180 individuals, 26 population tips"
        method_text = "PAUP* SVDquartets under the multispecies coalescent, population taxon partition, 100,000 randomly sampled lineage quartets, and 100 standard bootstrap replicates"
    summary = f"""# SVDquartets result summary

- Input: {input_text}, and 1,291 LD-pruned RAD SNPs.
- Method: {method_text}.
- Var. *kuoi* populations ({', '.join(kuoi)}) are monophyletic: **{kuoi_monophyletic}**.
- Bootstrap support for the complete var. *kuoi* population clade: **{support_text}**.
{outgroup_support_note.rstrip()}
- The combined var. *longinux* + var. *lativiolaciifolia* group is monophyletic after rooting on *Q. glauca*: **{non_kuoi_monophyletic}**.
- Var. *longinux* and var. *lativiolaciifolia* are not reciprocally monophyletic (see `svdquartets_variety_monophyly.tsv`).
"""
    (args.output_dir / "README.md").write_text(summary, encoding="utf-8")
    print(f"kuoi_monophyletic={kuoi_monophyletic}; bootstrap={kuoi_support}")


if __name__ == "__main__":
    main()
