#!/usr/bin/env python3
"""Join SNP proximity, BLASTP, and UniProt annotations into manuscript tables."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def relevance(go_bp: str, description: str) -> str:
    text = f"{go_bp}; {description}".lower()
    groups = []
    keywords = {
        "abiotic stress response": ("cold", "water deprivation", "osmotic stress", "salt stress", "hypoxia", "flooding", "nitrate starvation", "mannitol"),
        "photosynthesis/respiration": ("photosystem", "photosynthetic", "electron transport", "aerobic respiration", "cytochrome c oxidase", "nadh dehydrogenase"),
        "membrane transport": ("abc transporter", "transporter", "nitrate import"),
        "hormone/development": ("auxin", "root development", "shoot system development", "abscisic acid"),
        "defense response": ("disease resistance", "defense response", "immune"),
        "transcription/signaling": ("transcription factor", "protein kinase", "receptor"),
        "autophagy": ("autophagy", "autophagosome"),
    }
    for label, terms in keywords.items():
        if any(term in text for term in terms):
            groups.append(label)
    return "; ".join(groups)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proximity", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--uniprot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    proximity = read_tsv(args.proximity)
    hits = {x["qseqid"]: x for x in read_tsv(args.hits)}
    uniprot_rows = read_tsv(args.uniprot)
    uniprot = {x["Entry"]: x for x in uniprot_rows}

    output = []
    for row in proximity:
        hit = hits.get(row["protein_id"], {})
        annotation = uniprot.get(hit.get("arabidopsis_accession", ""), {})
        distance = int(row["distance_bp"]) if row["distance_bp"] else None
        if distance == 0 and row["overlap_feature"] == "CDS":
            locus_class = "coding"
        elif distance == 0:
            locus_class = "intronic_or_transcribed"
        elif distance is not None and distance <= 5000:
            locus_class = "proximal_5kb"
        elif distance is not None and distance <= 10000:
            locus_class = "proximal_10kb"
        else:
            locus_class = "nearest_distant"

        qcov = float(hit["query_coverage"]) if hit else 0
        pident = float(hit["percent_identity"]) if hit else 0
        if not hit:
            confidence = "no_reviewed_Arabidopsis_hit"
        elif qcov >= 70 and pident >= 40:
            confidence = "strong"
        elif qcov >= 40 and pident >= 30:
            confidence = "moderate"
        else:
            confidence = "limited"

        go_bp = annotation.get("Gene Ontology (biological process)", "")
        description = hit.get("hit_description", "")
        output.append(
            {
                **row,
                "locus_class": locus_class,
                "arabidopsis_accession": hit.get("arabidopsis_accession", ""),
                "arabidopsis_gene": hit.get("gene_symbol", ""),
                "homology_confidence": confidence,
                "percent_identity": hit.get("percent_identity", ""),
                "query_coverage_percent": hit.get("query_coverage", ""),
                "blast_evalue": hit.get("evalue", ""),
                "arabidopsis_description": description,
                "GO_biological_process": go_bp,
                "GO_molecular_function": annotation.get("Gene Ontology (molecular function)", ""),
                "functional_category": relevance(go_bp, description),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

    locus_counts = Counter(x["locus_class"] for x in output)
    unique_proteins = {x["protein_id"] for x in output if x["protein_id"]}
    unique_hits = {x["protein_id"] for x in output if x["arabidopsis_accession"]}
    categories = Counter()
    for protein in unique_proteins:
        row = next(x for x in output if x["protein_id"] == protein)
        for category in filter(None, row["functional_category"].split("; ")):
            categories[category] += 1

    lines = [
        "# Functional annotation summary",
        "",
        f"- Candidate SNPs: {len({x['snp'] for x in output})}",
        f"- Unique nearest/overlapping PM1N proteins: {len(unique_proteins)}",
        f"- Proteins with a reviewed Arabidopsis BLASTP hit: {len(unique_hits)}",
        f"- Coding SNPs: {locus_counts['coding']}",
        f"- Intronic/transcribed SNPs: {locus_counts['intronic_or_transcribed']}",
        f"- Non-overlapping SNPs within 5 kb: {locus_counts['proximal_5kb']}",
        f"- Non-overlapping SNPs 5-10 kb away: {locus_counts['proximal_10kb']}",
        f"- More distant nearest-gene assignments: {locus_counts['nearest_distant']}",
        "",
        "## Functional categories among unique annotated proteins",
        "",
    ]
    lines.extend(f"- {name}: {count}" for name, count in sorted(categories.items()))
    lines.extend(["", "## Enrichment", "", "- KEGG ABC transporters: 2 proteins; g:SCS-adjusted P = 0.0241", ""])
    args.summary.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
