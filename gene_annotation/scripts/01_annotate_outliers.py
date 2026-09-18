#!/usr/bin/env python3
"""Annotate FST/RDA outlier SNPs against the matching Q. robur PM1N GFF.

The SNP table uses GenBank accessions (LT*/OLKR*), while the official PM1N
annotation uses assembly sequence names (Qrob_Chr*/Qrob_H2.3_Sc*).  The NCBI
assembly report is therefore required to translate coordinates without a
liftover.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_attrs(text: str) -> dict[str, str]:
    out = {}
    for item in text.rstrip().split(";"):
        if not item:
            continue
        key, _, value = item.partition("=")
        out[key] = value
    return out


def transcript_to_protein(transcript: str) -> str:
    return re.sub(r"^Qrob_T", "Qrob_P", transcript)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outliers", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protein-output", type=Path, required=True)
    args = parser.parse_args()

    accession_to_sequence: dict[str, str] = {}
    with args.assembly_report.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            accession_to_sequence[fields[4]] = fields[0]

    transcripts: dict[str, list[dict[str, object]]] = defaultdict(list)
    child_features: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with open_text(args.gff) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            seqid, _, feature, start, end, _, strand, _, attr_text = fields
            attrs = parse_attrs(attr_text)
            start_i, end_i = int(start), int(end)
            if feature == "mRNA":
                tid = attrs.get("ID", "")
                if tid:
                    transcripts[seqid].append(
                        {
                            "id": tid,
                            "start": start_i,
                            "end": end_i,
                            "strand": strand,
                            "quality": attrs.get("gene_qual", ""),
                        }
                    )
            else:
                parent = attrs.get("Parent", "")
                if parent:
                    child_features[parent].append((start_i, end_i, feature))

    for values in transcripts.values():
        values.sort(key=lambda x: (int(x["start"]), int(x["end"]), str(x["id"])))

    protein_ids: set[str] = set()
    output_rows = []
    with args.outliers.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            accession, pos_text = row["snp"].rsplit("_", 1)
            pos = int(pos_text)
            seqid = accession_to_sequence.get(accession, "")
            candidates = transcripts.get(seqid, [])
            if not candidates:
                output_rows.append(
                    {**row, "accession": accession, "position": pos, "pm1n_sequence": seqid,
                     "relation": "no_annotation", "distance_bp": "", "transcript_id": "",
                     "protein_id": "", "strand": "", "gene_quality": "", "overlap_feature": ""}
                )
                continue

            scored = []
            for tx in candidates:
                start, end = int(tx["start"]), int(tx["end"])
                distance = start - pos if pos < start else pos - end if pos > end else 0
                scored.append((distance, start, end, str(tx["id"]), tx))
            min_distance = min(x[0] for x in scored)
            nearest = [x[-1] for x in scored if x[0] == min_distance]
            for tx in nearest:
                tid = str(tx["id"])
                protein_id = transcript_to_protein(tid)
                protein_ids.add(protein_id)
                feature_hits = sorted(
                    {feature for start, end, feature in child_features.get(tid, []) if start <= pos <= end}
                )
                relation = "overlap" if min_distance == 0 else "nearest"
                output_rows.append(
                    {**row, "accession": accession, "position": pos, "pm1n_sequence": seqid,
                     "relation": relation, "distance_bp": min_distance, "transcript_id": tid,
                     "protein_id": protein_id, "strand": tx["strand"], "gene_quality": tx["quality"],
                     "overlap_feature": ",".join(feature_hits) if feature_hits else ("intron_or_mRNA" if min_distance == 0 else "")}
                )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0])
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    found: set[str] = set()
    args.protein_output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.proteins) as source, args.protein_output.open("w", encoding="utf-8") as target:
        keep = False
        for line in source:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                keep = pid in protein_ids
                if keep:
                    found.add(pid)
            if keep:
                target.write(line)

    missing = protein_ids - found
    print(f"Wrote {len(output_rows)} SNP-gene relationships for {len(protein_ids)} proteins")
    print(f"Extracted {len(found)} proteins; missing protein records: {len(missing)}")
    if missing:
        print("Missing:", ",".join(sorted(missing)))


if __name__ == "__main__":
    main()
