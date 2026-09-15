#!/usr/bin/env python3
"""Audit target coverage using Table S2 identities and original guide UMIs only.

No survival phenotype or expression columns are used. Positive guide UMIs do
not establish a high-confidence perturbation assignment or biological replicate.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import warnings

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/v2/perturbation_feasibility"
WORKBOOK_SHA256 = "149ec75dda0adaae6238c00d3ac5ecd7a308a4d910e4f23359247eb0372e1247"
SHEET = "sgRNA_info_for_pooled_validatio"
FILES = {
    "GSM3543624_neuron_1_sgRNA_mapping.txt.gz": "bc2bb1718e9a6e884631c97a0750bf4c245dfbe8a2b198cd246511cf41cefd60",
    "GSM3543625_neuron_2_sgRNA_mapping.txt.gz": "3cfe2b3eaf8fac60c4eebe6b33e19754e7676706fa64bccc44e7002d223e760a",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mapping_directory", type=Path)
    parser.add_argument("--workbook", type=Path, default=ROOT /
                        "data_raw/perturbation_guides/NIHMS1535621-supplement-8.xlsx")
    parser.add_argument("--output-directory", type=Path, default=DOCS)
    args = parser.parse_args()
    if sha256(args.workbook) != WORKBOOK_SHA256:
        raise ValueError("Table S2 checksum differs from the uploaded source")
    manifest = json.loads((DOCS / "candidate_manifest.json").read_text())
    pathway = set(manifest["pathway_genes"])
    shared = set(manifest["prior_shared_leading_edge_genes"])
    if len(pathway) != 52 or len(shared) != 8 or not shared <= pathway:
        raise ValueError("Fixed candidate manifest is invalid")
    # The unsupported workbook extension is irrelevant to this read-only task.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unknown extension is not supported")
        wb = openpyxl.load_workbook(args.workbook, read_only=True, data_only=False)
        if wb.sheetnames != [SHEET]:
            raise ValueError("Unexpected workbook sheets")
        ws = wb[SHEET]
        # Only the first three identity columns are yielded. Never export the
        # seven subsequent survival-phenotype columns into an analysis table.
        rows = list(ws.iter_rows(min_col=1, max_col=3, values_only=True))
        dimensions = [ws.max_row, ws.max_column]
        wb.close()
    if list(rows[0]) != ["sgRNA_short name", "sgRNA_long name", "protospacer sequence"]:
        raise ValueError("Unexpected guide-identity header")
    guides = []
    for short, long_name, sequence in rows[1:]:
        if not all(isinstance(x, str) for x in (short, long_name, sequence)):
            raise ValueError("Missing or non-text identity")
        if not re.fullmatch("[ACGT]{20}", sequence):
            raise ValueError("Unexpected protospacer format")
        gene = "control" if short.startswith("non-targeting") else short.split("_")[0]
        guides.append({"guide": short, "gene": gene, "sequence": sequence})
    key = {r["sequence"]: r["gene"] for r in guides}
    if len(key) != len(guides):
        raise ValueError("Duplicated sequence identities")
    genes = set(key.values()) - {"control"}
    result = {
        "audit_date_utc": "2026-09-15",
        "source_paper": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6813890/",
        "source_workbook_sha256": WORKBOOK_SHA256,
        "source_dimensions": dimensions,
        "scope": "guide identity columns and guide-enrichment read/UMI counts only",
        "symbol_policy": "exact symbols; no inferred or unverified aliases",
        "guides": len(key), "noncontrol_gene_labels": len(genes),
        "nontargeting_guides": sum(v == "control" for v in key.values()),
        "workbook_pathway_targets": sorted(genes & pathway),
        "workbook_shared8_targets": sorted(genes & shared),
        "files": {},
        "decision": "NO_SUPPORTED_FIXED_TARGET_BENCHMARK",
        "limitations": [
            "Unidentified guide sequences are retained, not assigned guessed genes.",
            "No candidate targets are established in the mapped subset.",
            "This does not prove all unidentified sequences lack pathway targets.",
            "Guide detection is not expression QC or a biological replicate count.",
        ],
    }
    for name in FILES:
        path = args.mapping_directory / name
        if sha256(path) != FILES[name]:
            raise ValueError(f"Raw guide-mapping checksum mismatch: {name}")
        counts = Counter()
        detected = defaultdict(set)
        unknown_umis = Counter()
        with gzip.open(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["cell", "barcode", "read_count", "umi_count"]:
                raise ValueError("Unexpected raw mapping header")
            for row in reader:
                counts["rows"] += 1
                umi = int(row["umi_count"])
                if umi < 0 or int(row["read_count"]) < 0:
                    raise ValueError("Negative guide count")
                if umi == 0:
                    continue
                sequence = row["barcode"]
                group = ("key" if sequence in key else "unmapped_20nt"
                         if re.fullmatch("[ACGT]{20}", sequence) else "unprocessed")
                counts[group + "_rows_positive"] += 1
                counts[group + "_umis"] += umi
                detected[group].add(sequence)
                if group == "unmapped_20nt":
                    unknown_umis[sequence] += umi
        represented = {key[s] for s in detected["key"]} - {"control"}
        accession = name.split("_")[0]
        result["files"][name] = {
            "url": f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3543nnn/{accession}/suppl/{name}",
            "sha256": sha256(path),
            "counts": dict(sorted(counts.items())),
            "positive_distinct_sequences": {k: len(v) for k, v in sorted(detected.items())},
            "mapped_noncontrol_gene_labels": sorted(represented),
            "pathway_targets": sorted(represented & pathway),
            "shared8_targets": sorted(represented & shared),
            "unmapped_20nt_positive_umis": dict(sorted(unknown_umis.items())),
        }
    if sha256(args.workbook) != WORKBOOK_SHA256:
        raise ValueError("Source workbook changed during the audit")
    args.output_directory.mkdir(parents=True, exist_ok=True)
    (args.output_directory / "Tian2019_uploaded_guide_key_coverage.json").write_text(
        json.dumps(result, indent=2) + "\n")
    with (args.output_directory / "Tian2019_guide_identities.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["guide", "gene", "sequence"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(guides)
    print(json.dumps({k: v for k, v in result.items() if k != "files"}, indent=2))


if __name__ == "__main__":
    main()
