#!/usr/bin/env python3
"""Check GSE152988 target coverage without opening expression or phenotype data.

Uses the authors' sgRNA-to-cell assignment archive linked from
https://kampmannlab.ucsf.edu/crop-seq. This is a feasibility audit, not an
expression QC step, statistical test, or independent-replication result.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/v2/perturbation_feasibility"
ARCHIVE_SHA256 = "516d0b462e340bce4a08f0a5796b52c46d00f2af15cd0742dc1f857c22c20c45"
ARCHIVE_URL = "https://www.dropbox.com/scl/fi/za44p0krdgx2yqjn3ov7s/maps.zip?rlkey=t54xns72h7y7jb1a4xv5rlf27&dl=1"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Downloaded authors' maps.zip")
    parser.add_argument("--output", type=Path, default=DOCS / "GSE152988_target_coverage.json")
    args = parser.parse_args()
    digest = hashlib.sha256(args.archive.read_bytes()).hexdigest()
    if digest != ARCHIVE_SHA256:
        raise ValueError("Archive checksum differs from the audited version")
    manifest = json.loads((DOCS / "candidate_manifest.json").read_text())
    genes = set(manifest["pathway_genes"])
    shared = set(manifest["prior_shared_leading_edge_genes"])
    if len(genes) != 52 or len(shared) != 8 or not shared <= genes:
        raise ValueError("Candidate manifest membership is invalid")
    result = {
        "source_url": "https://kampmannlab.ucsf.edu/crop-seq",
        "archive_url": ARCHIVE_URL,
        "archive_sha256": digest,
        "retrieved_utc": "2026-09-13",
        "scope": "sgRNA-to-cell assignment labels only; no expression or screen phenotypes read",
        "files": {}, "modalities": {},
    }
    by_mode = defaultdict(Counter)
    with zipfile.ZipFile(args.archive) as archive:
        for name in archive.namelist():
            if not name.startswith("maps/") or not name.endswith(".txt"):
                continue
            content = archive.read(name)
            counts = Counter()
            for line in content.decode().splitlines():
                fields = line.split("\t")
                if len(fields) != 2:
                    raise ValueError(f"Unexpected assignment format: {name}")
                match = re.fullmatch(r"(.+)_[ia]\d+", fields[1])
                if match is None:
                    raise ValueError(f"Unrecognized guide label: {fields[1]}")
                counts[match.group(1)] += 1
            mode = "CRISPRa" if "CRISPRa" in name else "CRISPRi"
            by_mode[mode].update(counts)
            result["files"][name] = {
                "sha256": hashlib.sha256(content).hexdigest(),
                "assigned_cells": sum(counts.values()),
                "distinct_target_or_control_labels": len(counts),
                "candidate_assigned_cells": {g: counts[g] for g in sorted(genes)},
            }
    if len(result["files"]) != 6 or set(by_mode) != {"CRISPRi", "CRISPRa"}:
        raise ValueError("Unexpected archive file set")
    for mode, counts in by_mode.items():
        result["modalities"][mode] = {
            "assigned_cells": sum(counts.values()),
            "distinct_target_or_control_labels": len(counts),
            "pathway_covered": sorted(genes & counts.keys()),
            "shared8_covered": sorted(shared & counts.keys()),
            "candidate_assigned_cells": {g: counts[g] for g in sorted(genes)},
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    for mode, item in result["modalities"].items():
        print(f"{mode}: pathway {len(item['pathway_covered'])}/52; shared {len(item['shared8_covered'])}/8")


if __name__ == "__main__":
    main()
