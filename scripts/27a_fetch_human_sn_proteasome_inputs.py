"""Fetch checksum-locked inputs for script 27.

The tracked GSE178265 DA-cell metadata was frozen before target-expression
access. This fetcher obtains only the public cluster ordering and expression
values needed to pair with that metadata.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import time
from urllib.parse import urlencode, quote
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw" / "human_sn_proteasome"
GMT = ROOT / "data_raw" / "pathways" / "ReactomePathways.gmt"
KAMATH_META = ROOT / "metadata" / "human_sn_proteasome" / "GSE178265_DA_cell_metadata.tsv.gz"
SCP_API = "https://singlecell.broadinstitute.org/single_cell/api/v1"

GSE243639 = {
    "GSE243639_Clinical_data.csv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243639/suppl/GSE243639_Clinical_data.csv.gz",
        "ffea1163eb0c145d452fc08f8d4550986353c29b8b9496f1b69fa208296ebe24",
    ),
    "GSE243639_Filtered_count_table.csv.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243639/suppl/GSE243639_Filtered_count_table.csv.gz",
        "8a9fd3b08d4357ef2a667b28a9c7c0f658f57782df57f60d2145d4870a6d7f7f",
    ),
    "GSE243639_UMAP_coordinates.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243639/suppl/GSE243639_UMAP_coordinates.xlsx",
        "2fcd06645d8b5e4a32ec4462310b5571df5b657ce1e29c651c19b493a63035f4",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def request_bytes(url: str, attempts: int = 3) -> bytes:
    error = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"Accept": "application/json", "User-Agent": "pd-envtox/aan-v2"})
            with urlopen(request, timeout=180) as response:
                return response.read()
        except Exception as caught:  # noqa: BLE001
            error = caught
            if attempt + 1 < attempts:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"Download failed after {attempts} attempts: {url}") from error


def download_verified(path: Path, url: str, expected: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and sha256(path) == expected:
        print(f"Verified {path.name}")
        return
    partial = path.with_suffix(path.suffix + ".download")
    error = None
    for attempt in range(3):
        try:
            request = Request(url, headers={"User-Agent": "pd-envtox/aan-v2"})
            with urlopen(request, timeout=180) as response, partial.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
            break
        except Exception as caught:  # noqa: BLE001
            error = caught
            partial.unlink(missing_ok=True)
            if attempt < 2:
                time.sleep(2 ** attempt)
    else:
        raise RuntimeError(f"Download failed after 3 attempts: {url}") from error
    observed = sha256(partial)
    if observed != expected:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"SHA256 mismatch for {path.name}: {observed}")
    partial.replace(path)
    print(f"Downloaded and verified {path.name}")


def genes() -> list[str]:
    with GMT.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[:2] == ["Proteasome assembly", "R-HSA-9907900"]:
                result = [x.upper() for x in fields[2:]]
                if len(result) != 52:
                    raise RuntimeError("Frozen pathway does not contain 52 genes")
                return result
    raise RuntimeError("Frozen Reactome pathway is missing")


def valid_gene_response(path: Path, gene: str, n_cells: int) -> bool:
    try:
        payload = json.loads(path.read_text())
        return payload.get("genes") == [gene] and len(payload.get("data", {}).get("expression", [])) == n_cells
    except (OSError, ValueError):
        return False


def fetch_gene(gene: str, destination: Path, n_cells: int) -> None:
    path = destination / f"{gene}.json"
    if valid_gene_response(path, gene, n_cells):
        return
    cluster = quote("UMAP: Human DA Neurons", safe="")
    query = urlencode({"fields": "expression", "gene": gene})
    url = f"{SCP_API}/studies/SCP1768/clusters/{cluster}?{query}"
    content = request_bytes(url)
    partial = path.with_suffix(".json.download")
    partial.write_bytes(content)
    if not valid_gene_response(partial, gene, n_cells):
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Invalid SCP1768 expression response for {gene}")
    partial.replace(path)


def main() -> None:
    gse178265 = RAW / "gse178265"
    gse243639 = RAW / "gse243639"
    gene_dir = gse178265 / "genes"
    gene_dir.mkdir(parents=True, exist_ok=True)
    gse243639.mkdir(parents=True, exist_ok=True)

    for filename, (url, expected) in GSE243639.items():
        download_verified(gse243639 / filename, url, expected)

    cluster = quote("UMAP: Human DA Neurons", safe="")
    cluster_url = f"{SCP_API}/studies/SCP1768/clusters/{cluster}?{urlencode({'fields': 'cells,annotation'})}"
    cluster_path = gse178265 / "SCP1768_DA_cluster.json"
    cluster_path.write_bytes(request_bytes(cluster_url))
    payload = json.loads(cluster_path.read_text())
    cells = payload.get("data", {}).get("cells", [])
    annotations = payload.get("data", {}).get("annotations", [])
    tracked = pd.read_csv(KAMATH_META, sep="\t", usecols=["cell", "subtype"])
    if cells != tracked["cell"].tolist() or annotations != tracked["subtype"].tolist():
        raise RuntimeError("Current SCP1768 DA cluster order differs from the frozen metadata")

    target = genes()
    with ThreadPoolExecutor(max_workers=6) as executor:
        list(executor.map(lambda gene: fetch_gene(gene, gene_dir, len(cells)), target))
    failures = [gene for gene in target if not valid_gene_response(gene_dir / f"{gene}.json", gene, len(cells))]
    if failures:
        raise RuntimeError(f"Invalid gene responses after fetch: {failures}")
    print(f"Verified SCP1768 cluster and {len(target)} target-gene responses")


if __name__ == "__main__":
    main()
