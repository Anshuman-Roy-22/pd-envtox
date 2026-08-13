"""Frozen donor-level human SNpc proteasome-assembly validation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw" / "human_sn_proteasome"
GMT = ROOT / "data_raw" / "pathways" / "ReactomePathways.gmt"
OUT = ROOT / "results" / "v2" / "human_sn_proteasome"
OUT.mkdir(parents=True, exist_ok=True)
KAMATH_CELL_META = ROOT / "metadata" / "human_sn_proteasome" / "GSE178265_DA_cell_metadata.tsv.gz"

GMT_SHA256 = "89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b"
KAMATH_CELL_META_SHA256 = "6732fb4323cc36bad0a23bb53391d3242a2d68a1ebef5dd3d36857a72c2ab916"
GSE243639_SHA256 = {
    "GSE243639_Clinical_data.csv.gz": "ffea1163eb0c145d452fc08f8d4550986353c29b8b9496f1b69fa208296ebe24",
    "GSE243639_Filtered_count_table.csv.gz": "8a9fd3b08d4357ef2a667b28a9c7c0f658f57782df57f60d2145d4870a6d7f7f",
    "GSE243639_UMAP_coordinates.xlsx": "2fcd06645d8b5e4a32ec4462310b5571df5b657ce1e29c651c19b493a63035f4",
}
SHARED8 = ["POMP", "PSMA4", "PSMB4", "PSMC5", "PSMD1", "PSMD2", "PSMD4", "PSMG1"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, filename: str) -> None:
    frame.to_csv(OUT / filename, sep="\t", index=False, na_rep="NA")


def target_genes() -> list[str]:
    if sha256(GMT) != GMT_SHA256:
        raise RuntimeError("Tracked Reactome GMT failed its frozen SHA256 check")
    with GMT.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "Proteasome assembly" and fields[1] == "R-HSA-9907900":
                genes = [x.upper() for x in fields[2:]]
                if len(genes) != 52 or len(set(genes)) != 52:
                    raise RuntimeError("Frozen Proteasome assembly membership is not 52 unique genes")
                return genes
    raise RuntimeError("Frozen Proteasome assembly pathway is missing")


def load_kamath(genes: list[str]) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    base = RAW / "gse178265"
    cluster = json.loads((base / "SCP1768_DA_cluster.json").read_text())
    cells = cluster["data"]["cells"]
    subtypes = cluster["data"]["annotations"]
    if len(cells) != 22048 or len(subtypes) != len(cells):
        raise RuntimeError("SCP1768 DA cluster no longer contains the locked 22,048 cells")

    if sha256(KAMATH_CELL_META) != KAMATH_CELL_META_SHA256:
        raise RuntimeError("Tracked GSE178265 DA-cell metadata failed SHA256")
    meta = pd.read_csv(KAMATH_CELL_META, sep="\t", dtype={"donor": str})
    expected_columns = ["cell", "subtype", "donor", "status", "sex", "age", "pmi", "total_umi"]
    if list(meta.columns) != expected_columns or meta["cell"].tolist() != cells:
        raise RuntimeError("Tracked GSE178265 metadata does not match the SCP1768 DA cluster order")
    if meta["subtype"].tolist() != subtypes or meta[expected_columns].isna().any().any():
        raise RuntimeError("Tracked GSE178265 DA annotations are incomplete")
    meta["age"] = pd.to_numeric(meta["age"])
    meta["pmi"] = pd.to_numeric(meta["pmi"])
    meta["total_umi"] = meta["total_umi"].astype(np.int64)

    columns = []
    for gene in genes:
        path = base / "genes" / f"{gene}.json"
        payload = json.loads(path.read_text())
        values = payload.get("data", {}).get("expression", [])
        if payload.get("genes") != [gene] or len(values) != len(cells):
            raise RuntimeError(f"Invalid SCP1768 expression response for {gene}")
        array = np.asarray(values, dtype=np.int64)
        if np.any(array < 0):
            raise RuntimeError(f"Negative deposited counts for {gene}")
        columns.append(array)
    counts = np.column_stack(columns)
    return meta, counts, meta["total_umi"].to_numpy(np.int64)


def load_gse243639_cells() -> tuple[pd.DataFrame, dict[str, int]]:
    base = RAW / "gse243639"
    workbook = openpyxl.load_workbook(base / "GSE243639_UMAP_coordinates.xlsx", read_only=True, data_only=True)
    sheet = workbook["Neurons"]
    rows = sheet.iter_rows(values_only=True)
    header = {name: i for i, name in enumerate(next(rows))}
    records = []
    for row in rows:
        cell = str(row[header["CELL_ID"]])
        records.append({"cell": cell, "donor": cell.split("_", 1)[0], "subtype": str(row[header["IDENT"]])})
    meta = pd.DataFrame(records)
    if len(meta) != 4661 or meta["cell"].duplicated().any():
        raise RuntimeError("Unexpected GSE243639 neuron cell table")

    clinical_path = base / "GSE243639_Clinical_data.csv.gz"
    with gzip.open(clinical_path, "rt", newline="") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("N;"))
    reader = csv.DictReader(lines[start:], delimiter=";")
    clinical = []
    for row in reader:
        clinical.append({
            "donor": row["Sample ID"],
            "status": row["Clinical diagnosis"],
            "age": float(row["Age"]),
            "sex": row["Sex"],
            "pmi": float(row["PMI hours"]),
            "rin": float(row["RIN measure"]),
        })
    clinical_frame = pd.DataFrame(clinical)
    if len(clinical_frame) != 29 or clinical_frame["donor"].duplicated().any():
        raise RuntimeError("Unexpected GSE243639 clinical table")
    meta = meta.merge(clinical_frame, on="donor", validate="many_to_one")
    return meta, {cell: i for i, cell in enumerate(meta["cell"])}


def matrix_cell_name(name: str) -> str:
    return f"{name[:-2]}-1" if name.endswith(".1") else name


def load_gse243639_counts(meta: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    base = RAW / "gse243639"
    matrix_path = base / "GSE243639_Filtered_count_table.csv.gz"
    cache = base / "GSE243639_neuron_proteasome_cache.npz"
    expected_cells = meta["cell"].to_numpy(str)
    expected_genes = np.asarray(genes)
    if cache.exists():
        saved = np.load(cache, allow_pickle=False)
        if np.array_equal(saved["cells"], expected_cells) and np.array_equal(saved["genes"], expected_genes):
            return saved["counts"].astype(np.int64), saved["total_umi"].astype(np.int64)

    with gzip.open(matrix_path, "rb") as handle:
        header_line = handle.readline().decode("utf-8")
        header = next(csv.reader([header_line]))
        matrix_cells = [matrix_cell_name(x) for x in header[1:]]
        positions = {cell: i for i, cell in enumerate(matrix_cells)}
        try:
            selected = np.asarray([positions[x] for x in expected_cells], dtype=np.int64)
        except KeyError as error:
            raise RuntimeError(f"Neuron cell is absent from GSE243639 count matrix: {error}") from error
        total_umi = np.zeros(len(selected), dtype=np.int64)
        counts = np.zeros((len(selected), len(genes)), dtype=np.int64)
        gene_index = {gene: i for i, gene in enumerate(genes)}
        seen = set()
        for line_number, line in enumerate(handle, start=2):
            comma = line.find(b",")
            if comma < 3:
                raise RuntimeError(f"Malformed GSE243639 count row {line_number}")
            gene = line[:comma].strip().strip(b'"').decode("utf-8").upper()
            values = np.fromstring(line[comma + 1 :], dtype=np.int64, sep=",")
            if len(values) != len(matrix_cells):
                raise RuntimeError(f"Count-field mismatch on GSE243639 row {line_number}")
            selected_values = values[selected]
            total_umi += selected_values
            if gene in gene_index:
                counts[:, gene_index[gene]] += selected_values
                seen.add(gene)
    if seen != set(genes):
        raise RuntimeError(f"GSE243639 is missing frozen genes: {sorted(set(genes) - seen)}")
    np.savez_compressed(cache, cells=expected_cells, genes=expected_genes, counts=counts, total_umi=total_umi)
    return counts, total_umi


def aggregate_donors(
    cell_meta: pd.DataFrame,
    cell_counts: np.ndarray,
    cell_umi: np.ndarray,
    condition: pd.Series,
    min_cells: int,
    selection: np.ndarray | None = None,
) -> tuple[pd.DataFrame, np.ndarray]:
    if selection is None:
        selection = np.ones(len(cell_meta), dtype=bool)
    work = cell_meta.copy()
    work["pd"] = condition
    keep = selection & work["pd"].notna().to_numpy()
    donor_cell_counts = work.loc[keep].groupby("donor", sort=True).size()
    eligible = donor_cell_counts[donor_cell_counts >= min_cells].index.tolist()
    records = []
    aggregated = []
    for donor in eligible:
        idx = keep & (work["donor"].to_numpy() == donor)
        donor_rows = work.loc[idx]
        unique = donor_rows[["pd", "status", "age", "sex", "pmi"]].drop_duplicates()
        if len(unique) != 1:
            raise RuntimeError(f"Non-unique donor metadata for {donor}")
        record = unique.iloc[0].to_dict()
        record.update({"donor": donor, "n_cells": int(idx.sum()), "total_umi": int(cell_umi[idx].sum())})
        if "rin" in donor_rows:
            rin = donor_rows["rin"].drop_duplicates()
            if len(rin) != 1:
                raise RuntimeError(f"Non-unique RIN for {donor}")
            record["rin"] = float(rin.iloc[0])
        records.append(record)
        aggregated.append(cell_counts[idx].sum(axis=0))
    if not records:
        columns = ["pd", "status", "age", "sex", "pmi", "donor", "n_cells", "total_umi"]
        return pd.DataFrame(columns=columns), np.empty((0, cell_counts.shape[1]), dtype=np.int64)
    donor_meta = pd.DataFrame(records).sort_values(["pd", "donor"]).reset_index(drop=True)
    # Reorder the count rows to the already sorted donor metadata.
    count_map = {donor: row for donor, row in zip(eligible, aggregated)}
    donor_counts = np.vstack([count_map[x] for x in donor_meta["donor"]]).astype(np.int64)
    return donor_meta, donor_counts


def fit_ols(
    data: pd.DataFrame, outcome: np.ndarray, covariates: list[str]
) -> dict[str, float | int | str]:
    columns = [np.ones(len(data)), data["pd"].to_numpy(float)]
    names = ["intercept", "PD"]
    for covariate in covariates:
        if covariate == "sex":
            columns.append((data["sex"].str.lower() == "male").to_numpy(float))
            names.append("sex_male")
        else:
            values = data[covariate].to_numpy(float)
            columns.append(values - values.mean())
            names.append(f"centered_{covariate}")
    design = np.column_stack(columns)
    if not np.isfinite(design).all() or not np.isfinite(outcome).all():
        raise RuntimeError("Model contains non-finite values")
    rank = np.linalg.matrix_rank(design)
    if rank != design.shape[1]:
        raise RuntimeError("Frozen model design is not full rank")
    degrees = len(data) - rank
    if degrees <= 0:
        raise RuntimeError("Frozen model has no residual degrees of freedom")
    inverse = np.linalg.inv(design.T @ design)
    beta = inverse @ design.T @ outcome
    residual = outcome - design @ beta
    sigma2 = float(residual @ residual / degrees)
    standard_error = np.sqrt(np.diag(inverse) * sigma2)
    t_stat = beta / standard_error
    pd_index = names.index("PD")
    critical = stats.t.ppf(0.975, degrees)
    return {
        "beta_PD": float(beta[pd_index]),
        "se_PD": float(standard_error[pd_index]),
        "t_PD": float(t_stat[pd_index]),
        "df_residual": int(degrees),
        "p_one_sided_down": float(stats.t.cdf(t_stat[pd_index], degrees)),
        "p_two_sided": float(2 * stats.t.sf(abs(t_stat[pd_index]), degrees)),
        "ci95_low": float(beta[pd_index] - critical * standard_error[pd_index]),
        "ci95_high": float(beta[pd_index] + critical * standard_error[pd_index]),
        "model": "pathway_score ~ PD" + "".join(f" + {x}" for x in covariates),
    }


def label_result(beta: float, p_down: float, measured_n: int) -> str:
    if measured_n < 10:
        return "NOT_EVALUABLE"
    if beta < 0 and p_down < 0.05:
        return "CONFIRMED"
    if beta < 0 and p_down < 0.10:
        return "PARTIAL_SUPPORT"
    return "NOT_CONFIRMED"


def analyze(
    analysis: str,
    cohort: str,
    donor_meta: pd.DataFrame,
    donor_counts: np.ndarray,
    genes: list[str],
    covariates: list[str],
    confirmatory: bool,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    measured = (donor_counts >= 10).sum(axis=0) >= 3
    logcpm = np.log2(((donor_counts + 0.5) / (donor_meta["total_umi"].to_numpy()[:, None] + 1)) * 1e6)
    if measured.sum() < 10:
        raise RuntimeError(f"{analysis} has fewer than 10 measured genes")
    score = logcpm[:, measured].mean(axis=1)
    model = fit_ols(donor_meta, score, covariates)
    model.update({
        "analysis": analysis,
        "cohort": cohort,
        "confirmatory": confirmatory,
        "n_PD": int((donor_meta["pd"] == 1).sum()),
        "n_control": int((donor_meta["pd"] == 0).sum()),
        "measured_genes_n": int(measured.sum()),
    })
    model["outcome"] = label_result(model["beta_PD"], model["p_one_sided_down"], int(measured.sum()))

    donor_out = donor_meta.copy()
    donor_out.insert(0, "cohort", cohort)
    donor_out.insert(0, "analysis", analysis)
    donor_out["pathway_score"] = score

    measured_out = pd.DataFrame({
        "analysis": analysis,
        "cohort": cohort,
        "gene": genes,
        "measured": measured,
        "donors_with_at_least_10_counts": (donor_counts >= 10).sum(axis=0),
        "total_raw_count": donor_counts.sum(axis=0),
    })

    gene_rows = []
    for j, gene in enumerate(genes):
        if not measured[j]:
            continue
        result = fit_ols(donor_meta, logcpm[:, j], covariates)
        result.update({"analysis": analysis, "cohort": cohort, "gene": gene})
        gene_rows.append(result)
    gene_out = pd.DataFrame(gene_rows)
    gene_out["FDR_BH"] = stats.false_discovery_control(gene_out["p_two_sided"].to_numpy(), method="bh")
    return model, donor_out, measured_out, gene_out


def leave_one_out(
    analysis: str,
    cohort: str,
    donor_meta: pd.DataFrame,
    donor_counts: np.ndarray,
    genes: list[str],
    covariates: list[str],
) -> pd.DataFrame:
    measured = (donor_counts >= 10).sum(axis=0) >= 3
    logcpm = np.log2(((donor_counts + 0.5) / (donor_meta["total_umi"].to_numpy()[:, None] + 1)) * 1e6)
    rows = []
    for excluded in np.where(measured)[0]:
        use = measured.copy()
        use[excluded] = False
        result = fit_ols(donor_meta, logcpm[:, use].mean(axis=1), covariates)
        result.update({"analysis": analysis, "cohort": cohort, "excluded_gene": genes[excluded]})
        rows.append(result)
    return pd.DataFrame(rows)


def main() -> None:
    genes = target_genes()
    for filename, expected in GSE243639_SHA256.items():
        observed = sha256(RAW / "gse243639" / filename)
        if observed != expected:
            raise RuntimeError(f"GSE243639 input failed SHA256: {filename}")

    kamath_meta, kamath_counts, kamath_umi = load_kamath(genes)
    mart_meta, _ = load_gse243639_cells()
    mart_counts, mart_umi = load_gse243639_counts(mart_meta, genes)

    kamath_pure = kamath_meta["status"].map({"Ctrl": 0.0, "PD": 1.0})
    kamath_disease = kamath_meta["status"].map({"Ctrl": 0.0, "PD": 1.0, "LBD": 1.0})
    mart_condition = mart_meta["status"].map({"Control": 0.0, "Parkinson's": 1.0})

    selections: dict[str, tuple[str, pd.DataFrame, np.ndarray, np.ndarray, pd.Series, int, np.ndarray | None]] = {
        "GSE178265_DA_PD_vs_control_min5": (
            "GSE178265", kamath_meta, kamath_counts, kamath_umi, kamath_pure, 5, None
        ),
        "GSE178265_DA_PD_plus_LBD_vs_control_min5": (
            "GSE178265", kamath_meta, kamath_counts, kamath_umi, kamath_disease, 5, None
        ),
        "GSE243639_DA_PD_vs_control_min5": (
            "GSE243639", mart_meta, mart_counts, mart_umi, mart_condition, 5,
            (mart_meta["subtype"] == "neurons00").to_numpy(),
        ),
        "GSE243639_DA_PD_vs_control_min1": (
            "GSE243639", mart_meta, mart_counts, mart_umi, mart_condition, 1,
            (mart_meta["subtype"] == "neurons00").to_numpy(),
        ),
        "GSE243639_all_neurons_PD_vs_control_min5": (
            "GSE243639", mart_meta, mart_counts, mart_umi, mart_condition, 5, None
        ),
    }
    aggregated: dict[str, tuple[pd.DataFrame, np.ndarray, str]] = {}
    cell_audit = []
    for name, (cohort, meta, counts, umi, condition, minimum, selection) in selections.items():
        donor_meta, donor_counts = aggregate_donors(meta, counts, umi, condition, minimum, selection)
        aggregated[name] = (donor_meta, donor_counts, cohort)
        for row in donor_meta.itertuples(index=False):
            cell_audit.append({
                "selection": name,
                "cohort": cohort,
                "donor": row.donor,
                "pd": int(row.pd),
                "n_cells": row.n_cells,
                "total_umi": row.total_umi,
            })

    definitions = [
        ("GSE178265_primary_adjusted", "GSE178265_DA_PD_vs_control_min5", ["age", "sex", "pmi"], True),
        ("GSE178265_primary_unadjusted", "GSE178265_DA_PD_vs_control_min5", [], False),
        ("GSE243639_replication_adjusted", "GSE243639_DA_PD_vs_control_min5", ["age", "sex", "pmi"], True),
        ("GSE243639_replication_unadjusted", "GSE243639_DA_PD_vs_control_min5", [], False),
        ("GSE178265_PD_plus_LBD_adjusted", "GSE178265_DA_PD_plus_LBD_vs_control_min5", ["age", "sex", "pmi"], False),
        ("GSE243639_DA_min1_adjusted", "GSE243639_DA_PD_vs_control_min1", ["age", "sex", "pmi"], False),
        ("GSE243639_all_neurons_adjusted", "GSE243639_all_neurons_PD_vs_control_min5", ["age", "sex", "pmi"], False),
        ("GSE243639_replication_RIN_adjusted", "GSE243639_DA_PD_vs_control_min5", ["age", "sex", "pmi", "rin"], False),
    ]
    models, donors, measured, gene_results = [], [], [], []
    for analysis, selection_name, covariates, confirmatory in definitions:
        donor_meta, donor_counts, cohort = aggregated[selection_name]
        model, donor_out, measured_out, gene_out = analyze(
            analysis, cohort, donor_meta, donor_counts, genes, covariates, confirmatory
        )
        models.append(model)
        donors.append(donor_out)
        measured.append(measured_out)
        gene_results.append(gene_out)

    subtype_rows = []
    for subtype in sorted(kamath_meta["subtype"].unique()):
        selection = (kamath_meta["subtype"] == subtype).to_numpy()
        donor_meta, donor_counts = aggregate_donors(
            kamath_meta, kamath_counts, kamath_umi, kamath_pure, 5, selection
        )
        if (donor_meta["pd"] == 1).sum() < 3 or (donor_meta["pd"] == 0).sum() < 3:
            subtype_rows.append({
                "subtype": subtype,
                "n_PD": int((donor_meta["pd"] == 1).sum()),
                "n_control": int((donor_meta["pd"] == 0).sum()),
                "outcome": "NOT_EVALUABLE",
            })
            continue
        model, _, _, _ = analyze(
            f"GSE178265_subtype_{subtype}", "GSE178265", donor_meta, donor_counts,
            genes, ["age", "sex", "pmi"], False,
        )
        model["subtype"] = subtype
        subtype_rows.append(model)

    loo = []
    for analysis, selection_name in (
        ("GSE178265_primary_adjusted", "GSE178265_DA_PD_vs_control_min5"),
        ("GSE243639_replication_adjusted", "GSE243639_DA_PD_vs_control_min5"),
    ):
        donor_meta, donor_counts, cohort = aggregated[selection_name]
        loo.append(leave_one_out(analysis, cohort, donor_meta, donor_counts, genes, ["age", "sex", "pmi"]))

    shared_models = []
    shared_indices = np.asarray([genes.index(x) for x in SHARED8])
    for analysis, selection_name in (
        ("GSE178265_shared8_descriptive", "GSE178265_DA_PD_vs_control_min5"),
        ("GSE243639_shared8_descriptive", "GSE243639_DA_PD_vs_control_min5"),
    ):
        donor_meta, donor_counts, cohort = aggregated[selection_name]
        logcpm = np.log2(((donor_counts[:, shared_indices] + 0.5) /
                          (donor_meta["total_umi"].to_numpy()[:, None] + 1)) * 1e6)
        result = fit_ols(donor_meta, logcpm.mean(axis=1), ["age", "sex", "pmi"])
        result.update({"analysis": analysis, "cohort": cohort, "genes_n": len(SHARED8)})
        shared_models.append(result)

    model_frame = pd.DataFrame(models)
    primary = model_frame.set_index("analysis").loc["GSE178265_primary_adjusted"]
    replication = model_frame.set_index("analysis").loc["GSE243639_replication_adjusted"]
    if primary["outcome"] == "CONFIRMED" and replication["outcome"] == "CONFIRMED":
        overall = "INDEPENDENTLY_REPLICATED"
    elif primary["outcome"] == "CONFIRMED":
        overall = "PRIMARY_ONLY"
    elif primary["outcome"] == "PARTIAL_SUPPORT":
        overall = "PRIMARY_PARTIAL_SUPPORT"
    else:
        overall = "NOT_CONFIRMED"
    summary = pd.DataFrame([{
        "mechanism": "Reactome Proteasome assembly (R-HSA-9907900)",
        "declared_direction": "Down in PD",
        "primary_outcome": primary["outcome"],
        "replication_outcome": replication["outcome"],
        "overall_outcome": overall,
        "primary_beta": primary["beta_PD"],
        "primary_one_sided_p": primary["p_one_sided_down"],
        "replication_beta": replication["beta_PD"],
        "replication_one_sided_p": replication["p_one_sided_down"],
    }])

    write_tsv(pd.DataFrame(cell_audit), "donor_cell_count_audit.tsv")
    write_tsv(model_frame, "pathway_model_results.tsv")
    write_tsv(pd.concat(donors, ignore_index=True), "donor_pathway_scores.tsv")
    write_tsv(pd.concat(measured, ignore_index=True), "measured_gene_audit.tsv")
    write_tsv(pd.concat(gene_results, ignore_index=True), "gene_level_results.tsv")
    write_tsv(pd.DataFrame(subtype_rows), "GSE178265_subtype_sensitivity.tsv")
    write_tsv(pd.concat(loo, ignore_index=True), "leave_one_gene_out.tsv")
    write_tsv(pd.DataFrame(shared_models), "shared8_descriptive_results.tsv")
    write_tsv(summary, "validation_summary.tsv")

    manifest_rows = [
        {"input": str(GMT.relative_to(ROOT)), "sha256": sha256(GMT)},
        {"input": str(KAMATH_CELL_META.relative_to(ROOT)), "sha256": sha256(KAMATH_CELL_META)},
    ]
    for path in sorted((RAW / "gse243639").glob("GSE243639_*")):
        if path.is_file() and path.suffix != ".npz":
            manifest_rows.append({"input": str(path.relative_to(ROOT)), "sha256": sha256(path)})
    cluster_path = RAW / "gse178265" / "SCP1768_DA_cluster.json"
    manifest_rows.append({"input": str(cluster_path.relative_to(ROOT)), "sha256": sha256(cluster_path)})
    gene_hash = hashlib.sha256()
    for path in sorted((RAW / "gse178265" / "genes").glob("*.json")):
        gene_hash.update(path.name.encode())
        gene_hash.update(bytes.fromhex(sha256(path)))
    manifest_rows.append({"input": "data_raw/human_sn_proteasome/gse178265/genes/*.json", "sha256": gene_hash.hexdigest()})
    write_tsv(pd.DataFrame(manifest_rows), "input_manifest.tsv")

    print(summary.to_string(index=False))
    print(model_frame[["analysis", "beta_PD", "p_one_sided_down", "outcome"]].to_string(index=False))


if __name__ == "__main__":
    main()
