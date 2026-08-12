"""Outcome-blind validation of the frozen network-v2 panel in human PD SN cohorts."""

from __future__ import annotations

import gzip
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw" / "human_validation"
OUT = ROOT / "results" / "v2" / "human_validation"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = [
    "SNX3", "VPS26B", "VPS26A", "VPS29", "SNX2", "HSPA8", "SNX1",
    "MFN2", "CLTC", "SNCAIP", "EPS15", "UBC", "UBB", "TOMM20",
    "UBA52", "MFN1", "ITSN1", "RHOT1", "CLTCL1", "SH3GL2", "WLS",
]


def read_matrix(path: Path) -> tuple[pd.DataFrame, list[str], list[str]]:
    titles = accessions = None
    rows: list[list[str]] = []
    inside = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!Sample_geo_accession"):
                accessions = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!series_matrix_table_begin"):
                inside = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif inside:
                rows.append([x.strip('"') for x in line.rstrip("\n").split("\t")])
    if titles is None or accessions is None:
        raise RuntimeError(f"Missing sample metadata in {path}")
    header, body = rows[0], rows[1:]
    frame = pd.DataFrame(body, columns=header).set_index("ID_REF").apply(pd.to_numeric, errors="coerce")
    frame = frame.dropna(axis=0, how="any")
    # GEOquery's standard heuristic: these GPL570 matrices are linear-scale.
    if np.nanquantile(frame.to_numpy(), 0.99) > 100:
        frame = np.log2(frame.clip(lower=0) + 1)
    return frame, titles, accessions


def read_annotation(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        for i, line in enumerate(handle):
            if line.startswith("ID\t"):
                skip = i
                break
        else:
            raise RuntimeError("GPL annotation header not found")
    ann = pd.read_csv(path, sep="\t", compression="gzip", skiprows=skip, dtype=str)
    ann = ann.rename(columns={"ID": "probe", "Gene symbol": "gene"})
    ann = ann[["probe", "gene"]].dropna()
    # Freeze to unambiguous, single-symbol probes before examining outcomes.
    ann = ann.loc[~ann["gene"].str.contains(r"///|//|;", regex=True)].copy()
    ann["gene"] = ann["gene"].str.strip().str.upper()
    return ann.loc[ann["gene"].str.fullmatch(r"[A-Z0-9_.-]+", na=False)]


def analyze(accession: str, matrix: pd.DataFrame, titles: list[str], ann: pd.DataFrame) -> pd.DataFrame:
    pd_mask = np.array([("PD-" in x) if accession == "GSE20141" else (" PD " in f" {x} ") for x in titles])
    if accession == "GSE7621":
        pd_mask = np.array([" PD rep" in x for x in titles])
    probes = ann.loc[ann["probe"].isin(matrix.index)].copy()
    baseline = matrix.mean(axis=1).rename("mean_expression")
    probes = probes.join(baseline, on="probe")
    # One probe per gene: highest mean expression across all samples (outcome-blind).
    probes = probes.sort_values(["gene", "mean_expression", "probe"], ascending=[True, False, True])
    chosen = probes.drop_duplicates("gene", keep="first")
    records = []
    for row in chosen.itertuples(index=False):
        y = matrix.loc[row.probe].to_numpy(float)
        case, ctrl = y[pd_mask], y[~pd_mask]
        test = stats.ttest_ind(case, ctrl, equal_var=True)
        records.append({
            "gene": row.gene, "probe": row.probe, "logFC_PD_minus_control": case.mean() - ctrl.mean(),
            "mean_expression": row.mean_expression, "t": test.statistic, "p_value": test.pvalue,
            "n_PD": len(case), "n_control": len(ctrl),
        })
    out = pd.DataFrame(records)
    out["FDR"] = multipletests(out["p_value"], method="fdr_bh")[1]
    out.insert(0, "accession", accession)
    return out.sort_values("p_value")


def main() -> None:
    ann = read_annotation(RAW / "GPL570.annot.gz")
    cohorts = []
    for accession in ("GSE20141", "GSE7621"):
        matrix, titles, _ = read_matrix(RAW / f"{accession}_series_matrix.txt.gz")
        result = analyze(accession, matrix, titles, ann)
        result.to_csv(OUT / f"{accession}_gene_effects.tsv", sep="\t", index=False)
        cohorts.append(result)

    wide = cohorts[0].merge(cohorts[1], on="gene", suffixes=("_20141", "_7621"))
    for suffix in ("20141", "7621"):
        p = wide[f"p_value_{suffix}"].clip(np.finfo(float).tiny, 1)
        wide[f"signed_z_{suffix}"] = np.sign(wide[f"logFC_PD_minus_control_{suffix}"]) * stats.norm.isf(p / 2)
    wide["meta_z"] = (wide["signed_z_20141"] + wide["signed_z_7621"]) / np.sqrt(2)
    wide["meta_p"] = 2 * stats.norm.sf(np.abs(wide["meta_z"]))
    wide["meta_FDR"] = multipletests(wide["meta_p"], method="fdr_bh")[1]
    wide["direction_agreement"] = np.sign(wide["logFC_PD_minus_control_20141"]) == np.sign(wide["logFC_PD_minus_control_7621"])
    wide.to_csv(OUT / "human_PD_directional_meta_all_genes.tsv", sep="\t", index=False)

    panel = pd.DataFrame({"gene": PANEL}).merge(wide, on="gene", how="left", validate="one_to_one")
    panel.to_csv(OUT / "frozen21_human_PD_validation.tsv", sep="\t", index=False)

    summary = pd.DataFrame([{
        "panel_size": len(PANEL),
        "evaluable_both_cohorts": int(panel["meta_p"].notna().sum()),
        "direction_agree_n": int(panel["direction_agreement"].fillna(False).sum()),
        "cohort_FDR_sig_n_GSE20141": int((panel["FDR_20141"] <= 0.05).fillna(False).sum()),
        "cohort_FDR_sig_n_GSE7621": int((panel["FDR_7621"] <= 0.05).fillna(False).sum()),
        "meta_FDR_sig_n": int((panel["meta_FDR"] <= 0.05).fillna(False).sum()),
    }])
    summary.to_csv(OUT / "frozen21_human_PD_summary.tsv", sep="\t", index=False)

    # Matched-background panel test. Match jointly on mean-expression quintiles
    # in both cohorts; test prespecified consistency and effect metrics.
    universe = wide.dropna(subset=["meta_z", "mean_expression_20141", "mean_expression_7621"]).copy()
    for suffix in ("20141", "7621"):
        universe[f"expr_bin_{suffix}"] = pd.qcut(
            universe[f"mean_expression_{suffix}"], 5, labels=False, duplicates="drop"
        )
    universe["stratum"] = universe["expr_bin_20141"].astype(str) + "_" + universe["expr_bin_7621"].astype(str)
    evaluable_panel = universe.loc[universe["gene"].isin(PANEL)].copy()
    background = universe.loc[~universe["gene"].isin(PANEL)].copy()
    pools = {key: grp.index.to_numpy() for key, grp in background.groupby("stratum")}
    rng = np.random.default_rng(20260812)
    null_rows = []
    for draw in range(10_000):
        picked = [rng.choice(pools[s]) for s in evaluable_panel["stratum"]]
        sample = background.loc[picked]
        null_rows.append({
            "draw": draw + 1,
            "direction_agree_n": int(sample["direction_agreement"].sum()),
            "mean_abs_meta_z": float(sample["meta_z"].abs().mean()),
            "mean_meta_z": float(sample["meta_z"].mean()),
        })
    null = pd.DataFrame(null_rows)
    null.to_csv(OUT / "frozen21_matched_null_10000.tsv.gz", sep="\t", index=False, compression="gzip")
    observed = {
        "evaluable_n": len(evaluable_panel),
        "direction_agree_n": int(evaluable_panel["direction_agreement"].sum()),
        "mean_abs_meta_z": float(evaluable_panel["meta_z"].abs().mean()),
        "mean_meta_z": float(evaluable_panel["meta_z"].mean()),
    }
    observed.update({
        "empirical_p_direction_agreement_upper": (1 + (null["direction_agree_n"] >= observed["direction_agree_n"]).sum()) / (len(null) + 1),
        "empirical_p_mean_abs_meta_z_upper": (1 + (null["mean_abs_meta_z"] >= observed["mean_abs_meta_z"]).sum()) / (len(null) + 1),
        "empirical_p_mean_meta_z_two_sided": (1 + (null["mean_meta_z"].abs() >= abs(observed["mean_meta_z"])).sum()) / (len(null) + 1),
    })
    pd.DataFrame([observed]).to_csv(OUT / "frozen21_matched_null_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
