"""Prespecified factorial and module analysis for GSE46798."""

from __future__ import annotations

import gzip
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw" / "gse46798"
OUT = ROOT / "results" / "v2" / "gse46798"
OUT.mkdir(parents=True, exist_ok=True)
RNG_SEED = 20260812

PANEL = ["SNX3", "VPS26B", "VPS26A", "VPS29", "SNX2", "HSPA8", "SNX1",
         "MFN2", "CLTC", "SNCAIP", "EPS15", "UBC", "UBB", "TOMM20",
         "UBA52", "MFN1", "ITSN1", "RHOT1", "CLTCL1", "SH3GL2", "WLS"]
MODULES = {
    "retromer_endosomal": ["SNX1", "SNX2", "SNX3", "VPS26A", "VPS26B", "VPS29"],
    "mitochondrial_dynamics_transport": ["MFN1", "MFN2", "RHOT1", "TOMM20"],
    "endocytosis_vesicle_trafficking": ["CLTC", "CLTCL1", "EPS15", "ITSN1", "SH3GL2"],
    "ubiquitin_proteostasis": ["UBB", "UBC", "UBA52", "HSPA8"],
    "full_frozen21": PANEL,
}
CONTRASTS = {
    "exposure_corrected": np.array([0, 0, 1, 0]),
    "exposure_A53T": np.array([0, 0, 1, 1]),
    "genotype_by_exposure": np.array([0, 0, 0, 1]),
    "A53T_vehicle": np.array([0, 1, 0, 0]),
}


def write_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    payload = frame.to_csv(sep="\t", index=False).encode("utf-8")
    path.write_bytes(gzip.compress(payload, compresslevel=9, mtime=0))
    with gzip.open(path, "rb") as handle:
        if handle.read() != payload:
            raise RuntimeError(f"Compressed output verification failed: {path}")


def read_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    titles, rows, inside = None, [], False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = re.findall(r'"([^"]*)"', line)
            elif line.startswith("!series_matrix_table_begin"):
                inside = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif inside:
                rows.append([x.strip('"') for x in line.rstrip().split("\t")])
    if titles is None:
        raise RuntimeError("Sample titles missing")
    frame = pd.DataFrame(rows[1:], columns=rows[0]).set_index("ID_REF").apply(pd.to_numeric, errors="coerce")
    frame = frame.dropna(how="any")
    if np.nanquantile(frame.to_numpy(), 0.99) > 100:
        frame = np.log2(frame.clip(lower=0) + 1)
    return frame, titles


def read_annotation(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        for i, line in enumerate(handle):
            if line.startswith("ID\t"):
                skip = i
                break
        else:
            raise RuntimeError("Annotation header missing")
    ann = pd.read_csv(path, sep="\t", compression="gzip", skiprows=skip, dtype=str)
    ann = ann.rename(columns={"ID": "probe", "Gene symbol": "gene"})[["probe", "gene"]].dropna()
    ann = ann.loc[~ann["gene"].str.contains(r"///|//|;", regex=True)].copy()
    ann["gene"] = ann["gene"].str.strip().str.upper()
    return ann.loc[ann["gene"].str.fullmatch(r"[A-Z0-9_.-]+", na=False)]


def gene_matrix(matrix: pd.DataFrame, ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = ann.loc[ann["probe"].isin(matrix.index)].copy()
    eligible = eligible.join(matrix.mean(axis=1).rename("mean_expression"), on="probe")
    chosen = (eligible.sort_values(["gene", "mean_expression", "probe"], ascending=[True, False, True])
              .drop_duplicates("gene", keep="first"))
    expression = matrix.loc[chosen["probe"]].copy()
    expression.index = chosen["gene"].to_numpy()
    return expression, chosen


def fit_factorial(expression: pd.DataFrame, titles: list[str]) -> pd.DataFrame:
    genotype = np.array([int("A53T" in x) for x in titles], float)
    exposure = np.array([int("PQ/MB" in x) for x in titles], float)
    design = np.column_stack([np.ones(len(titles)), genotype, exposure, genotype * exposure])
    inv = np.linalg.inv(design.T @ design)
    beta = inv @ design.T @ expression.to_numpy().T
    residual = expression.to_numpy().T - design @ beta
    df = len(titles) - design.shape[1]
    sigma2 = np.sum(residual ** 2, axis=0) / df
    base = pd.DataFrame({"gene": expression.index, "mean_expression": expression.mean(axis=1).to_numpy()})
    frames = []
    for name, contrast in CONTRASTS.items():
        estimate = contrast @ beta
        se = np.sqrt(sigma2 * (contrast @ inv @ contrast))
        t = np.divide(estimate, se, out=np.full_like(estimate, np.nan), where=se > 0)
        p = 2 * stats.t.sf(np.abs(t), df)
        out = base.copy()
        out["contrast"] = name
        out["logFC"] = estimate
        out["t"] = t
        out["p_value"] = p
        out["FDR"] = np.nan
        finite = np.isfinite(p)
        out.loc[finite, "FDR"] = multipletests(p[finite], method="fdr_bh")[1]
        out["residual_df"] = df
        frames.append(out)
    return pd.concat(frames, ignore_index=True)


def matched_tests(results: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RNG_SEED)
    summaries, draws = [], []
    for contrast, universe in results.groupby("contrast", sort=False):
        universe = universe.dropna(subset=["t", "mean_expression"]).copy()
        universe["expr_bin"] = pd.qcut(universe["mean_expression"], 10, labels=False, duplicates="drop")
        for module, genes in MODULES.items():
            observed = universe.loc[universe["gene"].isin(genes)].copy()
            background = universe.loc[~universe["gene"].isin(PANEL)].copy()
            pools = {key: grp.index.to_numpy() for key, grp in background.groupby("expr_bin")}
            picked = np.column_stack([
                rng.choice(pools[b], size=10_000, replace=True) for b in observed["expr_bin"]
            ])
            sampled_t = universe["t"].reindex(picked.ravel()).to_numpy().reshape(picked.shape)
            null_abs = np.mean(np.abs(sampled_t), axis=1)
            null_signed = np.mean(sampled_t, axis=1)
            obs_abs, obs_signed = observed["t"].abs().mean(), observed["t"].mean()
            summaries.append({
                "contrast": contrast, "module": module, "declared_n": len(genes), "evaluable_n": len(observed),
                "mean_abs_t": obs_abs, "mean_signed_t": obs_signed,
                "empirical_p_mean_abs_t_upper": (1 + np.sum(null_abs >= obs_abs)) / 10001,
                "empirical_p_mean_signed_t_two_sided": (1 + np.sum(np.abs(null_signed) >= abs(obs_signed))) / 10001,
            })
            draws.extend({"contrast": contrast, "module": module, "draw": i + 1,
                          "mean_abs_t": a, "mean_signed_t": s}
                         for i, (a, s) in enumerate(zip(null_abs, null_signed)))
    summary = pd.DataFrame(summaries)
    summary["FDR_abs_across_20_tests"] = multipletests(summary["empirical_p_mean_abs_t_upper"], method="fdr_bh")[1]
    summary["FDR_signed_across_20_tests"] = multipletests(summary["empirical_p_mean_signed_t_two_sided"], method="fdr_bh")[1]
    return summary, pd.DataFrame(draws)


def main() -> None:
    matrix, titles = read_matrix(RAW / "GSE46798_series_matrix.txt.gz")
    ann = read_annotation(RAW / "GPL10558.annot.gz")
    expression, chosen = gene_matrix(matrix, ann)
    results = fit_factorial(expression, titles)
    write_tsv_gz(results, OUT / "GSE46798_factorial_all_genes.tsv.gz")
    chosen.to_csv(OUT / "GSE46798_frozen_probe_map.tsv", sep="\t", index=False)
    panel = (pd.MultiIndex.from_product([CONTRASTS, PANEL], names=["contrast", "gene"]).to_frame(index=False)
             .merge(results, on=["contrast", "gene"], how="left", validate="one_to_one"))
    panel.to_csv(OUT / "frozen21_GSE46798_contrasts.tsv", sep="\t", index=False)
    panel_probes = ann.loc[ann["gene"].isin(PANEL) & ann["probe"].isin(matrix.index)].copy()
    probe_expression = matrix.loc[panel_probes["probe"]].copy()
    probe_expression.index = [f"{g}|{p}" for g, p in zip(panel_probes["gene"], panel_probes["probe"])]
    probe_results = fit_factorial(probe_expression, titles)
    probe_results[["gene_symbol", "probe"]] = probe_results["gene"].str.split("|", n=1, expand=True)
    probe_results.drop(columns="gene").to_csv(OUT / "frozen21_all_unambiguous_probe_sensitivity.tsv", sep="\t", index=False)
    summary, draws = matched_tests(results)
    summary.to_csv(OUT / "GSE46798_module_tests.tsv", sep="\t", index=False)
    write_tsv_gz(draws, OUT / "GSE46798_module_matched_null_10000.tsv.gz")


if __name__ == "__main__":
    main()
