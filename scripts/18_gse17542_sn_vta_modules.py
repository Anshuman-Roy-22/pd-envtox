"""Prespecified regional-vulnerability analysis of GSE17542."""

from __future__ import annotations

import gzip
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw" / "gse17542"
OUT = ROOT / "results" / "v2" / "gse17542_sn_vta"
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
# Cell order: SN_C, SN_2d, SN_10d, VTA_C, VTA_2d, VTA_10d.
CONTRASTS = {
    "SN_MPTP10_vs_control": np.array([-1, 0, 1, 0, 0, 0]),
    "interaction10_SNminusVTA": np.array([-1, 0, 1, 1, 0, -1]),
    "VTA_MPTP10_vs_control": np.array([0, 0, 0, -1, 0, 1]),
    "SN_MPTP2_vs_control": np.array([-1, 1, 0, 0, 0, 0]),
    "interaction2_SNminusVTA": np.array([-1, 1, 0, 1, -1, 0]),
    "VTA_MPTP2_vs_control": np.array([0, 0, 0, -1, 1, 0]),
    "baseline_SN_minus_VTA": np.array([1, 0, 0, -1, 0, 0]),
}


def write_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    payload = frame.to_csv(sep="\t", index=False).encode()
    path.write_bytes(gzip.compress(payload, compresslevel=9, mtime=0))
    if gzip.decompress(path.read_bytes()) != payload:
        raise RuntimeError(f"gzip verification failed for {path}")


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
    frame = pd.DataFrame(rows[1:], columns=rows[0]).set_index("ID_REF").apply(pd.to_numeric, errors="coerce").dropna()
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


def collapse_genes(matrix: pd.DataFrame, ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ann = ann.loc[ann["probe"].isin(matrix.index)].drop_duplicates(["probe", "gene"])
    frames, meta = [], []
    for gene, group in ann.groupby("gene", sort=True):
        probes = group["probe"].tolist()
        frames.append(pd.Series(matrix.loc[probes].median(axis=0), name=gene))
        meta.append({"gene": gene, "eligible_probe_n": len(probes), "eligible_probes": ";".join(sorted(probes))})
    return pd.DataFrame(frames), pd.DataFrame(meta)


def fit(expression: pd.DataFrame) -> pd.DataFrame:
    design = np.repeat(np.eye(6), 3, axis=0)
    inv = np.linalg.inv(design.T @ design)
    beta = inv @ design.T @ expression.to_numpy().T
    residual = expression.to_numpy().T - design @ beta
    df = 18 - 6
    sigma2 = np.sum(residual ** 2, axis=0) / df
    frames = []
    for name, contrast in CONTRASTS.items():
        estimate = contrast @ beta
        se = np.sqrt(sigma2 * (contrast @ inv @ contrast))
        t = np.divide(estimate, se, out=np.full_like(estimate, np.nan), where=se > 0)
        p = 2 * stats.t.sf(np.abs(t), df)
        out = pd.DataFrame({"gene": expression.index, "contrast": name, "logFC": estimate,
                            "t": t, "p_value": p, "residual_variance": sigma2, "residual_df": df})
        out["FDR"] = np.nan
        finite = np.isfinite(p)
        out.loc[finite, "FDR"] = multipletests(p[finite], method="fdr_bh")[1]
        frames.append(out)
    return pd.concat(frames, ignore_index=True)


def module_tests(results: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng, summaries, draws = np.random.default_rng(RNG_SEED), [], []
    base_meta = metadata.copy()
    base_meta["probe_category"] = base_meta["eligible_probe_n"].clip(upper=3)
    for contrast, universe in results.groupby("contrast", sort=False):
        universe = universe.dropna(subset=["t", "residual_variance"]).merge(base_meta, on="gene")
        universe["variance_bin"] = pd.qcut(universe["residual_variance"], 10, labels=False, duplicates="drop")
        universe["stratum"] = universe["probe_category"].astype(str) + "_" + universe["variance_bin"].astype(str)
        background = universe.loc[~universe["gene"].isin(PANEL)].copy()
        pools = {key: grp.index.to_numpy() for key, grp in background.groupby("stratum")}
        for module, genes in MODULES.items():
            observed = universe.loc[universe["gene"].isin(genes)]
            missing_pools = sorted(set(observed["stratum"]) - set(pools))
            if missing_pools:
                raise RuntimeError(f"Empty matched strata: {missing_pools}")
            picked = np.column_stack([rng.choice(pools[s], 10_000, replace=True) for s in observed["stratum"]])
            sampled_t = universe["t"].reindex(picked.ravel()).to_numpy().reshape(picked.shape)
            null_abs, null_signed = np.abs(sampled_t).mean(1), sampled_t.mean(1)
            obs_abs, obs_signed = observed["t"].abs().mean(), observed["t"].mean()
            summaries.append({"contrast": contrast, "module": module, "declared_n": len(genes),
                              "evaluable_n": len(observed), "mean_abs_t": obs_abs, "mean_signed_t": obs_signed,
                              "empirical_p_mean_abs_t_upper": (1 + (null_abs >= obs_abs).sum()) / 10001,
                              "empirical_p_mean_signed_t_two_sided": (1 + (np.abs(null_signed) >= abs(obs_signed)).sum()) / 10001})
            draws.extend({"contrast": contrast, "module": module, "draw": i + 1,
                          "mean_abs_t": a, "mean_signed_t": s}
                         for i, (a, s) in enumerate(zip(null_abs, null_signed)))
    summary = pd.DataFrame(summaries)
    summary["FDR_abs_across_35_tests"] = multipletests(summary["empirical_p_mean_abs_t_upper"], method="fdr_bh")[1]
    summary["FDR_signed_across_35_tests"] = multipletests(summary["empirical_p_mean_signed_t_two_sided"], method="fdr_bh")[1]
    return summary, pd.DataFrame(draws)


def main() -> None:
    matrix, _ = read_matrix(RAW / "GSE17542_series_matrix.txt.gz")
    ann = read_annotation(RAW / "GPL1261.annot.gz")
    expression, metadata = collapse_genes(matrix, ann)
    results = fit(expression)
    write_tsv_gz(results, OUT / "GSE17542_SN_VTA_all_genes.tsv.gz")
    metadata.to_csv(OUT / "GSE17542_gene_probe_manifest.tsv", sep="\t", index=False)
    panel = (pd.MultiIndex.from_product([CONTRASTS, PANEL], names=["contrast", "gene"]).to_frame(index=False)
             .merge(results, on=["contrast", "gene"], how="left", validate="one_to_one"))
    panel.to_csv(OUT / "frozen21_GSE17542_SN_VTA_contrasts.tsv", sep="\t", index=False)
    summary, draws = module_tests(results, metadata)
    summary.to_csv(OUT / "GSE17542_SN_VTA_module_tests.tsv", sep="\t", index=False)
    write_tsv_gz(draws, OUT / "GSE17542_SN_VTA_module_null_10000.tsv.gz")


if __name__ == "__main__":
    main()
