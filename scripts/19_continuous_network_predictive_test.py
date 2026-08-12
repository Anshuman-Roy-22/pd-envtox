"""Secondary continuous evaluation of network proximity as response predictor."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
NETWORK = ROOT.parent / "pd-netprox"
OUT = ROOT / "results" / "v2" / "continuous_network_test"
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = {"SNCA", "LRRK2", "PRKN", "PINK1", "PARK7", "VPS35", "ATP13A2",
         "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1", "CHCHD2", "RAB32", "VPS13C"}


def network_table() -> pd.DataFrame:
    ranked = pd.read_csv(NETWORK / "results" / "ranked_candidates.tsv", sep="\t")
    edges = pd.read_csv(NETWORK / "data" / "string" / "ppi_edges.csv", usecols=["protein1", "protein2"])
    degree = pd.concat([edges["protein1"], edges["protein2"]]).value_counts().rename("string_degree")
    ranked = ranked.join(degree, on="ensp_id")
    ranked["gene_symbol"] = ranked["gene_symbol"].str.upper()
    ranked = (ranked.sort_values(["gene_symbol", "proximity", "sn_tpm", "ensp_id"],
                                 ascending=[True, True, False, True])
              .drop_duplicates("gene_symbol"))
    return ranked.loc[~ranked["gene_symbol"].isin(SEEDS)].copy()


def outcomes() -> dict[str, pd.DataFrame]:
    g46798 = pd.read_csv(ROOT / "results/v2/gse46798/GSE46798_factorial_all_genes.tsv.gz", sep="\t")
    g17542 = pd.read_csv(ROOT / "results/v2/gse17542_sn_vta/GSE17542_SN_VTA_all_genes.tsv.gz", sep="\t")
    human = pd.read_csv(ROOT / "results/v2/human_validation/human_PD_directional_meta_all_genes.tsv", sep="\t")
    return {
        "GSE46798_exposure_corrected": g46798.loc[g46798["contrast"] == "exposure_corrected", ["gene", "t"]].assign(response=lambda x: x.t.abs()),
        "GSE17542_SN_MPTP10": g17542.loc[g17542["contrast"] == "SN_MPTP10_vs_control", ["gene", "t"]].assign(response=lambda x: x.t.abs()),
        "GSE17542_interaction10": g17542.loc[g17542["contrast"] == "interaction10_SNminusVTA", ["gene", "t"]].assign(response=lambda x: x.t.abs()),
        "human_PD_postmortem_meta": human[["gene", "meta_z"]].assign(response=lambda x: x.meta_z.abs()),
    }


def partial_spearman(frame: pd.DataFrame) -> tuple[float, float]:
    y = stats.rankdata(frame["response"])
    x = stats.rankdata(-frame["proximity"])
    cov = np.column_stack([
        np.ones(len(frame)),
        stats.rankdata(np.log1p(frame["string_degree"])),
        stats.rankdata(np.log1p(frame["sn_tpm"])),
    ])
    xr = x - cov @ np.linalg.lstsq(cov, x, rcond=None)[0]
    yr = y - cov @ np.linalg.lstsq(cov, y, rcond=None)[0]
    return stats.pearsonr(xr, yr)


def main() -> None:
    network = network_table()
    network.to_csv(OUT / "network_v2_predictor_manifest.tsv", sep="\t", index=False)
    summaries, merged_tables = [], []
    for name, outcome in outcomes().items():
        outcome = outcome.dropna(subset=["response"]).copy()
        outcome["gene"] = outcome["gene"].str.upper()
        frame = network.merge(outcome[["gene", "response"]], left_on="gene_symbol", right_on="gene", validate="one_to_one")
        frame = frame.dropna(subset=["proximity", "sn_tpm", "string_degree", "response"])
        raw = stats.spearmanr(-frame["proximity"], frame["response"])
        partial_r, partial_p = partial_spearman(frame)
        summaries.append({"outcome": name, "n_genes": len(frame), "spearman_r_closeness_vs_abs_response": raw.statistic,
                          "spearman_p": raw.pvalue, "partial_spearman_r": partial_r,
                          "partial_spearman_p": partial_p})
        frame.insert(0, "outcome", name)
        merged_tables.append(frame)
    summary = pd.DataFrame(summaries)
    summary["partial_spearman_FDR_4_tests"] = multipletests(summary["partial_spearman_p"], method="fdr_bh")[1]
    summary.to_csv(OUT / "continuous_network_predictive_summary.tsv", sep="\t", index=False)
    pd.concat(merged_tables).to_csv(OUT / "continuous_network_merged_data.tsv.gz", sep="\t", index=False,
                                           compression={"method": "gzip", "mtime": 0})


if __name__ == "__main__":
    main()
