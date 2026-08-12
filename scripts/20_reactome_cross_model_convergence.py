"""Prespecified signed Reactome convergence across three PD model classes."""

from __future__ import annotations

from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
GMT = ROOT / "data_raw/pathways/ReactomePathways.gmt"
OUT = ROOT / "results/v2/reactome_convergence"
OUT.mkdir(parents=True, exist_ok=True)
TARGETED = [
    "Aerobic respiration and respiratory electron transport", "Respiratory electron transport",
    "Complex I biogenesis", "Mitochondrial biogenesis", "Cellular response to mitochondrial stress",
    "Mitophagy", "PINK1-PRKN Mediated Mitophagy", "Autophagy", "Macroautophagy",
    "Chaperone Mediated Autophagy", "Lysosome Vesicle Biogenesis",
    "Endosomal Sorting Complex Required For Transport (ESCRT)", "Vesicle-mediated transport",
    "Golgi Associated Vesicle Biogenesis", "Protein ubiquitination", "Proteasome assembly",
    "Neurotransmitter release cycle", "Dopamine Neurotransmitter Release Cycle",
    "Innate Immune System", "Cytokine Signaling in Immune system",
]


def rankings() -> dict[str, pd.DataFrame]:
    a = pd.read_csv(ROOT / "results/v2/gse46798/GSE46798_factorial_all_genes.tsv.gz", sep="\t")
    b = pd.read_csv(ROOT / "results/v2/gse17542_sn_vta/GSE17542_SN_VTA_all_genes.tsv.gz", sep="\t")
    c = pd.read_csv(ROOT / "results/v2/human_validation/human_PD_directional_meta_all_genes.tsv", sep="\t")
    return {
        "GSE46798_exposure_corrected": a.loc[a.contrast == "exposure_corrected", ["gene", "t"]].rename(columns={"t": "score"}),
        "GSE17542_SN_MPTP10": b.loc[b.contrast == "SN_MPTP10_vs_control", ["gene", "t"]].rename(columns={"t": "score"}),
        "human_PD_postmortem_meta": c[["gene", "meta_z"]].rename(columns={"meta_z": "score"}),
    }


def main() -> None:
    tables = []
    for dataset, ranking in rankings().items():
        output_path = OUT / f"{dataset}_Reactome_GSEA.tsv"
        if output_path.exists():
            tables.append(pd.read_csv(output_path, sep="\t"))
            continue
        ranking = ranking.dropna().copy()
        ranking["gene"] = ranking["gene"].str.upper()
        ranking = ranking.sort_values(["score", "gene"], ascending=[False, True]).drop_duplicates("gene")
        # Deterministic sub-machine-precision tie break; cannot cross a non-tied score.
        ranking["score"] = ranking["score"] + np.linspace(1e-9, 0, len(ranking), endpoint=False)
        pre = gp.prerank(rnk=ranking, gene_sets=str(GMT), min_size=10, max_size=500,
                         permutation_num=1000, seed=20260812, threads=1, verbose=False, outdir=None)
        result = pre.res2d.rename(columns={"Term": "pathway", "NOM p-val": "nominal_p",
                                           "FDR q-val": "dataset_FDR", "Lead_genes": "leading_genes"})
        result.insert(0, "dataset", dataset)
        result.to_csv(output_path, sep="\t", index=False)
        tables.append(result)

    all_results = pd.concat(tables, ignore_index=True)
    all_results.to_csv(OUT / "all_dataset_Reactome_GSEA.tsv", sep="\t", index=False)
    wide_nes = all_results.pivot(index="pathway", columns="dataset", values="NES").dropna()
    wide_p = all_results.pivot(index="pathway", columns="dataset", values="nominal_p").reindex(wide_nes.index)
    z_columns = []
    for dataset in wide_nes.columns:
        p = pd.to_numeric(wide_p[dataset], errors="coerce").clip(lower=1 / 1001, upper=1)
        z_columns.append(np.sign(pd.to_numeric(wide_nes[dataset])) * stats.norm.isf(p / 2))
    z = pd.concat(z_columns, axis=1)
    z.columns = wide_nes.columns
    meta = pd.DataFrame(index=wide_nes.index)
    meta["meta_z"] = z.sum(axis=1) / np.sqrt(z.shape[1])
    meta["meta_p"] = 2 * stats.norm.sf(meta["meta_z"].abs())
    meta["meta_FDR"] = multipletests(meta["meta_p"], method="fdr_bh")[1]
    meta["same_NES_direction"] = (wide_nes.gt(0).all(axis=1) | wide_nes.lt(0).all(axis=1))
    for dataset in wide_nes.columns:
        meta[f"NES_{dataset}"] = wide_nes[dataset]
        meta[f"nominal_p_{dataset}"] = wide_p[dataset]
    meta = meta.reset_index().sort_values(["meta_FDR", "meta_p", "pathway"])
    meta.to_csv(OUT / "Reactome_cross_model_signed_meta.tsv", sep="\t", index=False)

    targeted = meta.loc[meta["pathway"].isin(TARGETED)].copy()
    targeted["prespecified_order"] = targeted["pathway"].map({x: i for i, x in enumerate(TARGETED)})
    targeted.sort_values("prespecified_order").to_csv(OUT / "Reactome_targeted_PD_family.tsv", sep="\t", index=False)

    p_columns = [c for c in meta if c.startswith("nominal_p_")]
    strict = meta.loc[meta["same_NES_direction"] & meta[p_columns].lt(0.05).all(axis=1)].copy()
    leading = all_results.pivot(index="pathway", columns="dataset", values="leading_genes")
    intersections, unions = [], []
    for pathway in strict["pathway"]:
        gene_sets = [set(str(x).split(";")) for x in leading.loc[pathway].dropna()]
        intersections.append(";".join(sorted(set.intersection(*gene_sets))))
        unions.append(len(set.union(*gene_sets)))
    strict["shared_leading_genes_all_3"] = intersections
    strict["leading_gene_union_n"] = unions
    strict.to_csv(OUT / "Reactome_strict_three_model_recurrence.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
