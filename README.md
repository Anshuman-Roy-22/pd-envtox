# pd-envtox

R pipeline for integrating bulk microarray and single-nucleus RNA-seq data to rank environmentally responsive Parkinson's disease candidate genes.

## V2 validation status

The V2 branch records both positive and negative results under frozen analysis
plans. The original 21-gene network panel did not validate as a general
toxin-responsive panel. Proteasome assembly showed retrospective cross-model
convergence but did not meet the prespecified mature-LUHMES rotenone replication
threshold. A later held-out test of the nominated
cilium/Hedgehog/microtubule-trafficking mechanism was also `NOT_CONFIRMED`
(GSE196190 one-sided-down p = 0.1221), and two stronger external LUHMES
sensitivities were null.

A final frozen donor-level human substantia nigra test confirmed lower Reactome
Proteasome assembly expression in GSE178265 DA nuclei (adjusted beta = -0.4000,
one-sided p = 0.0352), but the independent GSE243639 DA-neuron replication was
null (beta = +0.0042, p = 0.5095). Its frozen overall label is `PRIMARY_ONLY`,
not independently replicated.

The defensible project conclusion is context-dependent pathway convergence, not
a universally replicated environmental-PD mechanism. See
`docs/v2/human_sn_proteasome_validation_report.md` for the latest formal result.

## Key Components of Repository

- A canonical, numbered R pipeline instead of multiple competing analysis paths.
- Automatic GEO bootstrapping when required seed inputs are missing from a fresh clone.
- Tracked raw snRNA-seq inputs plus tracked release-grade results and figures.
- A project-local `renv` environment with a committed lockfile for reproducible package restoration.
- SHA256 checksum manifests and session metadata for reproducibility-critical inputs and outputs.
- Detailed scripts, algorithms, etc.

## Quick start

```powershell
Rscript -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', repos = 'https://cloud.r-project.org')"
Rscript -e "renv::restore()"
Rscript run_pipeline.R
```

If the required GEO seed inputs are missing, `run_pipeline.R` will auto-enable `scripts/01_fetch_geo.R`.

To force a fresh GEO download pass:

```powershell
$env:RUN_FETCH_GEO = "true"
Rscript run_pipeline.R
```

To also fetch the optional `GSE187012` series-matrix metadata snapshot:

```powershell
$env:RUN_FETCH_GEO = "true"
$env:RUN_OPTIONAL_GSE187012_META = "true"
Rscript run_pipeline.R
```

## Canonical pipeline order

1. `scripts/01_fetch_geo.R` when GEO seed files are missing or a refresh is requested
2. `scripts/02_clean_gse17542_microarray.R`
3. `scripts/03_build_gse187012_signatures.R`
4. `scripts/03c_build_and_save_gse187012_seurat.R`
5. `scripts/04_deg_and_perturbation.R`
6. `scripts/05_deconvolution_gse17542.R`
7. `scripts/06_integrate_rank_visualize.R`
8. `scripts/06a_extract_original20.R`
9. `scripts/07_figures_original20.R`
10. `scripts/08_celltype_specific_DE_original20.R`
11. `scripts/10_GSEA_ranked_lists.R`
12. `scripts/11_fig4_celltype_heatmap_and_dotplot.R`
13. `scripts/12_compare_candidate20_vs_top20.R`
14. `scripts/13_permutation_matched_original20.R`
15. `scripts/13b_plot_permutation_matched_original20.R`
16. `scripts/14_plot_gsea_results.R`
17. `scripts/15_results_pack_original20.R`

Optional:

- `scripts/03a_fetch_gse187012_meta.R`

Not part of the runnable pipeline:

- `scripts/09_permutation_test_candidate20.R`

That file is a deprecated guard script kept only so older notes do not silently run stale logic. Not implementable.

## Tracked inputs and outputs

Tracked source or reduced-input files:

- `data_raw/GSE187012/GSE187012_RAW.tar`
- `data_raw/GSE187012/extracted/GSM5667021_CTL_barcodes.tsv.gz`
- `data_raw/GSE187012/extracted/GSM5667021_CTL_features.tsv.gz`
- `data_raw/GSE187012/extracted/GSM5667021_CTL_matrix.mtx.gz`
- `data_raw/GSE187012/extracted/GSM5667022_EXP_barcodes.tsv.gz`
- `data_raw/GSE187012/extracted/GSM5667022_EXP_features.tsv.gz`
- `data_raw/GSE187012/extracted/GSM5667022_EXP_matrix.mtx.gz`
- `metadata/raw_checksums.csv`

Tracked canonical outputs:

- `results/robust_sporadic_candidates_ranked.csv`
- `results/top200_candidates_for_writeup.csv`
- `results/original20_table1.csv`
- `results/original20_results_pack.csv`
- `results/original20_celltype_specific_logFC.csv`
- `results/deconv_gse17542_proportions.csv`
- `results/gse17542_limma_mptp2d_vs_control.csv`
- `results/gse17542_limma_mptp10d_vs_control.csv`
- `results/gse187012_pseudobulk_log2fc.csv`
- `results/permutation_matched_original20.csv`
- `results/compare_candidate20_vs_top20_metrics.csv`
- `results/fig1_original20_reactivity_bar.png`
- `results/fig2_original20_scatter_mptp10d_vs_pq.png`
- `results/fig3_original20_direction_tile.png`
- `results/fig4A_original20_celltype_logFC_heatmap.png`
- `results/fig4B_original20_dotplot_celltypes.png`
- `results/fig5_original20_norm_heatmap.png`
- `results/fig5_permutation_matched_original20.png`
- `results/fig6_candidate20_vs_top20_boxplots.png`
- `results/fig7A_gsea_mptp10d_reactome.png`
- `results/fig7B_gsea_mptp10d_gobp.png`
- `metadata/results_checksums.csv`

Not tracked in normal Git history:

- `data_intermediate/` rebuildable intermediates
- transient logs such as `results/_*.txt`
- local cache material produced during package restore or exploratory work

## Environment files

- `renv.lock`: Canonical R package lockfile for the repository.
- `.Rprofile` and `renv/activate.R`: Project-local activation hooks for `renv`.
- `metadata/sessionInfo.txt`: Recorded R session metadata for the tracked environment snapshot.
- `pd-envtox.Rproj`: RStudio project file with UTF-8 settings.

## Validation

Canonical rerun:

```powershell
Rscript run_pipeline.R
```

Refresh reproducibility metadata after a rerun:

```powershell
Get-ChildItem data_raw -Recurse -File |
  Get-FileHash -Algorithm SHA256 |
  Export-Csv metadata\raw_checksums.csv -NoTypeInformation

Get-ChildItem results -Recurse -File |
  Get-FileHash -Algorithm SHA256 |
  Export-Csv metadata\results_checksums.csv -NoTypeInformation
```

For environment restoration and metadata notes, see `metadata/README.md`.

## Repository layout

- `scripts/`: Canonical pipeline scripts.
- `data_raw/`: Tracked raw or reduced raw inputs used by the pipeline.
- `data_intermediate/`: Regenerated intermediate objects, not tracked in Git.
- `results/`: Canonical tabular outputs and figures tracked in Git.
- `metadata/`: Checksums, session metadata, and reproducibility notes.
- `renv/`: Project-local reproducible R environment bootstrap files.

## License

This repository is released under the MIT License. See `LICENSE`.

Reach out with questions
