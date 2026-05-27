# pd-envtox

`pd-envtox` is an R-based analysis pipeline for integrating bulk microarray and single-nucleus RNA-seq data to rank environmentally responsive Parkinson's disease candidate genes.

The current local folder can remain named `terra_pd_extension` if you want. The GitHub repository name should be `pd-envtox`.

## Datasets

- `GSE17542`: bulk microarray MPTP mouse substantia nigra dataset
- `GSE187012`: single-nucleus RNA-seq paraquat/maneb exposure dataset

## Repository Layout

- `scripts/`: numbered analysis scripts
- `data_raw/`: downloaded or extracted raw inputs used directly by the pipeline
- `data_intermediate/`: regenerated intermediate files; do not track in Git
- `results/`: final tables and figures intended to be frozen for a release
- `metadata/`: reproducibility metadata such as session info and checksums

## Reproducible Setup

Before the first public push, initialize `renv` from the repository root and commit the resulting `renv.lock`:

```r
install.packages("renv")
renv::init(bare = TRUE)

install.packages(c(
  "data.table",
  "Matrix",
  "Seurat",
  "nnls",
  "ggplot2",
  "ggrepel",
  "msigdbr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GEOquery", "limma", "clusterProfiler"))

renv::snapshot()
```

Anyone reproducing the project should then be able to run:

```r
renv::restore()
```

## Running The Pipeline

From the repository root:

```powershell
Rscript run_pipeline.R
```

By default, `run_pipeline.R` uses the local inputs already present in the repository and only falls back to `scripts/01_fetch_geo.R` when required GEO seed inputs are missing.

If the required GEO seed inputs are missing, `run_pipeline.R` now auto-enables `scripts/01_fetch_geo.R` so a fresh clone can bootstrap itself more reliably.

To include live GEO downloads:

```powershell
$env:RUN_FETCH_GEO = "true"
Rscript run_pipeline.R
```

`RUN_FETCH_GEO=true` forces a refresh pass even when the seed inputs are already present.

To also fetch and save the optional `GSE187012` series matrix metadata:

```powershell
$env:RUN_FETCH_GEO = "true"
$env:RUN_OPTIONAL_GSE187012_META = "true"
Rscript run_pipeline.R
```

## Script Order

Core execution order:

1. `scripts/01_fetch_geo.R` when raw GEO downloads are needed
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

Currently excluded from the documented pipeline:

- `scripts/09_permutation_test_candidate20.R`

That file is now a deprecated guard script so older notes do not silently run stale logic.

## Versioning Policy

Recommended Git policy for this project:

- track `scripts/`
- track `data_raw/` inputs that are actually consumed by the scripts
- do not track `data_intermediate/`
- track release-grade `results/`
- store `sessionInfo.txt` and checksum manifests in `metadata/`

## Reproducibility Metadata

Suggested metadata files:

- `metadata/sessionInfo.txt`
- `metadata/raw_checksums.csv`
- `metadata/results_checksums.csv`

Instructions are in [metadata/README.md](metadata/README.md).

## Publishing

Use [PUBLISHING_CHECKLIST.md](PUBLISHING_CHECKLIST.md) for the exact Git and GitHub procedure, and [CODE_READABILITY_NOTES.md](CODE_READABILITY_NOTES.md) for code cleanup items worth doing before the first public release.
