repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
scripts_dir <- file.path(repo_root, "scripts")

if (!dir.exists(scripts_dir)) {
  stop("run_pipeline.R must be run from the repository root.")
}

run_fetch_geo <- identical(tolower(Sys.getenv("RUN_FETCH_GEO", "false")), "true")
run_optional_meta <- identical(tolower(Sys.getenv("RUN_OPTIONAL_GSE187012_META", "false")), "true")

download_scripts <- c(
  "scripts/01_fetch_geo.R"
)

optional_scripts <- c(
  "scripts/03a_fetch_gse187012_meta.R"
)

core_scripts <- c(
  "scripts/02_clean_gse17542_microarray.R",
  "scripts/03_build_gse187012_signatures.R",
  "scripts/03c_build_and_save_gse187012_seurat.R",
  "scripts/04_deg_and_perturbation.R",
  "scripts/05_deconvolution_gse17542.R",
  "scripts/06_integrate_rank_visualize.R",
  "scripts/06a_extract_original20.R",
  "scripts/07_figures_original20.R",
  "scripts/08_celltype_specific_DE_original20.R",
  "scripts/10_GSEA_ranked_lists.R",
  "scripts/11_fig4_celltype_heatmap_and_dotplot.R",
  "scripts/12_compare_candidate20_vs_top20.R",
  "scripts/13_permutation_matched_original20.R",
  "scripts/13b_plot_permutation_matched_original20.R",
  "scripts/14_plot_gsea_results.R",
  "scripts/15_results_pack_original20.R"
)

scripts_to_run <- c(
  if (run_fetch_geo) download_scripts,
  "scripts/02_clean_gse17542_microarray.R",
  if (run_optional_meta) optional_scripts,
  core_scripts[core_scripts != "scripts/02_clean_gse17542_microarray.R"]
)

missing_scripts <- scripts_to_run[!file.exists(file.path(repo_root, scripts_to_run))]
if (length(missing_scripts) > 0) {
  stop("Missing script(s): ", paste(missing_scripts, collapse = ", "))
}

dir.create(file.path(repo_root, "metadata"), showWarnings = FALSE)
dir.create(file.path(repo_root, "results"), showWarnings = FALSE)

cat("Repository root:", repo_root, "\n")
cat("RUN_FETCH_GEO =", run_fetch_geo, "\n")
cat("RUN_OPTIONAL_GSE187012_META =", run_optional_meta, "\n")
cat("Scripts queued:", length(scripts_to_run), "\n\n")

for (script in scripts_to_run) {
  cat("=== Running", script, "===\n")
  source(file.path(repo_root, script), echo = TRUE, chdir = FALSE)
  cat("=== Completed", script, "===\n\n")
}

cat("Pipeline completed successfully.\n")

