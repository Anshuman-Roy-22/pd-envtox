# Code Readability And Publication Notes

These are the cleanup items worth addressing before the first public `pd-envtox` release.

## Highest Priority

1. Decide whether to keep or delete `scripts/09_permutation_test_candidate20.R`.

- It is now a deprecated guard script rather than a runnable analysis step.
- Keeping it preserves historical numbering in notes and slide decks.
- Deleting it later is reasonable once no external references depend on it.

2. Keep future 10x input changes centralized in `scripts/helpers_10x_io.R`.

The duplicated helper logic from these scripts has been consolidated:

- `scripts/03_build_gse187012_signatures.R`
- `scripts/03c_build_and_save_gse187012_seurat.R`
- `scripts/04_deg_and_perturbation.R`

Any future fixes to feature parsing, barcode parsing, matrix orientation, or Seurat object creation should be made in the shared helper file first.

3. Replace magic constants with named parameters.

Examples:

- `scripts/06_integrate_rank_visualize.R`: `0.66`, `0.34`, `0.20`
- `scripts/03_build_gse187012_signatures.R`: marker cutoffs and cluster labeling assumptions
- `scripts/04_deg_and_perturbation.R`: pseudocount and scaling choices

Recommended change:

- define the constants near the top of each script
- add a one-line comment explaining the rationale

4. Cache live annotation lookups that affect results.

In `scripts/02_clean_gse17542_microarray.R`, GPL annotation may be fetched live with `getGEO(gpl_id)`.

That is convenient, but it makes exact reproduction depend on a network call and upstream availability. A stronger publication version would:

- save the retrieved GPL table into `data_raw/` or `data_intermediate/`
- record the platform ID in the README
- avoid recomputing the mapping from a live request when a local cache is present

## Medium Priority

5. Clean encoding artifacts before publication.

Examples already visible:

- `scripts/04_deg_and_perturbation.R`
- `scripts/03_build_gse187012_signatures.R`
- `scripts/05_deconvolution_gse17542.R`

The files appear to contain mojibake in comments. Standardize the repo on UTF-8 and resave those scripts cleanly.

6. Add short script headers everywhere.

Each script should open with 4 lines:

- purpose
- required inputs
- outputs written
- key package dependencies

7. Make optional steps explicit.

`scripts/03a_fetch_gse187012_meta.R` looks optional because the current downstream scripts do not consume `gse187012_seriesmatrix.rds`.

Document it as optional in both the README and the script header.

8. Consolidate hard-coded paths.

The current pipeline is already better than average because it uses repo-relative paths and avoids `setwd()`.

Still, readability would improve if common paths were standardized:

- define `raw_dir`, `intermediate_dir`, and `results_dir`
- avoid repeating long path literals across scripts

## Low Priority

9. Capture assumptions in comments where they drive interpretation.

Examples:

- why the ranking weights favor `mptp10d`
- why `shared_loose` uses a threshold of `0.20`
- why the analysis focuses on the original 20-gene panel at later stages

10. Add a single orchestration point and keep it authoritative.

`run_pipeline.R` is now the intended top-level entrypoint. If script order changes, update that file first and treat it as the source of truth.
