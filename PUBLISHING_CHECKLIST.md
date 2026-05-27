# Publishing Checklist For `pd-envtox`

## Goal

Turn this local folder into a GitHub repository named `pd-envtox` with enough project metadata that another analyst can restore the environment and re-run the workflow.

## Recommended Sequence

1. Keep or rename the local folder.

If you want the local folder name to match the GitHub repository:

```powershell
Set-Location C:\Users\das_t\Documents
Rename-Item terra_pd_extension pd-envtox
Set-Location C:\Users\das_t\Documents\pd-envtox
```

If you keep the current local folder name, that does not block creating a GitHub repository called `pd-envtox`.

2. Initialize and snapshot the R environment.

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

3. Freeze reproducibility metadata.

```r
dir.create("metadata", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "metadata/sessionInfo.txt")
```

```powershell
Get-ChildItem data_raw -Recurse -File |
  Get-FileHash -Algorithm SHA256 |
  Export-Csv metadata\raw_checksums.csv -NoTypeInformation

Get-ChildItem results -Recurse -File |
  Get-FileHash -Algorithm SHA256 |
  Export-Csv metadata\results_checksums.csv -NoTypeInformation
```

4. Review code cleanup items in [CODE_READABILITY_NOTES.md](CODE_READABILITY_NOTES.md).

5. Initialize local Git from the repository root.

```powershell
git init -b main
git status --short
```

6. Stage and commit the scaffold plus project contents you want versioned.

Recommended tracking policy:

- include `scripts/`
- include `data_raw/` files that the scripts actually consume
- include release-grade `results/`
- exclude `data_intermediate/`

```powershell
git add .
git commit -m "Initial reproducible pd-envtox pipeline"
```

7. Create the GitHub repository named `pd-envtox`.

GitHub web UI:

- create a new repository named `pd-envtox`
- do not let GitHub auto-create a README, license, or `.gitignore`

GitHub CLI alternative:

```powershell
gh repo create pd-envtox --source . --private --remote origin --push
```

If you want it public, change `--private` to `--public`.

8. If you created the repository in the web UI, connect and push manually.

```powershell
git remote add origin https://github.com/<YOUR-USERNAME>/pd-envtox.git
git push -u origin main
```

9. Tag the first frozen release.

```powershell
git tag -a v0.1.0 -m "First reproducible baseline"
git push origin v0.1.0
```

10. Validate reproducibility from a clean clone.

In a fresh clone:

```powershell
Rscript -e "renv::restore()"
Rscript run_pipeline.R
```

## Time Estimate

- scaffolding and repo metadata review: `20-40 minutes`
- `renv` initialization and package snapshot: `30-90 minutes`
- first Git commit and GitHub push: `15-30 minutes`
- readability cleanup before public release: `60-120 minutes`

Realistic total: `2-5 hours`

