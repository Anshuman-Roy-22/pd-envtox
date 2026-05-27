# scripts/03a_fetch_gse187012_meta.R
suppressPackageStartupMessages(library(GEOquery))

dir.create("data_raw", showWarnings = FALSE)
dir.create("data_intermediate", showWarnings = FALSE)

# This pulls the SOFT / series matrix content (often has sample/annotation tables)
gse <- getGEO("GSE187012", GSEMatrix = TRUE)

saveRDS(gse, "data_intermediate/gse187012_seriesmatrix.rds")
cat("Saved: data_intermediate/gse187012_seriesmatrix.rds\n")
