# scripts/01_fetch_geo.R

library(GEOquery)

dir.create("data_raw", showWarnings = FALSE)
dir.create("data_intermediate", showWarnings = FALSE)

### GSE17542 (Bulk Microarray, MPTP) ###
gse17542 <- getGEO("GSE17542", GSEMatrix = TRUE)
gse17542 <- gse17542[[1]]

saveRDS(gse17542, "data_intermediate/gse17542_eset.rds")

cat("GSE17542 downloaded and saved\n")


### GSE187012 (snRNA-seq, Paraquat/Maneb) ###
# GEO stores this as a supplementary TAR file
getGEOSuppFiles("GSE187012", baseDir = "data_raw")

cat("GSE187012 supplementary files downloaded\n")
