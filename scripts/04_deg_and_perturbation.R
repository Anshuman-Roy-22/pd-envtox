# scripts/04_deg_and_perturbation.R
suppressPackageStartupMessages({
  library(Matrix)
  library(limma)
  library(data.table)
})

dir.create("data_intermediate", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
source("scripts/helpers_10x_io.R")

# Helpers
minmax01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# PART A: GSE17542 (bulk microarray) - limma DEG + perturbation
bulk <- readRDS("data_intermediate/gse17542_clean.rds")
expr <- bulk$expr_gene             # genes x samples
meta <- bulk$meta
stopifnot("group" %in% colnames(meta))

# design for 3 groups (control, mptp_2d, mptp_10d)
group <- droplevels(meta$group)
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("^group", "", colnames(design))

fit <- lmFit(expr, design)

# contrasts vs control
contr <- makeContrasts(
  mptp_2d_vs_ctl  = mptp_2d - control,
  mptp_10d_vs_ctl = mptp_10d - control,
  levels = design
)

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

tab2d  <- topTable(fit2, coef = "mptp_2d_vs_ctl",  number = Inf, sort.by = "none")
tab10d <- topTable(fit2, coef = "mptp_10d_vs_ctl", number = Inf, sort.by = "none")

tab2d$gene  <- rownames(tab2d)
tab10d$gene <- rownames(tab10d)

# perturbation = |log2FC| * -log10(FDR)
tab2d$perturb  <- abs(tab2d$logFC)  * (-log10(pmax(tab2d$adj.P.Val, 1e-300)))
tab10d$perturb <- abs(tab10d$logFC) * (-log10(pmax(tab10d$adj.P.Val, 1e-300)))

# Normalize within each comparison (0-1)
tab2d$perturb_norm  <- minmax01(tab2d$perturb)
tab10d$perturb_norm <- minmax01(tab10d$perturb)

# Save
fwrite(as.data.table(tab2d),  "results/gse17542_limma_mptp2d_vs_control.csv")
fwrite(as.data.table(tab10d), "results/gse17542_limma_mptp10d_vs_control.csv")

saveRDS(list(
  mptp2d_vs_ctl  = tab2d,
  mptp10d_vs_ctl = tab10d
), "data_intermediate/gse17542_deg_tables.rds")

cat("Saved GSE17542 DEG tables.\n")

# PART B: GSE187012 (10x CTL vs EXP) - pseudo-bulk log2FC + perturbation
in_dir <- "data_raw/GSE187012/extracted"

ctl_feat <- file.path(in_dir, "GSM5667021_CTL_features.tsv.gz")
ctl_bc   <- file.path(in_dir, "GSM5667021_CTL_barcodes.tsv.gz")
ctl_mtx  <- file.path(in_dir, "GSM5667021_CTL_matrix.mtx.gz")

exp_feat <- file.path(in_dir, "GSM5667022_EXP_features.tsv.gz")
exp_bc   <- file.path(in_dir, "GSM5667022_EXP_barcodes.tsv.gz")
exp_mtx  <- file.path(in_dir, "GSM5667022_EXP_matrix.mtx.gz")

stopifnot(file.exists(ctl_feat), file.exists(ctl_bc), file.exists(ctl_mtx))
stopifnot(file.exists(exp_feat), file.exists(exp_bc), file.exists(exp_mtx))

# read matrices
mat_ctl <- read_10x_counts(ctl_feat, ctl_bc, ctl_mtx)
mat_exp <- read_10x_counts(exp_feat, exp_bc, exp_mtx)

# pseudo-bulk: sum counts across all cells per condition
pb_ctl <- Matrix::rowSums(mat_ctl)
pb_exp <- Matrix::rowSums(mat_exp)

# align genes
genes_common <- intersect(names(pb_ctl), names(pb_exp))
pb_ctl <- pb_ctl[genes_common]
pb_exp <- pb_exp[genes_common]

# log2FC with pseudocount
pseudo <- 1
log2fc <- log2((pb_exp + pseudo) / (pb_ctl + pseudo))

gse187012 <- data.table(
  gene = genes_common,
  log2FC = as.numeric(log2fc),
  perturb = abs(as.numeric(log2fc))
)
gse187012[, perturb_norm := minmax01(perturb)]

fwrite(gse187012, "results/gse187012_pseudobulk_log2fc.csv")
saveRDS(gse187012, "data_intermediate/gse187012_pseudobulk.rds")

cat("Saved GSE187012 pseudobulk log2FC.\n")

cat("\nDONE: 04_deg_and_perturbation\n")


