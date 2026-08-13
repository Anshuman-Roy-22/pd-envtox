# Limma empirical-Bayes confirmation

All frozen contrasts were rerun with limma 3.66.0 under R 4.5.2.

- GSE46798: no fixed candidate is FDR-significant for either pesticide effect or
  the genotype-by-exposure interaction. SNX2, WLS, and VPS26A remain significant
  for the prespecified A53T-versus-corrected vehicle contrast.
- GSE17542: no fixed candidate is FDR-significant in any SN, VTA, interaction,
  or baseline-region contrast.
- Human postmortem meta-analysis: no fixed candidate passes meta-FDR. TOMM20 is
  the strongest lead (meta-p = 0.00240, meta-FDR = 0.133).

These results confirm the prior gene-level interpretation while replacing the
ordinary-model statistics for final use.
