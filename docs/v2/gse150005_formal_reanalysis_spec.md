# GSE150005 formal multitoxicant reanalysis specification

## Inferential status

This specification records the final reproducible analysis used for GSE150005.
Provisional gene-set outcomes had already been inspected before this document was
written. Therefore, GSE150005 is an exploratory/retrospective analysis and is not
described as a held-out or preregistered confirmation. Its role is to characterize
effect heterogeneity and nominate mechanisms for a later, genuinely unseen test.

## Dataset and contrasts

The analysis uses the GEO raw-count matrix
`GSE150005_Tong_et_al_DopaNeuron_PD_Tox_RAW_COUNT_MATRIX.txt.gz`. It contains 24
differentiated human LUHMES samples, with three replicates for DMSO and each of
seven perturbations.

The six PD- or environmental-toxicant contrasts are:

- rotenone versus DMSO
- MPP+ versus DMSO
- paraquat versus DMSO
- 6-hydroxydopamine versus DMSO
- ziram versus DMSO
- methylmercury versus DMSO

MS-275 versus DMSO is retained as a mechanistically distinct HDAC-inhibitor
comparator and is not counted as independent PD-toxicant replication.

## Preprocessing and differential expression

1. Remove rows without an HGNC symbol.
2. Sum raw counts for duplicated HGNC symbols.
3. Retain genes with at least 10 counts in at least three samples.
4. Apply TMM library normalization.
5. Fit a no-intercept cell-means model with `voom`, `lmFit`, contrasts against
   DMSO, and robust empirical Bayes moderation.
6. Control gene-level false discovery rate separately within each contrast using
   Benjamini-Hochberg adjustment.

## Competitive pathway testing

- Use the same official Reactome GMT snapshot used in the earlier cross-model
  analysis.
- Retain pathways with 10 to 500 measured genes.
- Run limma CAMERA on the voom expression matrix for every contrast.
- Estimate pathway-specific intergene correlation by setting
  `inter.gene.cor = NA`.
- Report two-sided competitive p-values, direction, estimated correlation, and
  Reactome-wide Benjamini-Hochberg FDR separately for every contrast.
- Interpret highly overlapping Reactome terms as pathway families rather than
  independent biological discoveries.

## Targeted summaries and sensitivity analysis

The output includes transparent retrospective summaries for proteasome assembly
and the mechanism families nominated after the full screen. For the previously
proposed downregulated proteasome direction, a one-sided p-value is derived from
the CAMERA two-sided result. Benjamini-Hochberg correction across the six
PD/environmental toxicants is reported; MS-275 remains separate.

As a preprocessing sensitivity analysis, proteasome assembly is retested after
edgeR `filterByExpr` filtering. No sensitivity result changes the exploratory
status of GSE150005.

## Interpretation boundary

A significant pathway in one GSE150005 contrast is evidence of a transcriptional
response in that condition. Nominal recurrence across several contrasts in this
single experiment is within-study support. Independent replication requires a
new dataset whose accession and outcomes were not inspected before freezing a
separate held-out validation plan.
