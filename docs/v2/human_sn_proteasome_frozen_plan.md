# Frozen human substantia nigra proteasome-validation plan

Frozen on 2026-08-13 before requesting expression for any member of the target
pathway or opening the GSE243639 count matrix. A technical API check used `TH`,
which is not a member of the target pathway. No proteasome-expression value was
inspected before this lock.

## Objective and fixed hypothesis

Test whether the exact Reactome pathway `Proteasome assembly` is transcriptionally
downregulated in dopaminergic nuclei from postmortem Parkinson disease substantia
nigra pars compacta. This is the final disease-tissue validation of a mechanism
nominated before these two single-nucleus datasets were selected.

Pathway membership is the 52-gene set `R-HSA-9907900` in the tracked
`data_raw/pathways/ReactomePathways.gmt` file, SHA256
`89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b`.
Symbols are uppercased. No gene may be added or removed on the basis of either
cohort's result. The declared disease direction is down.

## Cohorts and biological units

### Primary: GSE178265 / SCP1768

- Human postmortem SNpc nuclei from Kamath et al. 2022.
- Use the authors' `UMAP: Human DA Neurons` membership. All ten author-defined
  DA subtypes are included.
- Compare donors labeled `PD` with donors labeled `Ctrl`. Donors labeled `LBD`
  are excluded from the primary contrast and retained for a declared sensitivity.
- A donor is eligible with at least five author-annotated DA nuclei. The metadata
  lock yields 6 PD and 8 control donors; all 14 pass this rule.
- The donor, never the nucleus or library, is the inferential unit.

The cell list, DA-subtype labels, donor status, age, sex, and postmortem interval
come from the public SCP1768 visualization API. Per-cell total UMI counts come
from the authors' official GitHub intermediate metadata at commit
`cae72e19290efe1ea8756e2032507905bdea52ba`.

### Independent replication: GSE243639

- Human postmortem SNpc nuclei from Martirosyan et al. 2024.
- Use the authors' `neurons00` cluster, which the paper identifies by `TH`,
  `SLC6A3`, `SLC18A2`, `ALDH1A1`, `KCNJ6`, and `AGTR1` as the dopaminergic
  neuronal population.
- Compare clinical diagnosis `Parkinson's` with `Control`.
- A donor is eligible with at least five `neurons00` nuclei. The metadata lock
  yields 7 PD and 8 control donors.
- The donor is the inferential unit.

Clinical covariates come from `GSE243639_Clinical_data.csv.gz`; cell membership
comes from the `Neurons` sheet in `GSE243639_UMAP_coordinates.xlsx`.

## Locked count processing and pathway score

For each cohort and cell selection:

1. Use the deposited raw integer counts for the 52 fixed pathway genes.
2. Sum each gene's counts across selected nuclei within each donor.
3. Sum total UMI counts across the same nuclei within each donor. GSE178265 uses
   the authors' per-cell `nUMI`; GSE243639 sums all rows of the deposited filtered
   count matrix.
4. A pathway gene is measured if it has at least 10 raw counts in at least three
   eligible donors. At least 10 measured pathway genes are required.
5. Calculate donor-by-gene expression as
   `log2(((gene_count + 0.5) / (total_UMI + 1)) * 1e6)`.
6. Define the donor pathway score as the unweighted arithmetic mean across all
   measured pathway genes. This is a log2 geometric-mean expression score. It
   keeps the donor sample size and naturally retains intergene correlation.

No cell, donor, or gene is removed using the disease effect, pathway score,
principal components, or expression-based outlier status.

## Locked models and decisions

The primary model in each cohort is ordinary least squares:

`pathway_score ~ PD + centered_age + sex + centered_PMI`

`PD` is 1 for Parkinson disease and 0 for control. Sex is a two-level factor.
All terms must be complete and the design must be full rank. The estimand is the
adjusted PD coefficient in log2 score units. The fixed one-sided-down p-value is
the lower-tail t probability for that coefficient.

The primary cohort is `CONFIRMED` only if all conditions hold:

- at least 10 pathway genes are measured
- adjusted PD coefficient is negative
- one-sided-down p < 0.05

The mechanism is `INDEPENDENTLY_REPLICATED` only if both cohorts independently
meet those conditions. A primary pass with a replication miss is `PRIMARY_ONLY`.
A down coefficient with one-sided p from 0.05 to less than 0.10 is
`PARTIAL_SUPPORT` for that cohort. Every other evaluable result is
`NOT_CONFIRMED`. Fewer than three eligible donors per group, fewer than 10
measured genes, a non-full-rank design, or unusable deposited data is
`NOT_EVALUABLE`.

The replication is run and reported regardless of the primary result. It cannot
rescue a failed primary result by pooling nuclei or by replacing the target
pathway. No additional biological mechanism will replace this final target after
outcome inspection.

## Declared sensitivity and descriptive analyses

These analyses cannot replace either confirmatory result:

- unadjusted `pathway_score ~ PD` in both cohorts
- GSE178265 PD plus LBD versus control
- GSE178265 DA-subtype summaries when at least three donors per group contribute
  five or more cells to a subtype
- GSE243639 all author-annotated neurons with the same five-cell rule
- GSE243639 `neurons00` with a one-cell minimum
- GSE243639 primary model with centered RIN added
- the eight previously shared leading-edge genes (`POMP`, `PSMA4`, `PSMB4`,
  `PSMC5`, `PSMD1`, `PSMD2`, `PSMD4`, `PSMG1`) as a descriptive score
- gene-level coefficients and leave-one-gene-out pathway scores

Two-sided p-values and Benjamini-Hochberg-adjusted gene-level results are
descriptive. No individual gene is required to pass FDR.

## Metadata artifacts inspected before the lock

- Reactome GMT, tracked in the repository
- SCP1768 public study, cluster, and annotation API responses retrieved
  2026-08-13
- Kamath et al. official GitHub metadata archive, SHA256
  `97b7cdd379507e279ba7168dc32ea6d4e78271a820bcf4348dd36bcf2a5fed55`
- GSE243639 clinical file, SHA256
  `ffea1163eb0c145d452fc08f8d4550986353c29b8b9496f1b69fa208296ebe24`
- GSE243639 UMAP workbook, SHA256
  `2fcd06645d8b5e4a32ec4462310b5571df5b657ce1e29c651c19b493a63035f4`
- GEO family SOFT records and supplementary-file listings for GSE178265 and
  GSE243639

Confirmation would establish a replicated coordinated RNA-expression decrease
in surviving postmortem DA nuclei. It would not establish reduced proteasome
activity, impaired protein degradation, causality, or the expression state of
DA neurons already lost before tissue collection.
