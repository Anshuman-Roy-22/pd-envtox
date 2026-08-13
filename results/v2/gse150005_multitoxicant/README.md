# GSE150005 formal multitoxicant reanalysis

## Inferential status

This is a reproducible retrospective analysis. Provisional outcomes had been
inspected before the formal specification was recorded, so none of these results
is labeled held-out confirmation. The analysis characterizes multitoxicant
heterogeneity and nominates hypotheses for a later unseen dataset.

## Analysis

- 24 differentiated human LUHMES RNA-seq samples, three per condition.
- Raw counts collapsed to HGNC symbols.
- 16,322 genes retained after requiring at least 10 counts in at least three
  samples.
- TMM normalization, voom/limma differential expression, and robust empirical
  Bayes moderation.
- Competitive Reactome testing by CAMERA with pathway-specific intergene
  correlation and separate within-contrast FDR control.
- MS-275 retained as an HDAC-inhibitor comparator rather than counted as a
  PD/environmental toxicant.

## Full-screen results

| Contrast | FDR-significant genes | FDR-significant Reactome pathways |
|---|---:|---:|
| Rotenone | 965 | 18 |
| MPP+ | 2,151 | 0 |
| Paraquat | 0 | 0 |
| 6-hydroxydopamine | 1 | 2 |
| Ziram | 155 | 8 |
| Methylmercury | 47 | 0 |
| MS-275 comparator | 8,181 | 211 |

No exact Reactome pathway reached FDR <= 0.05 in two different
PD/environmental-toxicant contrasts. Therefore, nominal same-direction recurrence
inside GSE150005 is described as within-study support rather than independent
replication.

## Exploratory mechanism families

Rotenone produced a coherent downregulated cilium, microtubule, and trafficking
family. Significant terms included connexon trafficking (FDR 0.00395),
intraflagellar transport (FDR 0.0112), periciliary cargo trafficking (FDR 0.0151),
tubulin post-translational modification (FDR 0.0186), and Hedgehog "off-state"
machinery (FDR 0.0461). Intraflagellar transport, periciliary trafficking, and
Hedgehog "off state" also showed nominal same-direction results in MPP+ and
ziram, but did not pass their contrast-wide FDR thresholds.

Aggrephagy was downregulated in rotenone (FDR 0.0272), nominally down in ziram
(p = 0.00459), and borderline in MPP+ (p = 0.0535). This is a heterogeneous
candidate rather than a universal mechanism.

ATF6 chaperone-gene activation was significant in 6-hydroxydopamine (FDR 0.0161)
and nominally up in rotenone and methylmercury. The HRI/heme-deficiency response
was significant in ziram (FDR 0.00516) and nominally up in MPP+, paraquat, and
6-hydroxydopamine. These results nominate an exposure-dependent integrated/ER
stress response but do not provide external replication.

## Proteasome result

Proteasome assembly did not reproduce in GSE150005 rotenone:

- direction: down
- CAMERA correlation: 0.0943
- two-sided p = 0.480
- prespecified one-sided-down p = 0.240
- Reactome-wide FDR = 0.919

MPP+ was also null (one-sided-down p = 0.116). The edgeR `filterByExpr`
sensitivity analysis produced similar results: rotenone p = 0.238 and MPP+
p = 0.117. A provisional matched-gene test that ignored intergene correlation is
not used for inference.

## Reproduction

From the repository root after restoring the recorded R environment:

```powershell
Rscript scripts/23_gse150005_multitoxicant_camera.R
```

The formal specification is
`docs/v2/gse150005_formal_reanalysis_spec.md`. Exact pathway, gene, sample,
filtering, and sensitivity outputs are stored in this directory.
