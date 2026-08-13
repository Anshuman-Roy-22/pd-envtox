# Human substantia nigra proteasome-validation report

Analysis date: 2026-08-13 EDT.

This analysis implements
`docs/v2/human_sn_proteasome_frozen_plan.md`. The protocol was committed as
`ad8d0b7` before any expression value for a target-pathway gene was requested.
The exact hypothesis was lower expression of Reactome Proteasome assembly
(`R-HSA-9907900`) in Parkinson disease dopaminergic nuclei.

## Confirmatory results

| Cohort | Donors | Measured genes | Adjusted PD beta | One-sided-down p | Frozen outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| GSE178265 primary DA nuclei | 6 PD, 8 control | 48 of 52 | -0.4000 | 0.0352 | `CONFIRMED` |
| GSE243639 independent DA replication | 7 PD, 8 control | 43 of 52 | +0.0042 | 0.5095 | `NOT_CONFIRMED` |

Overall frozen outcome: **PRIMARY_ONLY**.

The primary coefficient is a 0.400 log2-unit decrease in the unweighted
geometric-mean pathway score after adjustment for age, sex, and postmortem
interval. Its two-sided p-value is 0.0704 and its 95% two-sided confidence
interval is -0.8412 to +0.0411. The one-sided result is confirmatory because the
down direction and one-sided decision rule were committed before target
expression was accessed.

The independent replication estimate is centered almost exactly on zero. Its
95% interval is -0.3787 to +0.3871, so the replication excludes neither a
moderate decrease nor increase but provides no directional evidence for the
frozen mechanism.

## Primary-cohort audit

- All 6 PD and 8 control donors passed the five-DA-nucleus minimum.
- Forty of 48 measured pathway genes had negative adjusted PD coefficients.
- No individual gene passed two-sided Benjamini-Hochberg FDR 0.05. The result is
  a prespecified pathway-level signal, not a claim about individually significant
  genes.
- The unadjusted pathway estimate remained negative at -0.4310, with
  one-sided p = 0.0745. Covariate adjustment changed the estimated effect only
  slightly and reduced residual variance.
- Every leave-one-gene-out score remained negative and passed the frozen
  one-sided threshold. Betas ranged from -0.4211 to -0.3809 and p-values from
  0.0262 to 0.0440.
- A separate base-R `lm` calculation reproduced the primary coefficient,
  standard error, t statistic, degrees of freedom, and p-value exactly.

These checks show that the primary result is coordinated across genes and is
not caused by one pathway member. They do not supply independent replication.

## Replication and declared sensitivities

| Analysis | Adjusted PD beta | One-sided-down p | Interpretation |
| --- | ---: | ---: | --- |
| GSE243639 DA nuclei, five-cell minimum | +0.0042 | 0.5095 | Confirmatory replication miss |
| GSE243639 DA nuclei, one-cell minimum | -0.0172 | 0.4832 | Null with 13 PD and 11 controls |
| GSE243639 all neurons | -0.2419 | 0.2498 | Negative but imprecise |
| GSE243639 DA nuclei, RIN added | +0.0139 | 0.5251 | Null |
| GSE178265 PD plus LBD versus control | -0.1367 | 0.3308 | Null |

The GSE243639 unadjusted DA-neuron result was also null (beta = -0.0091,
one-sided p = 0.4776). Only 22 of 43 measured genes had negative adjusted
coefficients, which is consistent with the near-zero pathway estimate.

The previously shared eight-gene leading-edge score was negative in both
cohorts, but neither result was significant: GSE178265 beta = -0.2583,
one-sided p = 0.1935; GSE243639 beta = -0.1469, p = 0.1957.

Several GSE178265 subtype analyses had nominal one-sided p-values below 0.05.
Most contained only three PD donors, all were declared descriptive, and no
subtype multiplicity correction was prespecified. They do not alter the
`PRIMARY_ONLY` label.

## Interpretation

The added analysis establishes statistically significant, coordinated lower
Proteasome assembly RNA expression in one independent postmortem DA-neuron
cohort under a frozen donor-level test. It does not establish independent human
replication because a second, separately collected cohort was null under the
same direction, score, donor unit, and covariate structure.

This result is compatible with context-dependent proteostasis disruption. It is
not evidence that proteasome activity is universally reduced in PD, and it does
not prove impaired protein degradation or causality. Postmortem single-nucleus
data also measure surviving nuclei, not DA neurons lost before tissue
collection. The replication cohort contains few DA nuclei for many donors,
which limits precision, but the prespecified one-cell and all-neuron
sensitivities did not reveal a hidden confirmatory effect.

The project conclusion therefore remains: Proteasome assembly is the strongest
recurrent mechanism in the current evidence stack, with one frozen human
disease-tissue confirmation, but it is not independently replicated across the
two single-nucleus cohorts.

## Reproduction and tracked artifacts

From the repository root with Python packages `numpy`, `pandas`, `scipy`, and
`openpyxl` installed:

```powershell
python scripts/27a_fetch_human_sn_proteasome_inputs.py
python scripts/27_human_sn_proteasome_validation.py
```

The fetcher verifies all three GSE243639 input checksums, validates the current
SCP1768 cluster against the frozen tracked cell metadata, and validates all 52
target-gene responses. Raw downloads and the derived GSE243639 count cache are
ignored by Git. The frozen target-independent GSE178265 cell metadata is tracked
under `metadata/human_sn_proteasome/`.

Complete model results, donor scores, donor and cell-count audits, gene-level
results, measured-gene decisions, subtype sensitivities, leave-one-gene-out
results, and an input manifest are tracked under
`results/v2/human_sn_proteasome/`.
