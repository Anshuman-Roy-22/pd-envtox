# Held-out mechanism-validation report

Analysis date: 2026-08-13 EDT.

This report implements the frozen protocol in
`docs/v2/heldout_mechanism_validation_frozen_plan.md`. The candidate search and
the exact GSE196190 contrast were committed at `ba47aa2` before any expression
values were opened. All evaluable results are retained regardless of direction
or significance.

## Primary result

The locked comparison used day-54 human iPSC-derived dopaminergic neurons after
the sustained 100 micromolar MPP+ exposure versus its matched control endpoint,
with three biological samples per group.

| Criterion | Result | Pass |
| --- | ---: | :---: |
| Measured genes in the frozen union | 202 | Yes |
| CAMERA direction | Down | Yes |
| One-sided-down p | 0.1221 | No |
| Constituent pathways with negative mean moderated t | 5 of 5 | Yes |

Frozen outcome: **NOT_CONFIRMED**.

The directional pattern is internally coherent, but it does not meet the
precommitted inferential threshold. The simple count-filter sensitivity was
similar (Down, one-sided p = 0.1170). Every leave-one-pathway-out union remained
Down, with one-sided p values from 0.0851 to 0.1366.

The sample filenames also encode three consistent replicate series
(S6/S7, S24/S25, and S42/S43). Because GEO does not expose a formal pairing
field, the locked unpaired model remains primary. A declared alternative-valid-
handling sensitivity that blocked these filename-derived series was still above
the confirmatory threshold (Down, one-sided p = 0.0692; 5 of 5 negative
constituent means).

Two constituents were nominally Down before correction: intraflagellar
transport (one-sided p = 0.0188) and microtubule-dependent connexon trafficking
(one-sided p = 0.0352). Their two-sided Benjamini-Hochberg values across the five
constituents were both 0.1758, so neither is a multiplicity-adjusted replicated
pathway.

## Gated secondary result

The frozen Aggrephagy test was not run. Its hierarchical gate remained closed
because the primary mechanism was not confirmed. No secondary p-value was
generated.

## Declared external sensitivities

The same union test was applied to all stronger eligible external sensitivities
found by the frozen search. Neither passed the primary-style criteria:

| Dataset | Model and contrast | Direction | One-sided-down p | Negative constituent means |
| --- | --- | :---: | ---: | ---: |
| GSE229460 | Differentiated LUHMES, 10 micromolar MPP+ for 48 hours, 3 vs 3 | Down | 0.4005 | 1 of 5 |
| GSE287941 | Wild-type differentiated LUHMES, 6-OHDA, 4 vs 4 | Down | 0.3924 | 1 of 5 |

The lower-confidence GSE4773 SK-N-MC sensitivity also did not support the frozen
direction. At one week the union was significantly **Up** rather than Down
(two-sided p = 0.0241; one-sided-down p = 0.9880). The two- and four-week
one-sided-down p values were 0.3540 and 0.9334, respectively, with zero of five
negative constituent means at every time point.

These sensitivity results cannot change the locked primary label. They also do
not provide a consistent external directional pattern.

## Interpretation

The cilium/Hedgehog/microtubule-trafficking family remains a strong retrospective
rotenone-associated pattern in GSE150005, but this prospective held-out exercise
does not establish it as a generalized toxicant mechanism. The primary dataset
shows directional consistency without statistical confirmation, while the
additional LUHMES and SK-N-MC sensitivities do not reproduce that consistency.

No new independently replicated biological mechanism was established by this
test. Existing significant results elsewhere in the repository remain
retrospective, within-study, or cross-model convergence findings and should not
be relabeled as held-out confirmation.

## Reproduction

After restoring `renv.lock`, run from the repository root:

```powershell
Rscript scripts/24a_fetch_heldout_inputs.R
Rscript scripts/24_gse196190_heldout_mechanism_validation.R
Rscript scripts/25_heldout_external_sensitivity.R
Rscript scripts/26_gse4773_borderline_sensitivity.R
```

Exact sample manifests, input checksums, annotation-database checksums, package
versions, complete gene statistics, pathway statistics, and sensitivity outputs
are stored under `results/v2/gse196190_heldout/` and
`results/v2/heldout_external_sensitivity/`.
