# Requested specificity follow-up: completed human checks and blocked mouse metric

Plan committed before new calculations at `ee992a0`, parent `1c7983e`, on 18 September 2026. Known prior outcomes remain known; these analyses are sensitivity checks, not an independent validation cohort.

## Status

1. **Global disruption and realized-knockdown regression: BLOCKED_INPUT.** No such regression has been fitted and no result is implied. The available committed mouse tables contain selected readouts, not broad expression centroids. The previous 245-vector cache was both absent after workspace reset and selected around target/endpoint genes. Reusing that selected panel as a global-disruption measure would be circular. Public-source file metadata identifies `analysis/preprocessed_gex.zarr.tar.gz` (35,059,391,830 bytes) and batch h5ad files at the frozen Hugging Face revision. They were not downloaded during this follow-up. A broad expression acquisition or user-supplied compatible matrix is required. The full metric, covariates and interpretation are fixed in PLAN.md. The existing matched-panel result is unchanged.
2. **Permutation-consistent human intervals: complete.** All four pointwise intervals remain below zero under the constant additive-shift/exchangeability model.
3. **CAMERA fixed-correlation curve: complete.** Sensitivity is substantial; the estimated-correlation result remains nonsignificant. The installed implementation reveals a degrees-of-freedom difference in addition to differing correlation estimates.
4. **Dependency-versus-synaptic-effect plot: complete from preserved aggregate tables.** This is the requested descriptive subpart of item 1. It does not stand in for the blocked disruption check.

The code was restored directly from the user's already-pushed `1c7983e` GitHub branch. No PowerShell command was needed for that recovery. The missing untracked cloud cache cannot be recovered by `git fetch` or a Git bundle. No further broad search or large expression download was undertaken after diagnosing that gap.

## Human intervals

The original exact p values and four-test BH correction are unchanged. The new confidence sets invert the same absolute difference-of-means permutation statistic over all 252 assignments after subtracting each candidate constant shift from the treated cultures. Boundaries are calculated from the exact piecewise-linear breakpoints; no Monte Carlo sampling or endpoint grid approximation is used. Ties have tolerance 1e-12. The observed confidence sets are connected. The previous Welch intervals remain in the table for explicit comparison, and the new figure uses the inverted intervals.

| Contrast | Effect, log2 units | Permutation-inverted 95% interval | Exact two-sided p | BH4 |
|---|---:|---|---:|---:|
| Chemical synapse minus all non-synaptic | -0.318928 | [-0.432478, -0.217127] | 0.007937 | 0.015873 |
| Chemical synapse minus matched background | -0.287148 | [-0.384266, -0.185484] | 0.007937 | 0.015873 |
| Neurotransmitter release minus all non-synaptic | -0.384609 | [-0.655458, -0.148280] | 0.015873 | 0.015873 |
| Neurotransmitter release minus matched background | -0.335816 | [-0.591878, -0.101422] | 0.015873 | 0.015873 |

These are pointwise intervals, not simultaneous intervals adjusted across four contrasts. Their model assumes that removing a constant additive effect makes cultures exchangeable. This does not remove cell-health or normalization/composition limitations of the underlying experiment. With equal groups, each assignment has a complementary assignment giving the same absolute statistic. The minimum attainable **two-sided** p is therefore **2/252=0.0079365**, not 1/252. The latter is attainable for a directional extreme. More Monte Carlo draws cannot improve the finite assignment resolution. The common BH value follows the observed discrete p values and four-test correction; it is not a general minimum BH value for every possible four-test dataset.

## CAMERA sensitivity

All 862 eligible pathways are retained at every grid point, with a separate full-family BH adjustment at each point. The primary complete matrix still contains 3,409 proteins and ten cultures. No cutoff, pathway membership or contrast was changed.

| Pathway | Last passing tested fixed correlation | FDR there | First failing tested fixed correlation | FDR there | Estimated correlation | Estimated-mode FDR |
|---|---:|---:|---:|---:|---:|---:|
| Chemical synapse | 0.03 | 0.02707 | 0.05 | 0.11340 | 0.013414 | 0.62654 |
| Neurotransmitter release | 0.02 | 0.04500 | 0.03 | 0.07254 | 0.178004 | 0.62654 |

The curve reports the tested grid. No precise interpolated transition is asserted. At the originally used fixed correlation 0.01, FDRs reproduce 0.00067525 and 0.0237531.

**The estimated mode is not merely another point on this common-correlation curve.** In the installed limma 3.58.1 `camera.default`, fixed-correlation, non-rank inference uses `df.camera = G - 2`, or 3,407 here. Estimated-correlation inference uses `min(df.residual, G - 2)`, or 8 here. Each pathway also receives its own estimated correlation and the resulting full-family p distribution changes BH adjustment. The small correlation change from 0.01 to 0.0134 alone does not explain the chemical-synapse result's large change. The exact installed function is preserved in `installed_camera_default.R.txt` for inspection. This observation clarifies the methods; it does not establish that either favorable setting is the correct one or erase the estimated-mode null.

## Comparator dependency plot

The union contains 51 distinct comparator genes. Each point averages its eligible mouse-level synaptic effects equally; repeated cells or mice are not counted as additional gene observations. Comparator-only Pearson r is **0.2636** (descriptive p=0.0617), and Spearman rho is **0.0735** (descriptive p=0.6083). The focal genes lie below the comparator cloud in this plot. Psmb4 is shown descriptively and retains its original matched-test eligibility failure.

Do not write "there is no dependency trend." The Pearson estimate is compatible with a moderate positive relationship and the panels have restricted, bimodal dependency coverage by construction. Genes also share mice and controls, so these association p values are descriptive. Neither the scatter nor the nonsignificant associations excludes global disruption or a moderate dependence trend. Measured global disruption and realized target suppression remain untested under the new model.

## Reproduction and verification

Use the preserved primary matrices and current committed comparator tables. Python 3.12 with numpy, pandas, scipy and matplotlib as recorded in `verification.json`; R 4.3.3 with limma 3.58.1 and statmod 1.5.0. From the repository root:

```bash
R --vanilla --slave -f scripts/specificity_followup_20260918.R --args .
python scripts/specificity_followup_20260918.py
```

The R command writes the full grid and estimated-mode results. Python computes exact inverted intervals, saves every evaluated breakpoint, creates the dependency table and renders both figures as 300-DPI PNG and vector PDF. Run the commands in this order. No new mouse expression download is needed for these completed portions. The main mouse disruption model has no executable fitting script yet because its input matrix is unavailable; its documented plan must not be mistaken for a completed test.

Intervals were independently checked for common-offset invariance, treatment-shift equivariance and label-reversal symmetry on all four contrasts. Direct permutations on shifted cultures verify the null, effect estimate and interval endpoints. The fixed-0.01 and estimated CAMERA outputs reproduce the earlier results. Both figures passed visual inspection. Input/output SHA256 values and rerun status are recorded with this follow-up. This is a targeted sensitivity verification, not the repository-wide clean-clone gate.
