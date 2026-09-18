# Mechanistic reconstruction and matched-perturbation follow-up

Date: 18 September 2026. Surviving parent: `f86baa8`. Reconstruction plan: `232f19e`. Exact comparator lock, before comparator responses: `2a220e1`.

## Main finding

Naa20, Pomp and Psmc5 perturbations produced larger decreases in a frozen synaptic RNA score than their respective panels of other gene perturbations matched on an external neuronal-dependency score, baseline expression and mouse coverage. All three evaluable comparisons passed two-sided BH correction across four planned tests and retained their direction and significance in litter-level sensitivity analyses. Psmb4 had insufficient matched comparators and remains **NOT_EVALUABLE** for this test.

This is a material improvement in perturbation specificity evidence. Its scope is comparison with these fixed panels. It does not demonstrate uniquely selective synaptic dysfunction, eliminate every form of cellular stress, or establish a toxin-to-NatB causal chain.

## Provenance and what was reconstructed

The later human limma audit, mouse bridge, overview and human specificity files reported in the conversation did not survive the workspace reset and were not in the user's local branch. This addition reconstructs their analyses from preserved inputs and public sources. It does not recreate their lost Git identities or claim their outcomes were unseen. The surviving branch was verified remotely at `f86baa8`; remote `main` remained `04df4ed`.

Human and focal mouse results were already known. Comparator selection rules were committed before computing the new comparator synaptic responses, followed by an exact gene-list checkpoint. This internal sequence does not turn a retrospectively selected data source into prospective external validation. A metadata-only clarification excluded Safe_target_N controls because the browser collapses their labels; molecular targets and non-targeting controls then agreed exactly across releases.

## New mouse comparison

| Target | Matched genes | Paired mice | Additional synaptic decrease vs matched panel | 95% CI | Two-sided BH4 |
|---|---:|---:|---:|---|---:|
| Naa20 | 15 | 38 | -0.012189 | [-0.018540, -0.005837] | 0.000540 |
| Pomp | 15 | 19 | -0.035166 | [-0.045989, -0.024343] | 0.00000434 |
| Psmb4 | 1 | Unevaluable | NA | NA | NA |
| Psmc5 | 20 | 25 | -0.039533 | [-0.051363, -0.027704] | 0.00000157 |

Units are differences in mean natural-log(1 + 10,000 * gene UMI / total RNA UMI) scores. They are not percentage changes in synaptic function. Naa20's effect remains small. NA comparisons enter the four-test correction as p=1 without an estimated effect. There is no claim for all eight toxin-derived leads: Pomp and Psmc5 are two of those eight; Naa20 is the NatB perturbation. Psmb4 is a third lead with a reconstructed non-targeting contrast, but lacks adequate matched controls.

Each target is first compared with non-targeting cells within mouse and author neuronal group, weighting group contrasts by target-cell counts. Primary eligibility requires >=20 target neurons, >=80% retention and >=20 reference neurons per stratum. The paired comparison subtracts the equal mean of eligible comparator-gene effects in the same mouse, requiring >=5 comparators per mouse, >=8 mice and >=10 selected genes. Genes and cells are never treated as independent animals.

Controls were selected from independent Tian2021 plus-antioxidant CRISPRi neuronal-survival percentiles, unique author human/mouse orthology and metadata eligibility. All proteasome-pathway genes, NatB subunits and the original eight leads were excluded. The plan fixes dependency, expression and coverage distance bounds and nearest-neighbor selection. The exact list, every covariate and rejected shortlist members are retained in this directory. Expression matching used non-targeting cells only.

**Important matching limitations:** the dependency score comes from another cellular context. Naa20 and Psmc5 have percentiles about 0.127 and 0.108, whereas Pomp and Psmb4 are about 0.865 and 0.745. Calling all four essential genes in that reference is inaccurate. Matching did not equalize realized knockdown efficacy, guide toxicity, selective cell loss or every neuronal subtype's abundance. Candidate screening required retained-cell coverage, so the test concerns surviving observed neurons. Comparator panels are finite and partly overlapping. Matching limits were fixed; they were not widened after Psmb4 failed.

Litter-level paired results retain eight litters per evaluable target: Naa20 difference -0.015672, BH4 0.00558; Pomp -0.032046, BH4 0.00921; Psmc5 -0.044973, BH4 0.00400. Leave-one-comparator-out means are saved. The prespecified broad dependency-bottom-quintile challenge also retains negative differences, but is a descriptive sensitivity drawn from the shortlist, not an expression-matched or genome-wide essential-gene comparison. It does not rescue Psmb4's primary eligibility failure.

## Reconstructed human proteomic specificity

The original source contains five non-targeting and five NAA25-knockdown iNeuron cultures. Robust limma tests 4,511 proteins; 994 pass FDR<0.05. SYP decreases by -1.306 log2 units (FDR 0.00367) and SYT1 by -1.214 (FDR 0.00864). Strict-backbone and sample-omission models are retained separately, including eligibility failures. SNCA's primary proteome-wide FDR is 0.0627.

The complete-case synaptic score decreases more than all non-synaptic proteins and more than detectability-matched non-synaptic proteins:

| Panel | Background | Relative log2 difference | 95% CI | Exact permutation BH4 |
|---|---|---:|---|---:|
| Chemical synapse, 87 proteins | All non-synaptic | -0.318928 | [-0.417302, -0.220554] | 0.015873 |
| Chemical synapse, 87 proteins | Matched | -0.287148 | [-0.384352, -0.189944] | 0.015873 |
| Neurotransmitter release, 18 proteins | All non-synaptic | -0.384609 | [-0.643752, -0.125466] | 0.015873 |
| Neurotransmitter release, 18 proteins | Matched | -0.335816 | [-0.592889, -0.078742] | 0.015873 |

Backgrounds exclude both synaptic annotations, NAA25 and SNCA. Five nearest proteins per member are chosen using label-blind pooled abundance and peptide-backbone count ranks. Reused matched proteins are not extra biological replicates. Inference uses the ten culture scores and all 252 group assignments. The displayed CIs are Welch intervals; permutation p values answer the corresponding label-exchangeability test. The overlapping panels and backgrounds are dependent, and the four-test family is explicit.

Within-culture subtraction cancels a shared additive log-scale offset. Persistence after this subtraction argues against a uniform offset being the sole explanation. It does not exclude biological cell-health effects or proteome composition changes that differ by functional category.

**The competitive enrichment limitation remains:** chemical-synapse CAMERA FDR is 0.000675 at fixed correlation 0.01 and 0.626543 when correlation is estimated. No claim of robust proteome-wide competitive synaptic enrichment is supported across both settings. The positive relative-score tests do not replace that null. The full seven-panel score tests, all eligible Reactome CAMERA pathways and respiratory/translation results remain available.

## Reconstruction checks

- Mouse release reconciliation: 3,390,123 eligible cells across all retained molecular targets/non-targeting controls; 3,387,094 browser matches, 3,029 missing, zero analyzed mouse/target/group disagreements. Focal-plus-control subset: 219,610 eligible and 219,422 matched cells.
- Unique orthology retains 85/87 human synaptic genes and 37/39 proteasome genes. Exact membership and missing symbols are saved.
- An independently read raw-count shard verifies the browser scale: 4,207 shared cells, 3,679 nonzero Xkr4 measurements, maximum absolute error about 6.05e-7.
- Reconstructed focal non-targeting synaptic effects: Naa20 -0.013029 (41 mice), Pomp -0.032593 (20), Psmb4 -0.077027 (13), Psmc5 -0.039219 (25). All pass the original BH8 two-module family; core previously reported values agree within 1e-10.
- Proteasome RNA scores increase in all four perturbations. This direction differs from the earlier toxin-associated assembly decrease. Compensation is a hypothesis, not a demonstrated mechanism.
- Ninety-nine four-guide non-targeting pseudo-groups yield three nominal positives; binomial excess p=0.87755. This checks calibration only. These tests share controls and are dependent.
- Minimum-ten-cell, Syt1-omission, single-mouse omission, individual diagnostic genes, guide coverage and litter analyses are retained. Naa20 perturbation increases Snca RNA, so it does not replicate the human SNCA protein decrease at the RNA level.

## Reproduce this addition

Use Python 3.12 with `requirements.txt` in this directory. The human R audit requires **R 4.3.3, limma 3.58.1 and statmod 1.5.0**; the script checks these exactly. R must be available as `R`, or pass `--r /absolute/path/to/R`. This historical audit environment is distinct from the older project's separate R 4.5.2 reproduction gate. Use Linux/WSL for a consistent runtime; allow several GB of disk and memory.

```bash
python -m pip install -r docs/v2/reconstruction_20260918/requirements.txt
python scripts/run_reconstruction_20260918.py --cache /absolute/path/to/reconstruction-cache --r R
```

The command downloads and verifies the public metadata, one raw calibration shard and only required expression vectors; reconstructs the human results; rechecks and recreates comparator selection; runs mouse inference; and verifies scientific artifact checksums. Existing sources that disagree with their hashes stop the run. UCSC source URLs are mutable, so the recorded hash is required in addition to the URL. The Hugging Face revision is immutable. Use `--offline` with a complete verified cache. Use `--verify-only` to check packaged code, frozen selection and outputs without analysis or downloads. `R_session.txt` is retained as a runtime record and excluded from cross-platform byte matching.

The input manifest includes the original NatB acquisition manifest and frozen Reactome/neuronal-survival tables. The public supplementary zip supplies the author litter metadata. This command reproduces this addition from those preserved upstream results; it does not rerun the entire earlier toxin and survival-analysis project.

## Figure and files

`results/v2/reconstruction_20260918/specificity_comparison.png` and `.pdf` show the four human relative contrasts and the mouse non-targeting/matched-perturbation comparisons. Bars are unadjusted 95% confidence intervals. Human tests use five cultures per group; mouse non-targeting sample sizes are 41/20/13/25 and matched sample sizes 38/19/unevaluable/25 in Naa20/Pomp/Psmb4/Psmc5 order. The two panels use different measurements and scales and must not be interpreted as directly comparable magnitudes. Human CAMERA's null and Psmb4's missing specificity test are displayed. The generated TSVs supply exact numbers, membership, animal-level values and all sensitivities.

## Sources and permissible interpretation

The NatB rationale and experiments originate in [Santhosh Kumar et al., Science Advances, DOI 10.1126/sciadv.adj4767](https://www.science.org/doi/10.1126/sciadv.adj4767). The mouse experiments and broad proteasome/synaptic observations originate in [Shi et al., DOI 10.64898/2026.03.16.711480](https://doi.org/10.64898/2026.03.16.711480), an unreviewed preprint used here. The atlas concerns young mouse brain and is not a Parkinson's disease or toxicant model. External dependency scores originate in the already-documented Tian2021 analysis; see `results/v2/neuronal_survival/README.md`.

The new contribution is the explicit cross-source hypothesis, auditable reconstruction, frozen gene panels, mouse-level estimation, relative proteomic contrasts and prespecified matched-perturbation comparison. The source experiments and their established findings must be credited to their authors. These results support convergent synaptic transcriptional/proteomic responses to NatB/proteasome perturbation with partial specificity evidence. They do not establish proteasome activity, direct regulation of all eight leads, α-synuclein-mediated causality, synaptic physiology, a rescue or a treatment effect.

## Verification completed

An isolated staged-index snapshot reran the full command using the verified public-input cache. All 25 locked code/input artifacts and all 51 scientific output artifacts matched byte-for-byte, including PNG and PDF. `verification.json` records the snapshot, runtime and scope. This verifies this addition with cached inputs; it does not close the older whole-project R clean-clone gate or establish fresh-environment portability.
