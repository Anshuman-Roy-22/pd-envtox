# AAN non-writing revision and completion audit

Audit date: 2026-08-13

Audited branch: `aan-v2-overhaul`

Comparison base: `origin/main` at `04df4ed`

Audited head before this document: `bad812c`

## Executive verdict

The scientific analysis and validation expansion is complete. The entire
non-writing submission package is not yet complete.

The branch now contains the full planned sequence of panel validation,
pathway-level reanalysis, independent toxicant testing, held-out testing, and
donor-level human substantia nigra testing. It also preserves negative results
and uses frozen analysis plans for the major confirmatory steps.

Three technical deliverables remain before the non-writing work can be called
submission-ready:

1. Create final V2 figures that show the actual revised evidence. There are
   currently no figure files under `results/v2/`.
2. Connect scripts 16 through 27 to a master V2 runner and add acquisition for
   the older untracked GEO inputs. The current `run_pipeline.R` stops at script
   15 and therefore reproduces the original analysis only.
3. Complete release verification with checksums for the newest result
   directories and a clean-clone end-to-end rerun.

No further mechanism search or additional dataset testing is required by the
current scientific plan. Any new confirmatory dataset would be a new extension,
not completion of the existing overhaul.

## Completion matrix

| Non-writing component | Status | Evidence or remaining requirement |
| --- | --- | --- |
| Original 21-gene panel validation | Complete | Human PD, pesticide, regional-vulnerability, and matched-null tests are tracked. |
| Empirical-Bayes differential expression | Complete | Limma-confirmed gene-level analyses are tracked for the major microarray cohorts. |
| Continuous network-score evaluation | Complete | Four prespecified genome-wide predictive tests are tracked. |
| Pathway-level convergence analysis | Complete | Signed Reactome convergence and targeted PD-family results are tracked. |
| Independent rotenone replication | Complete | Frozen GSE116280 primary contrast and declared sensitivities are tracked. |
| Formal multitoxicant analysis | Complete | Seven GSE150005 contrasts were analyzed with voom/limma and CAMERA. |
| Held-out mechanism validation | Complete | GSE196190 primary and three external sensitivities were completed. |
| Human disease-tissue validation | Complete | GSE178265 primary and GSE243639 independent replication were completed at donor level. |
| Frozen-plan and metadata-lock discipline | Complete | Six frozen plans, a formal reanalysis specification, and a primary metadata lock are tracked. |
| Analysis result tables and audits | Complete | Gene, pathway, donor, mapping, QC, sensitivity, and manifest tables are tracked. |
| Final V2 scientific figures | Incomplete | `results/v2/` contains zero PNG, PDF, SVG, or TIFF figure artifacts. |
| One-command V2 reproduction | Incomplete | `run_pipeline.R` runs scripts 01 through 15 only. |
| Clean-clone raw-input acquisition | Partial | Fetchers exist for held-out and final human-SN inputs. Older V2 GEO inputs still depend on local untracked directories. |
| Checksum coverage | Partial | Checksum manifests are missing from `gse196190_heldout`, `heldout_external_sensitivity`, and `human_sn_proteasome`. |
| Final clean-clone verification | Incomplete | The complete V2 sequence has not been rerun from a fresh clone through one master command. |
| Branch isolation and traceability | Complete | Work is isolated on `aan-v2-overhaul`; `origin/main` was not modified. |

## Scope of the overhaul

Relative to `origin/main`, the audited branch contains:

- 18 analytical and protocol commits before this audit document.
- 152 changed files.
- 14 new analysis or input-acquisition scripts.
- 12 new protocol, metadata-lock, search-log, or validation-report files.
- 118 files under `results/v2/`.
- Approximately 115 MB of newly tracked inputs, outputs, code, and metadata.
- 590,095 inserted lines, most of which are full gene-level and pathway-level
  result tables rather than handwritten source code.

## Detailed revision record

### 1. Recovered human PD validation of the frozen network panel

Commit: `b48cfb3` (`Add recovered human PD validation checkpoint`)

What changed:

- Added `scripts/16_human_pd_validation.py`.
- Reconstructed outcome-blind probe mapping and preprocessing for GSE20141 and
  GSE7621 substantia nigra cohorts.
- Evaluated the frozen 21-gene panel in each cohort and with a directional
  Stouffer meta-analysis.
- Added a 10,000-draw expression-matched null comparison.
- Tracked cohort-level gene effects, meta-analysis effects, panel summaries,
  matched-null results, checksums, and a result README.

Result:

- 20 of 21 panel genes were evaluable in both cohorts.
- 14 of 20 had the same PD-control direction.
- Zero genes passed cohort-level or meta-analysis FDR.
- Direction agreement was not enriched against matched genes, p = 0.1974.
- Mean absolute meta-z was not enriched, p = 0.4579.
- Mean panel meta-z was negative but exploratory, two-sided matched-null
  p = 0.0522.

Impact on the project:

The frozen individual-gene panel did not receive convincing human PD support.
This prevented the project from presenting the original candidate list as a
validated disease signature.

### 2. GSE46798 factorial pesticide and genotype analysis

Commits:

- `c97a317` (`Freeze GSE46798 factorial validation plan`)
- `e7d7f34` (`Add GSE46798 factorial and module validation`)

What changed:

- Froze the design and evaluation rules before the final analysis.
- Added `scripts/17_gse46798_factorial_validation.py`.
- Separated pesticide exposure, SNCA A53T genotype, and
  genotype-by-exposure effects.
- Tested individual panel genes and fixed functional modules.
- Used expression-matched 10,000-draw null distributions for response
  magnitude and signed coordination.
- Added probe-mapping, all-gene, panel-contrast, module-test, sensitivity,
  checksum, and interpretation files.

Result:

- No fixed-panel gene was FDR-significant for pesticide exposure or the
  genotype-by-exposure interaction.
- No fixed module or full panel showed significant response magnitude or
  coordinated direction relative to matched genes.
- SNX2, WLS, and VPS26A were significant in the prespecified A53T-versus-
  corrected vehicle contrast. This was a genotype observation and did not
  validate pesticide responsiveness.

Impact on the project:

The analysis removed a confounding interpretation in which a genotype effect
could have been presented as a pesticide effect. It also supplied a formal
negative external test of the original panel.

### 3. GSE17542 substantia nigra versus VTA vulnerability analysis

Commits:

- `1641d4b` (`Freeze GSE17542 SN VTA module plan`)
- `bc1ad7c` (`Add GSE17542 SN VTA module analysis`)

What changed:

- Froze region, time, module, and matched-null decisions.
- Added `scripts/18_gse17542_sn_vta_modules.py`.
- Modeled MPTP responses in substantia nigra and VTA dopamine neurons,
  region-by-treatment interactions, and baseline regional differences.
- Added individual-gene, fixed-module, probe-manifest, and 10,000-draw
  matched-null outputs.

Result:

- 19 of 21 frozen candidates mapped unambiguously.
- No candidate passed genome-wide FDR in a declared contrast.
- Neither the full panel nor a fixed module was enriched for response magnitude
  or signed coordination.
- The smallest nominal module result was a baseline retromer/endosomal regional
  difference, p = 0.0618. It was not an exposure response.

Impact on the project:

Regional vulnerability did not rescue the original panel hypothesis.

### 4. Genome-wide continuous network predictive test

Commits:

- `bb15357` (`Freeze continuous network predictive test`)
- `ed5dd29` (`Evaluate continuous network predictive validity`)

What changed:

- Froze a secondary analysis asking whether continuous network proximity had
  predictive value even though the 21-gene cutoff failed.
- Added `scripts/19_continuous_network_predictive_test.py`.
- Tested more than 10,000 eligible genes per outcome.
- Adjusted rank relationships for STRING degree and GTEx substantia nigra
  expression.
- Evaluated four outcomes spanning pesticide exposure, MPTP response,
  SN-versus-VTA interaction, and human PD.

Result:

- Unadjusted network closeness was unrelated to absolute response in all four
  outcomes.
- One adjusted association was statistically nonzero in GSE46798:
  partial Spearman r = 0.03075, p = 0.00106, four-test FDR = 0.00425.
- The effect explained less than 0.1% of rank variation and did not replicate
  in the other three outcomes.

Impact on the project:

The analysis established that the minimum-distance network score had
negligible and non-generalizing transcriptomic predictive value. The branch
therefore does not use a weak genome-wide association to rescue the fixed
panel.

### 5. Reactome cross-model convergence and limma confirmation

Commits:

- `10de0c9` (`Freeze Reactome convergence plan`)
- `6a9058d` (`Add signed Reactome cross-model convergence`)
- `974a368` (`Confirm differential expression and pathway convergence with limma`)

What changed:

- Froze the pathway database, ranking rules, model classes, targeted PD-family
  set, and recurrence criteria.
- Added `scripts/20_reactome_cross_model_convergence.py`.
- Added `scripts/21_limma_confirmatory_analysis.R`.
- Replaced ordinary gene-level uncertainty estimates with limma
  empirical-Bayes statistics for GSE46798, GSE17542, GSE20141, and GSE7621.
- Recomputed Reactome enrichment from moderated rankings.
- Combined corrected-neuron pesticide exposure, mouse SN MPTP exposure, and
  human postmortem PD using signed cross-model meta-analysis.
- Tracked full enrichment tables, strict recurrence, targeted PD-family
  results, all-gene limma tables, checksums, and R session information.

Result:

- Proteasome assembly was the strongest interpretable PD-relevant convergence:
  - GSE46798 corrected DA neurons: NES = -1.496, nominal p = 0.0237.
  - GSE17542 SN dopamine neurons: NES = -1.422, nominal p = 0.0383.
  - Human PD substantia nigra meta-analysis: NES = -3.155, nominal p < 0.001.
  - Cross-model meta-FDR = 0.00203.
- Eight leading-edge genes recurred in all three rankings: POMP, PSMA4, PSMB4,
  PSMC5, PSMD1, PSMD2, PSMD4, and PSMG1.
- Neither toxin dataset produced a within-dataset FDR-significant proteasome
  result.

Impact on the project:

The central project claim shifted from individual network-prioritized genes to
context-dependent pathway convergence. The result was retrospective and
required independent confirmation.

### 6. Independent GSE116280 rotenone replication

Commits:

- `6a02457` (`Freeze GSE116280 rotenone replication plan`)
- `c090bb7` (`Add independent rotenone replication analysis`)

What changed:

- Froze a primary mature-LUHMES contrast before analysis.
- Added `scripts/22_gse116280_rotenone_limma.R`.
- Used limma empirical-Bayes differential expression and Reactome enrichment.
- Evaluated all declared differentiation-day, concentration, and exposure-time
  contrasts as labeled sensitivities.
- Tracked all-gene statistics, probe mapping, pathway results, checksums, and R
  session information.

Result:

- The primary day-15, 50 nM, 24-hour contrast was negative in direction but did
  not meet the nominal threshold: NES approximately -1.31, p approximately
  0.09.
- Day-8, 24-hour contrasts at 50 nM and 100 nM were nominally negative:
  p = 0.0172 and p = 0.0210 in the committed branch result.
- Direction and significance were not stable across all dose-time cells.

Impact on the project:

The independent replication threshold was not met. Earlier-culture nominal
signals were retained as evidence of context dependence.

### 7. Formal GSE150005 multitoxicant reanalysis

Commit: `e010709` (`Add formal GSE150005 multitoxicant analysis`)

What changed:

- Added a formal reanalysis specification recording which outcomes had already
  been inspected.
- Added `scripts/23_gse150005_multitoxicant_camera.R`.
- Tracked the raw count matrix and Reactome GMT used by the analysis.
- Analyzed 24 differentiated human LUHMES RNA-seq samples across rotenone,
  MPP+, paraquat, 6-OHDA, ziram, methylmercury, vehicle, and an MS-275
  comparator.
- Used TMM normalization, voom/limma, robust empirical Bayes, CAMERA pathway
  testing, within-contrast FDR, and a filtering sensitivity analysis.
- Added full gene tables, full pathway tables, sample and analysis manifests,
  pathway summaries, checksums, and R session information.

Result:

- Toxicants produced heterogeneous transcriptomic programs.
- No exact Reactome pathway passed FDR in two different environmental-toxicant
  contrasts within this study.
- Rotenone produced downregulated cilium, microtubule, trafficking, and
  aggrephagy signals.
- Proteasome assembly was null in rotenone, one-sided-down p = 0.240, and MPP+,
  p = 0.116.

Impact on the project:

The analysis identified a plausible cilium/Hedgehog/microtubule-trafficking
mechanism for a later held-out test. It also supplied another clear failure of
universal proteasome suppression.

### 8. Held-out cilium/Hedgehog/microtubule-trafficking validation

Commits:

- `53d23c1` (`Freeze held-out mechanism validation plan`)
- `ba47aa2` (`Lock primary held-out GSE196190 contrast`)
- `ca18391` (`Run held-out mechanism validation`)

What changed:

- Defined a frozen union of five Reactome pathways and a Down direction.
- Logged the held-out dataset search and exclusion decisions.
- Locked the GSE196190 primary contrast from metadata before target-expression
  testing.
- Added `scripts/24a_fetch_heldout_inputs.R` with fixed URLs and SHA256 checks.
- Added `scripts/24_gse196190_heldout_mechanism_validation.R`.
- Added `scripts/25_heldout_external_sensitivity.R`.
- Added `scripts/26_gse4773_borderline_sensitivity.R`.
- Tracked primary, constituent, leave-one-pathway-out, replicate-block,
  filtering, mapping, sample-QC, PCA, and external-sensitivity outputs.
- Updated the R environment lockfile for the required packages.

Result:

- GSE196190 primary: 202 analyzed genes, Down direction, one-sided p = 0.1221,
  five of five constituent mean statistics negative, frozen outcome
  `NOT_CONFIRMED`.
- Inferred replicate-block sensitivity: p = 0.0692, still nonsignificant.
- GSE229460 MPP+ sensitivity: p = 0.4005.
- GSE287941 6-OHDA sensitivity: p = 0.3924.
- GSE4773 did not support the declared Down mechanism at one, two, or four
  weeks.

Impact on the project:

The strongest newly nominated mechanism did not validate in held-out data. The
negative result was preserved instead of being relabeled after inspection.

### 9. Frozen donor-level human substantia nigra proteasome validation

Commits:

- `ad8d0b7` (`Freeze human SN proteasome validation`)
- `bad812c` (`Validate human SN proteasome mechanism`)

What changed:

- Froze Reactome Proteasome assembly `R-HSA-9907900` and the declared Down-in-PD
  direction.
- Defined donor, cohort, cell-type, covariate, pathway-score, missing-gene,
  sensitivity, and interpretation rules before target-expression analysis.
- Added `scripts/27a_fetch_human_sn_proteasome_inputs.py` with checksum
  verification.
- Added `scripts/27_human_sn_proteasome_validation.py`.
- Added frozen GSE178265 dopamine-neuron metadata and provenance documentation.
- Tested author-annotated GSE178265 DA nuclei as the primary cohort and
  author-annotated GSE243639 dopamine neurons as the independent replication.
- Used donor-level pathway scores and adjusted models rather than treating
  individual nuclei as independent samples.
- Added donor cell-count audits, donor scores, measured-gene audits, gene-level
  results, leave-one-gene-out results, subtype sensitivities, shared-leading-
  edge descriptions, input checksums, and a frozen outcome summary.

Result:

- GSE178265 primary: 6 PD donors and 8 controls, 48 measured pathway genes,
  adjusted beta = -0.4000, one-sided-down p = 0.0352, frozen outcome
  `CONFIRMED`.
- GSE243639 replication: 7 PD donors and 8 controls, 43 measured pathway genes,
  adjusted beta = +0.0042, one-sided-down p = 0.5095, outcome
  `NOT_CONFIRMED`.
- The frozen overall label is `PRIMARY_ONLY`.
- The primary effect was stable to leave-one-gene-out analysis, but no
  individual gene passed gene-level FDR.

Impact on the project:

The pathway gained one outcome-frozen human disease-tissue confirmation. It did
not gain independent cross-cohort human replication.

## Reproducibility and research-integrity revisions

The overhaul added the following safeguards beyond new datasets:

- Frozen analysis plans before major confirmatory analyses.
- A formal record distinguishing retrospective analyses from held-out tests.
- A held-out dataset search log containing inclusion and exclusion decisions.
- A metadata lock before the GSE196190 target-expression analysis.
- Explicit primary, gated-secondary, exploratory, and sensitivity labels.
- Gene-level FDR and pathway-level multiplicity control where applicable.
- CAMERA tests that account for intergene correlation.
- Empirical-Bayes limma confirmation for the major microarray analyses.
- Matched-gene null models for fixed-panel and module tests.
- Donor-level human models that avoid pseudoreplication across nuclei.
- Probe and identifier mapping manifests.
- Sample manifests, sample-QC tables, PCA outputs, and cell-count audits.
- Input SHA256 verification in the newer acquisition scripts.
- R session information and an updated `renv.lock` for major R analyses.
- Tracked null and sensitivity results, including failed confirmations.

## Current scientific conclusion supported by the branch

The original 21-gene network panel does not validate as a general toxicant-
responsive or human-PD panel. Continuous network proximity also lacks useful,
generalizing predictive strength.

Reactome Proteasome assembly is the strongest pathway-level observation. It
shows retrospective negative convergence across corrected human DA neurons,
mouse SN DA neurons, and human postmortem PD, with cross-model meta-FDR =
0.00203. Its frozen mature-LUHMES rotenone replication was nonsignificant. A
later frozen human SN donor analysis confirmed the declared direction in one
cohort and failed in an independent cohort.

The supported conclusion is context-dependent proteasome-pathway convergence
with one frozen human cohort confirmation. The evidence does not establish a
universal environmental-PD mechanism or independent cross-cohort human
replication.

## Remaining non-writing work

### A. Final V2 figure set

Status: required.

At minimum, the technical figure set should contain:

1. A study-design and evidence-flow schematic separating retrospective,
   frozen-primary, and sensitivity analyses.
2. A cross-model proteasome convergence panel with model-specific NES and
   nominal p values.
3. The frozen replication panel showing the GSE116280 primary and declared
   dose-time sensitivities.
4. A donor-level human substantia nigra panel showing GSE178265 and GSE243639
   effect estimates with confidence intervals and individual donor scores.
5. A compact validation-outcome matrix showing which candidate and mechanism
   claims were confirmed, partially supported, or not confirmed.

Every plotted number should be generated directly from tracked V2 tables.

### B. Master V2 runner

Status: required.

The existing `run_pipeline.R` runs the original scripts through script 15. A
new V2 orchestrator should:

- Restore or verify the R environment.
- Verify Python dependencies.
- Acquire every required untracked input.
- Run scripts 16 through 27 in dependency order.
- Keep retrospective, confirmatory, and sensitivity stages visibly separated.
- Stop immediately on missing inputs, checksum mismatches, or failed commands.
- Regenerate the final V2 figures and checksum manifest.

### C. Complete input acquisition

Status: required for clean-clone reproducibility.

Fixed acquisition exists for the held-out inputs and final human SN inputs.
The following locally present raw directories are still untracked and lack one
consolidated acquisition path:

- `data_raw/human_validation/`
- `data_raw/gse46798/`
- `data_raw/gse17542/`
- `data_raw/gse116280/`

The V2 runner should download and checksum-lock the required GEO series
matrices and platform annotations for these analyses.

### D. Release checksums and environment capture

Status: required.

Add output checksum manifests for:

- `results/v2/gse196190_heldout/`
- `results/v2/heldout_external_sensitivity/`
- `results/v2/human_sn_proteasome/`

Also record a Python environment manifest for the Python analyses and verify
that all checksum manifests match the release files.

### E. Fresh-clone verification

Status: required last step.

The technical release should be tested in a new clone using only tracked files
and documented downloads. Completion requires:

- Every required input downloads or is already tracked.
- Every checksum passes.
- Every V2 script exits successfully.
- Key frozen result labels and numeric assertions reproduce.
- Final figures regenerate from the tracked tables.
- `git status` remains clean after the verification workflow, except for files
  explicitly documented as generated and ignored.

## Worktree preservation and tracking

The audit excludes the user's preexisting modified and untracked research
files. They were not staged into the AAN commits. At the audited head, those
paths included local raw-data directories, a modified GSE187012 archive, a
modified GSE116280 result table, an untracked ordinary-statistics Reactome
directory, and Python cache files.

All AAN overhaul commits listed above are on `aan-v2-overhaul`. No commit in
this revision sequence was pushed to `origin/main`.

## Definition of non-writing completion

The non-writing work can be marked complete when all five conditions are true:

1. The existing analyses and frozen outcomes remain unchanged.
2. The final V2 figures are generated from tracked tables.
3. A master V2 command can acquire inputs and reproduce the analyses.
4. Release checksums and environment records cover all final outputs.
5. A clean-clone verification passes and the verified branch state is committed.
