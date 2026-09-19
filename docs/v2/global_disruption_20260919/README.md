# Rebuilt Item 1 runner

**The observed global-disruption regression is NOT RUN.** This package restores the missing runner. It is a new implementation based on the committed original plan, not a recovery of the lost 7e2c86b Git object. Its history starts at the published 0d93272 commit. Keep aan-v2-overhaul separate from main.

## What was recovered and verified

- Original prespecified analysis, 51-comparator union, four focal perturbations, and original mouse-level synaptic results.
- Fourteen original H5AD source files, total 77,099,607,247 bytes, pinned to the immutable source revision and publisher SHA256 values in batch_sources.json.
- Complete source metadata, browser cell labels and author litter workbook. The restored workbook matches the previously recorded inner-file checksum. The historical supplementary ZIP hash discrepancy is preserved in litter_source.json.
- All 1,591 eligible target-mouse observations, with original target-cell counts and retention requirements unchanged.
- 19,070 biological genes, 85 synaptic genes and 18,931 global-background genes. The browser `_range` metadata entry is excluded.
- This implementation drops controls in strata with no eligible target observation, which have zero weight in every retained contrast. This reduces the selected-cell set from the prior runner's 264,688 to 245,249, without changing the inferential sample or weighting. Control centroids in every used stratum retain all eligible control cells.
- First complete raw batch downloaded again, publisher SHA256 verified, selected cells reduced, and Xkr4 independently compared with the original browser expression vector. Raw-to-browser discrepancy <=5.97e-7; weighted contribution discrepancy <=5.94e-8.
- Profiled REML agrees with statsmodels on synthetic data (maximum fixed-prediction discrepancy approximately 1.32e-9). Repeated-litter bootstrap identifiers and prediction path were checked on synthetic data.

These checks validate components. The remaining 13 real batches, full synaptic reconstruction, observed fit, 1,000 observed-data bootstrap refits and actual-result figure have not run here. The download package contains software and verification records, not a biological regression result or a reusable first-batch cache.

## Windows execution

Requires Python 3.12; recommend 16 GB RAM and at least 30 GB free disk on the cache drive. Source downloads total about 77 GB. Raw batches are processed individually; the largest is about 9.74 GB. Download assembly temporarily keeps chunk files and the assembled batch. The full dataset need not occupy disk at once.

After importing the replacement bundle and pushing aan-v2-overhaul, run these commands one at a time:

```powershell
cd "$env:USERPROFILE\pd-envtox"

py -3.12 -m venv "$env:USERPROFILE\pd-item1-env"
& "$env:USERPROFILE\pd-item1-env\Scripts\python.exe" -m pip install -r docs/v2/global_disruption_20260919/requirements.txt

& "$env:USERPROFILE\pd-item1-env\Scripts\python.exe" scripts/run_global_disruption.py --cache "$env:USERPROFILE\pd-item1-rebuilt-cache" --export "$env:USERPROFILE\Downloads\pd-envtox-item1-results.zip" --verify-only

& "$env:USERPROFILE\pd-item1-env\Scripts\python.exe" scripts/run_global_disruption.py --cache "$env:USERPROFILE\pd-item1-rebuilt-cache" --export "$env:USERPROFILE\Downloads\pd-envtox-item1-results.zip"
```

If a command fails, stop and send its output. Do not use the older Python 3.10 research environment. No R installation is needed. The runner checks pinned package versions and source/code checksums, including normalization of Windows Git line endings.

Use an initially empty dedicated cache folder outside the project. Keep PowerShell open and prevent the computer sleeping during execution. Interrupted work resumes by rerunning the same final command. Only checksum-verified reduced batches and verified download chunks are reused. A changed code/input fingerprint stops reuse. The bootstrap checkpoints every ten draws.

The runner deletes each raw H5AD from its own cache after verified reduction. It preserves reduced checkpoints until the results have been audited. Do not delete the cache during execution or before reviewing a failure. After final verification and backup, the dedicated cache and environment can be removed to recover disk space; those cleanup commands should be based on the actual folder paths.

Return pd-envtox-item1-results.zip for review. A ZIP may contain an incomplete run: consult execution_status.json and model_status.json. NOT_EVALUABLE is a valid scientific outcome and must remain reported. A completed run still needs scientific and visual audit before submission.

## Analysis and reporting

The original plan is docs/v2/specificity_followup_20260918/PLAN.md, committed at ee992a0. Recovery-specific implementation details were committed at d92db36 before reduction or model execution. No new hypothesis or success threshold was introduced.

Calculate within-mouse/group contrasts from all eligible cells, combine group contrasts using the original target-cell weights, then calculate RMS across the frozen background. Save realized target-expression contrast, matched control mean and detection. Use the fixed measurability thresholds of control mean >=0.1 and detection >=10%.

Model comparator observations only: synaptic effect ~ global RMS + realized target-expression difference + log target-cell count + mouse fixed effects, with perturbed-gene random intercept and REML. Continuous predictors are centered/scaled for numerical conditioning, preserving the same model. Require >=20 comparator genes with >=3 measurable mice each, full-rank fixed design and converged finite fit. Validate the observed REML fit against statsmodels. Failed gates produce NOT_EVALUABLE without a simpler replacement model.

Predict focal Naa20, Pomp and Psmc5 using population gene intercept zero. Preserve unmeasurable observations and absent mouse levels as unavailable; flag predictor extrapolation. Psmb4's metrics are retained descriptively and do not restore its failed matched-comparator eligibility. Residual summaries average within mice, then litters. Whole-litter resampling uses 1,000 draws, seed 20260918, distinct labels for repeated litter copies and unchanged gene identities; >=900 successful fits and available focal summaries are required for descriptive percentile intervals.

The zero-random-variance boundary is allowed within the same REML model. No single nominal p value or residual interval establishes specificity, excludes sickness/viability confounding, or demonstrates causal mediation. Compare all available and in-range residual summaries, and retain failure/exclusion records. Previous matched-panel results remain unchanged.

## Source provenance

Raw counts: perturbai/wholebrain_crispr_atlas, immutable Hugging Face revision a7c65dc0a64da4bd47cf6ef5f4dec6c7ef745e87. Original browser scale: UCSC whole-brain-perturb/combined, checksums inherited from the prior reconstruction. Source experiments are the public Shi et al. whole-brain mouse CRISPR atlas, DOI 10.64898/2026.03.16.711480, an unreviewed preprint. This secondary analysis does not constitute a new animal experiment or a direct PD/toxin model.

## Windows HTTPS certificate fix

The first Windows execution successfully acquired all_obs.parquet, then failed certificate-chain validation at the UCSC barcode URL. Downloads now use truststore 0.10.4 with the native system trust mechanism (Windows CryptoAPI). Hostname and certificate verification remain mandatory. No unverified-TLS fallback or HTTP substitution is provided. This transport-only change preserves source hashes, statistical code, comparator selection, thresholds, prepared inputs and reduced-checkpoint compatibility. Existing verified downloads are reused. Install the updated requirements after importing the certificate-fix commit, then resume with the same cache path.

The Linux transport and application integrity can be checked here; Windows-specific behavior must be confirmed by the local resumed run. A continued failure should be reported with the log, without altering certificate-verification settings. Implementation reference: [Truststore documentation](https://truststore.readthedocs.io/en/latest/).
