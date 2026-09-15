# Molecular selectivity of a proteasome-dependent intervention

Frozen 2026-09-15 before calculating intensity differences. This is an explicitly retrospective reanalysis: the authors' positive SNCA result and mechanistic experiments are known. Source headers, gene identities, peptide sequences and detection counts have been inspected. The intended contribution is a quantitative test of selectivity relative to the original proteasome-assembly set and neuronal identity, not discovery of the NatB mechanism.

## Experimental basis and input

Santhosh Kumar et al., Science Advances (2024), DOI 10.1126/sciadv.adj4767, supplementary Table S5. The source contains ten named MS runs, NTC_01 through NTC_05 and NAA25_06 through NAA25_10. Methods identify biological samples as separately processed neuronal cell pellets. They do not establish ten donors or independent differentiations. Treat the runs as five culture samples per condition, with independence across culture samples an explicit assumption. No paired design is inferred from the run numbers. This is not a toxin-exposure or survival experiment.

A gene-target coverage audit has established that Table S1 covers all eight original leads, but Tables S2 and S3 do not contain those eight. Do not treat the targeted neuronal screen as a negative result for untested genes. Source S5 has one internal SNCA peptide, with valid quantities in 5/5 control and 4/5 knockdown runs. NAA25 has seven precursor entries and at least one valid quantity in 5/5 control and 4/5 knockdown runs. POMP, UBE2W and PLK2 are not quantified. These omissions cannot be imputed as null effects.

## Fixed biological contrasts

The three directional primary comparisons are:

1. NAA25 knockdown lowers NAA25 abundance (manipulation check).
2. SNCA decreases relative to the fixed 52-gene Reactome Proteasome assembly set.
3. SNCA decreases relative to the fixed neuronal panel: MAP2, TUBB3, RBFOX3, SNAP25, SYP, SYN1, DLG4, SYT1.

The proteasome set is copied exactly from the earlier human-validation gene inventory, with no selection by this proteomics result. Require at least 20 quantified module genes. The neuronal panel includes the absent RBFOX3; require at least six quantified genes. SNCG is a predefined family-member specificity comparator, reported separately. PSMG1 and all eight original leads receive descriptive coverage/effect records, without promoting a favorable lead to a primary hypothesis.

## Preprocessing, fixed before intensity inspection

- Preserve the source CSV and checksums. No raw-MS reprocessing or claim of independently reidentifying spectra.
- Use only unambiguous single-gene rows: split semicolon-separated Genes labels, and retain rows whose nonempty unique label set has size one.
- Match Qvalue and TotalQuantity columns by exact run identifier, not incidental order. Valid observations require finite positive quantity and numeric precursor Qvalue <= 0.01. `Filtered`, zero and missing entries remain missing.
- Keep a precursor if it has valid observations in at least 8/10 runs, using a group-blind rule. Deduplicate exact gene/precursor identities only if every source value agrees; otherwise stop for a design audit.
- Log2-transform valid quantities. Subtract each precursor's median across its observed runs. For each gene and run, take the median centered precursor value, then center that gene by its median across observed runs before any sample-loading adjustment. Reuse these fixed gene intercepts in all normalization variants. Retain the source's distinct charge/modification features, and report unique unmodified peptide-backbone counts separately. Peptides are measurement features, never independent biological replicates.
- Estimate sample loading offsets from the median gene score among genes observed in all ten runs with at least two distinct peptide backbones, excluding SNCA and NAA25. Subtract each run's offset. This is relative-abundance normalization, not absolute protein per cell.
- A module uses the fixed subset of its genes observed in every run. Use the fixed, already-centered gene values and take the per-run median across module genes. This prevents changes in which genes are present from changing module composition between samples.
- A primary contrast requires at least four finite observations per condition. Analyze the observed SNCA runs; do not fill the missing NAA25_08 SNCA quantity. Do not drop NAA25_07 from other endpoints merely because its NAA25 quantity is missing.

## Inference and operational decision

For each contrast compute the knockdown-minus-control mean log2 difference and a Welch 95% confidence interval. Enumerate every assignment of the observed sample scores to groups of the same sizes and calculate a one-sided exact permutation p for a decrease, including the observed assignment. Use a tolerance of 1e-12 for ties. Adjust the three primary p values with Holm's method. This culture-level exchangeability calculation does not establish donor generalization, and nonrandom missingness remains a limitation.

To test preservation rather than infer it from a nonsignificant difference, calculate a Welch 90% confidence interval for each module's mean difference. A module is equivalent only if that entire interval lies within +/-log2(1.25). The symmetric log margin denotes ratios from 0.8 to 1.25 and is a prespecified practical abundance margin, not a clinical safety threshold. Module preservation is relative to the measured proteome and does not demonstrate proteasome assembly, catalytic activity, neuronal survival or lack of toxicity.

Call this `SELECTIVITY_SUPPORTED_SINGLE_PEPTIDE` only if all three Holm p values are <=0.05, both SNCA-relative contrasts are <=-log2(1.25), and both module equivalence criteria pass. If any calculable criterion fails, use `SELECTIVITY_NOT_ESTABLISHED`. Missing primary coverage yields `NOT_EVALUABLE`. The single-peptide qualifier is mandatory even for a positive result. Report each criterion separately; do not replace an unsuccessful endpoint.

## Fixed supporting analyses and sensitivities

- Report all measurable gene effects with Welch two-sided p values and BH correction across the measured proteome, with >=4 observed samples in each group. These are exploratory protein estimates.
- Repeat the primary contrasts without extra loading normalization. The two relative contrasts should be invariant to a common sample offset; verify this algebraically.
- Repeat with the per-run mean rather than median across module genes, retaining the original membership.
- Repeat excluding each run in turn, report effect directions and magnitudes only. Do not demand or select significant leave-one-out p values.
- Require two peptide backbones for SNCA as an evidence-quality sensitivity. This is known to be unevaluable from the source and will remain so.
- Report the original eight genes' existing human-cohort results alongside coverage in the intervention data. These are preexisting observational results, not new validation or pooled with the intervention p values.
- Extract all eight results from the genome-wide SKMEL30 screen as descriptive published source statistics. Do not combine melanoma and neuronal p values or interpret failed targeted-screen inclusion as neuronal evidence.

## Mechanistic interpretation and stop conditions

The authors' published NAA25/UBE2W double knockdown and proteasome-inhibitor experiments provide a dependency anchor. Our analysis can assess relative molecular selectivity, not reproduce that factorial dependency test from unavailable replicate-level data. The POMP/NRF2/PLK2 experiments in Bi et al. (2021) are external prior work and must be credited separately. Their existence does not validate disease causality of the project's observational proteasome signature.

Preserve the result and all failed criteria. Any further route must address a distinct measured biological question, retain the previous results, and be qualified for actual input access before execution. Do not alter this analysis to secure a positive label. Changes needed to correct coding or source interpretation must be documented.
