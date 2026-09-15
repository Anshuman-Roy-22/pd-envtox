# Mechanistic redesign after the neuronal guide-key audit

Date: 2026-09-15. Status: **data acquisition and feasibility completed;
biological reanalysis and model evaluation not completed**.

## Decision

The uploaded Tian 2019 guide key does not supply the missing proteasome
perturbations. It contains none of the fixed 52 genes or shared eight leads.
Two detected raw guide sequences remain unidentified. We cannot validate a
model on experimental target responses that have not been established.
The earlier antioxidant-withdrawal survival comparison remains nonsignificant
(two-sided P = 0.14558), and human validation remains PRIMARY_ONLY.
The upload does not materially strengthen the present mechanistic claim.

The next justified question is:

**Which interventions produce coherent functional rescue in human
SNCA-triplication dopaminergic neurons, and which merely change an
alpha-synuclein-associated readout?**

The proposed original contribution is an independently specified analysis
of rescue coherence and its limits across measured endpoints. Its novelty
and strength still need to be demonstrated. It is not a claim to have
discovered a new neuroprotective drug, the published uncoupling mechanism,
or a universal proteasome defect in PD.

This is a new exploratory extension authorized by the user's request to
reapproach. It does not replace the fixed targets, endpoints, thresholds or
outcomes of earlier analyses. It shifts from testing a particular transcript
set to evaluating the consequences of actual interventions in a PD model.

## Why this is closer to a mechanistic test

The selected study uses patient-derived SNCA-triplication and genetically
corrected isogenic human midbrain dopaminergic neurons. Researchers actually
applied compounds and measured cellular responses. Public measurements include
alpha-synuclein staining, dopaminergic and neurite features, live-cell counts,
oxidative stress and mitochondrial membrane potential. Other experiments in
the study measured respiration, ATP, neuronal activity and protein levels.

These interventions and orthogonal assays supply a stronger experimental
foundation than a disease-versus-control RNA association alone. However,
compound responses do not identify a unique molecular target. Shared endpoints
do not establish mediation, and different assays on separate cultures are not
paired measurements of the same cells. A direct genetic or pharmacological
blockade/rescue experiment would still be needed to establish necessity of a
proposed causal link.

The authors already report Tyrphostin A9 and analogue effects, mitochondrial
uncoupling, and alpha-synuclein reduction. Those are their findings. Repeating
them is reproduction, not a new discovery or an unseen confirmation.

## Public-data audit

Primary source: Gorgogietas et al., *Scientific Reports* (2025),
DOI https://doi.org/10.1038/s41598-025-14735-0.
Author repository: https://github.com/Ksilink/Notebooks/tree/main/Neuro/PD_MorphProfileScreening.

Sixteen files totaling 8,659,066 bytes were obtained directly from the author
repository and locked to commit
`1ed064f3d7881e7f520ae6afefcc8aed75ffc6ba` (2025-06-25).
Their individual immutable URLs and SHA256 values are in `input_manifest.json`.
`metadata_inventory.json` records actual file headers and design labels.

| Deposited component | Verified structure | Permitted use / limitation |
|---|---|---|
| Initial isogenic/control imaging | 144 well rows, four named plates, mutant vehicle, mutant Prostratin and corrected vehicle | Establish phenotype axes; keep wells grouped by plate. Article describes four independent experiments, but identifiers must be reconciled before claiming a donor/experiment count from filenames |
| Primary compound screen | 4,800 well rows on 20 plates; 4,080 compound-map rows | Broad intervention screen. Article describes 1,020 compounds; map concentrations include both 2.5 and 5 micromolar, so do not assume every compound was tested at 5 |
| Dose-response follow-up | 930 rows, four plates, compound identities and actual concentrations | Within-study follow-up of author-selected compounds, not an unbiased independent cohort; contains processed features and prior classifier distances |
| ROS compound workbooks | 20 named compound sheets plus an empty sheet; corresponding tidy export has 19 named compounds | Real oxidative-stress measurements. Replicate identifiers are absent from the tidy export; TIC10 occurs in the workbook but not that export |
| JC-1 compound workbooks | 20 compound sheets; corresponding tidy export has 19 named compounds | Membrane-potential assay. Same identifier problem and TIC10 discrepancy. Lower JC-1 is not automatically beneficial |
| Analogue viability file | 30 rows, six compound codes at three concentrations, genotype controls | Must resolve compound-code identities, normalization and experimental blocks before inference |
| Repository proteasome assays | Two activity files with genotype, MG132 and three Tyrphostin dose columns | Not a drug-by-MG132 factorial design: there are no combined Tyrphostin-plus-MG132 columns. Cannot establish proteasome dependence of drug rescue |

File row counts are not counts of independent biological replicates. The
biochemical workbooks have repeated control columns across compound sheets.
The author conversion notebooks concatenate sheets, melt values and strip
asterisk annotations. The exact experimental blocks and annotations need
reconciliation; copied controls must not be counted repeatedly as new samples.
No biological replication is inferred from every three adjacent rows.

The repository snapshot predates the final publication. It contains proteasome
assay files under Fig6, while the final paper's Figure 6 describes additional
protein/autophagy measurements. Treat these as versioned author-deposited data,
not an automatically exact reconstruction of the final article's figure panels.

## Concrete next implementation

1. **Reconstruct experimental units and compound identities.** Join by stable
   catalogue/compound identifiers and plate/well keys. Reconcile renamed salts,
   compound spellings, actual concentrations, shared controls, and the missing
   TIC10 export. Preserve source exclusions and annotation flags. Resolve raw
   biochemical biological-replicate membership from deposited metadata if
   possible; otherwise retain those assays as descriptive evidence.

2. **Specify separate phenotype axes before scoring compounds.** Candidate
   deposited fields are alpha-synuclein intensity/localization;
   `Cell_Neurites_LengthPerNuclei_MAP2`;
   `Cell_Neurites_LengthPerNuclei_TH`;
   `Cell_Intensity_MeanIntensity_TH`; and
   `Nuclei_Number_Living` / `Nuclei_Ratio_Living`.
   A fall in total fluorescence alone is insufficient evidence of rescue.
   Evaluate per-cell measures, live-cell abundance and neurite integrity
   separately. Do not silently discard toxic treatments before assessing
   whether a marker-only score misclassifies them.

3. **Use an interpretable control-trained baseline.** Learn the disease-to-
   corrected direction from untreated genotype controls, with scaling and
   feature selection inside training folds. Compare an alpha-synuclein-only
   score with a multivariate phenotype score. Keep plate/experiment groups
   separate during validation; random well splitting risks leakage. Do not
   use the authors' precomputed classifier distances as raw independent
   outcomes for evaluating a replacement classifier.

4. **Test follow-up performance before mechanistic interpretation.** Carry
   the fixed scores into the separate dose-response data. Match doses
   explicitly and aggregate technical wells within the proper experiment.
   Report performance on all eligible matched compounds, including failures.
   Author-selected follow-up compounds introduce selection bias and cannot
   estimate performance across the whole chemical library without adjustment
   or an appropriately qualified estimand.

5. **Compare independent assay modalities.** Ask whether phenotype rescue
   corresponds to ROS and mitochondrial-function changes while preserving
   viability. Do not equate any decrease in membrane potential with rescue:
   uncoupling can reduce ROS and ATP together. Assay directions must follow
   their control calibration. Dose, timing, normalization and missingness
   must match the question. An association across separate assays remains
   cross-assay concordance, not established mediation or a causal chain.

6. **Freeze the evaluable statistical question after the unit audit.** Lock
   eligible samples, compounds, doses, folds, primary endpoint, comparator,
   effect thresholds, uncertainty method and multiplicity policy in a
   separate commit before deriving new target effects. The source literature
   is already known, so label this a prespecified reanalysis of published
   experiments, not a prospectively unseen discovery. If biological units
   cannot be recovered, do not manufacture inferential P values.

This specification deliberately does not invent a final biological success
threshold or claim a frozen confirmatory analysis before the experimental
units have been resolved. The next executable scientific task is the unit
and compound reconciliation in step 1, not another gene-set significance scan.

## Tools and compute

- Python, pandas and openpyxl: source reconciliation, assay tables and provenance.
- scikit-learn: control-trained regularized linear baselines with grouped
  validation, fitted only after the analysis specification is locked.
- SciPy or statsmodels: effect estimation, dose-response comparisons and
  experimental-unit uncertainty where the design supports it.
- Matplotlib: dose-response plots, separate rescue axes and validation figures.

The downloaded data fit comfortably on the available CPU environment. No paid
GPU or State/GEARS training is required for this stage. More complex models
would need to outperform the simple baseline on the fixed validation task.
Molecular docking or predicted protein structures alone would not resolve
the neuronal rescue question.

## Alternatives considered and why they are not the selected next stage

| Source | Mechanistic value | Current decision |
|---|---|---|
| NERINE / CiS-CN SNCA-toxicity CRISPRi study, 2026 | Perturbation-by-SNCA-context design is closer to a specificity test | The article says complete screen details/hits are forthcoming. No complete usable screen-count resource was established in this audit; do not budget the month around unavailable counts |
| Sequential CRISPR / NatB-UBE2W study, 2024 | Published proteasome-inhibitor and double-knockdown experiments demonstrate the kind of dependency evidence required | Useful mechanistic reference and possible model benchmark. Its reported epistasis is already known and cannot be presented as our discovery or validation of the fixed 52-gene transcript set |
| Bidirectional alpha-synuclein aggregation screen, 2026 | Direct aggregation endpoint and neuronal follow-ups | Broad primary screens were in HEK-derived cells with technical duplicates; the article reports no primary genes at FDR <0.05. Do not assume this is a ready independent DA-neuron benchmark |

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC13347950/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10857481/
- https://doi.org/10.1002/2211-5463.70233

## Exposure and attribution record

For Tian 2019, only guide identities and guide-enrichment counts were used;
no survival-phenotype or RNA-expression values were analyzed. For the new
functional-rescue source, publication text, repository paths, headers,
design labels and author notebook code cells were inspected. Notebook output
cells were not used. A preliminary type-based workbook preview also exposed
some numeric control entries stored as strings with asterisks; this prevents
describing the entire audit as strictly phenotype-blind. No drug-effect
estimates, fitted model, ranked candidates or statistical tests were computed.

Published positive rescue findings were necessarily visible in the source
articles. They are attributed to their authors. The new analytical question
and its evaluation must be distinguished from those known results, and the
student must understand and personally defend any submitted interpretation.

## Reproduction of the completed stage

From the repository root, with the existing pinned openpyxl dependency:

```bash
python scripts/acquire_functional_rescue_metadata.py --fetch
```

This downloads missing sources into the ignored `data_raw/functional_rescue/`
cache, verifies every checksum, and regenerates `metadata_inventory.json`.
It never executes downloaded notebooks or computes new biological effects.
The input manifest and metadata inventory are tracked separately from raw data.

## What would materially strengthen the project

A useful result would show that a fixed computational measure predicts
coherent functional responses in data excluded from model fitting, retains
performance across experimental groups, and distinguishes marker changes
from impaired viability or loss of neuronal features. Stronger molecular
claims additionally require intervention-specific dependency or rescue
evidence. Positive training fit, a more elaborate model, or reproduction
of the paper's best-known drug finding alone would not meet that standard.

This approach has better mechanistic relevance to the user's goal, but it
has not yet demonstrated a new mechanism or an AAN-winning result. Winning
probability cannot be responsibly quantified from the present evidence.
The remaining competition month should prioritize this bounded evaluation,
interpretation and reproducibility over repeated post-null mechanism searches.
