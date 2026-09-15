# A synaptic-protein selectivity concern in a proteasome-dependent intervention

**Result:** public human-neuron NAA25-knockdown proteomics supports a decrease in synaptic proteins despite preservation of the measured proteasome-assembly panel. The strongest robust protein result is SYT1. The frozen test for selective alpha-synuclein lowering failed: `SELECTIVITY_NOT_ESTABLISHED`.

This is a new analysis of someone else's intervention experiment. It adds a measured molecular consequence and a falsifiable mechanistic constraint to this project. It does **not** establish synaptic dysfunction, neurotoxicity, PD rescue or a newly discovered NatB clearance mechanism.

## What was analyzed

[Santhosh Kumar et al., Science Advances (2024)](https://doi.org/10.1126/sciadv.adj4767), supplementary Table S5, contains 48,028 precursor rows from five control and five NAA25-knockdown human iNeuron culture samples. These are not ten donors, and independence across differentiations is not documented. The published study supplies the experimental context; all estimates below are our reanalysis.

The authors' proteasome-inhibitor and NAA25/UBE2W double-knockdown experiments motivate a proteasome-dependent clearance mechanism. Those experiments were not numerically reanalyzed here, and their discovery belongs to the authors. The present question is whether this intervention preserves other neuronal proteins.

The analysis plan was committed at `a19b720` before calculating S5 intensity differences, after source metadata and the published findings had been inspected. Accordingly, this is a frozen retrospective reanalysis, **not** an unseen-dataset confirmation. The later peptide checks and mechanistic follow-up are explicitly post hoc.

## New quantitative evidence

| Readout | Knockdown/control estimate | Statistical evidence | Interpretation |
|---|---:|---|---|
| NAA25 | 0.442, about 56% lower | Exact one-sided p = 0.00794; primary Holm p = 0.02381 | Manipulation check passes |
| SYP | 0.404, about 60% lower | Proteome-wide BH FDR = 0.02389 | Supporting synaptic-protein result |
| SYT1 | 0.431, about 57% lower | Proteome-wide BH FDR = 0.03427 | Supporting synaptic-protein result |
| SYT1, complete and equally weighted backbones | 0.368, about 63% lower | BH FDR = 0.03753 across 2,183 proteins | Strongest coverage-robust result |
| Fixed neuronal panel, 7 measured genes | 0.697, about 30% lower | Welch p = 0.00263; exploratory exact two-sided p = 0.00794 | Broader neuronal-protein response |
| Fixed proteasome-assembly panel, 39 measured genes | 0.999 | 90% CI for ratio about 0.898–1.113, within the fixed 0.8–1.25 margin | Relative abundance preservation; not measured enzyme activity |

Ratios are powers of estimated differences on the log2 scale, not ratios of raw arithmetic protein means. Protein-level FDR values are from two-sided Welch tests corrected over the whole eligible proteome. The original analysis tested 4,511 proteins, of which 520 met FDR <0.05, including 323 decreases. SYP and SYT1 are therefore not the only changing proteins, and these data do not establish a uniquely synaptic effect.

The primary SYP log2 effect is −1.306 (95% CI −1.759 to −0.853); SYT1 is −1.214 (95% CI −1.793 to −0.635). See `proteome_effects.tsv` for full precision and every protein, including null and opposite-direction results.

## The original selectivity hypothesis did not pass

| Frozen primary contrast | Log2 effect | Exact one-sided p | Holm p |
|---|---:|---:|---:|
| NAA25 manipulation | −1.1776 | 0.00794 | 0.02381 |
| SNCA relative to proteasome panel | −0.8036 | 0.05556 | 0.11111 |
| SNCA relative to neuronal panel | −0.3938 | 0.15079 | 0.15079 |

Both relative-effect point estimates exceed the prespecified reduction margin, but neither corrected selectivity p value passes, and the neuronal panel is not preserved. The source measures SNCA through one internal peptide with five control and four knockdown observations. Its absolute abundance estimate is lower, but this reanalysis alone does not establish a statistically secure SNCA decrease. The published immunoblot and other assays remain separate evidence.

No missing SNCA quantity was imputed. A two-peptide SNCA analysis is `NOT_EVALUABLE`. This limitation is retained rather than solved by counting charge states or peptides as biological replicates.

## Robustness and unsuccessful checks

- **Peptide weighting:** after collapsing charge/modification features to equal-weight peptide backbones, SYP remains lower (ratio 0.348, FDR 0.01399) and SYT1 remains lower (ratio 0.422, FDR 0.03561).
- **No missing measurements:** SYT1 remains lower using only precursors quantified in all ten runs and requiring at least two backbones per protein. All 14 complete SYT1 backbones point down. These are technical measurements of the same cultures, not 14 independent experiments.
- **SYP strict-coverage limitation:** only one SYP backbone is complete in all ten runs, so the strict two-backbone SYP analysis is not evaluable. Its other four eligible backbones all point down, but this does not remove the strict-coverage limitation.
- **Every sample omission:** SYP effects range from −1.413 to −1.185 log2, SYT1 from −1.348 to −1.057, and the neuronal panel from −0.584 to −0.459. Every effect stays negative. Omitting the low-coverage NAA25_07 run also retains all three directions.
- **Not solely the highlighted pair:** after dropping both SYP and SYT1, the remaining five-gene panel is about 16% lower (log2 effect −0.250, 95% CI −0.400 to −0.099). This is a post-hoc sensitivity, not a new confirmatory endpoint.
- **Panel comparison:** neuronal-minus-proteasome abundance is lower with both median and mean aggregation. The median contrast is −0.520 log2 (95% CI −0.876 to −0.164; exploratory exact two-sided p = 0.00794).
- **Normalization:** the frozen SNCA-relative contrasts are mathematically invariant to a common sample loading offset; the implementation verifies this to 1e−12.

All original and alternative protein tables, not only the highlighted results, are included. Repeated significance across these correlated analyses is evidence of computational robustness, not independent biological replication.

## Mechanistic interpretation

The intervention targets NAA25, an auxiliary NatB component, rather than the original eight proteasome leads. The useful connection is a substrate-clearance strategy that operates through the proteasome. The new analysis asks whether retaining bulk proteasome protein abundance is enough to retain a neuronal protein program. In this experiment it is not: the relative-abundance panel is preserved while multiple neuronal proteins fall.

SYP and SYT1 have canonical N termini beginning ML and MV, outside the common NatB recognition motif MD/ME/MN/MQ. This makes a simple canonical-substrate explanation incomplete, but it does not prove that either effect is indirect. N-terminal acetylation and degradation rates were not measured here.

Several explanations remain consistent with the data: altered neuronal maturation, altered synaptic protein synthesis or stability, downstream effects of alpha-synuclein loss, NatB-independent NAA25 effects, or another aspect of the perturbation protocol. The experiment does not isolate those pathways. Viability and electrophysiological function cannot be inferred from relative protein abundance, and MAP2 is not decreased in the original analysis, arguing against treating the panel as a direct neuron-count assay.

The existing Tian survival-screen phenotypes are recorded separately in `natb_existing_neuronal_survival.tsv`. They use different outcomes and cannot validate this synaptic-protein effect. Their rank percentiles are not p values. No independent synaptic replication was found in the qualified inputs; a failed metadata qualification is not a biological null result.

The most direct future discriminating test would measure SYP/SYT1 and synaptic function in a NAA25-by-UBE2W intervention, alongside alpha-synuclein. Whether blocking the known clearance route also restores synaptic proteins would separate shared-pathway from parallel effects. That experiment has not been performed by this project.

## Relationship to the original project and AAN

The earlier proteasome convergence, mixed human validation, failed ciliary validation and unsuccessful prediction benchmarks retain their original labels. This result does not retroactively validate those hypotheses. The project now has an intervention-based molecular selectivity finding, supported by multiple peptides and sample-level analysis, rather than only observational pathway overlap.

A defensible central claim is: **Reanalysis of a published alpha-synuclein clearance intervention identifies a robust decrease in SYT1 and a broader neuronal-protein response despite preserved relative proteasome-assembly abundance, revealing a selectivity constraint that warrants functional testing.**

This is a substantive mechanistic addition with a specific experimental consequence. It remains a single-study computational reanalysis. Award competitiveness cannot be quantified from these data, and a replicated new causal mechanism has not been established. Do not claim that we discovered NatB-mediated alpha-synuclein clearance, proved synaptic toxicity, or developed a validated therapy.

## Reproduce

From the repository root, using CPython 3.12.14 and the tested environment:

```bash
python -m pip install -r docs/v2/natb_selectivity/environment_freeze.txt
python scripts/run_natb_selectivity.py --fetch
```

For verification of checked-in artifacts only:

```bash
python scripts/run_natb_selectivity.py --verify-only
```

The runner checks the hashes of code, specifications and existing dependencies, acquires the four source CSVs when absent, verifies their individual SHA256 values, runs the primary analysis and post-hoc checks, renders both figures, and checks the output manifest. The outer Europe PMC ZIP has variable wrapper metadata; frozen hashes apply to its inner source files. Raw source CSVs are kept under ignored `data_raw/natb_selectivity/`. No source credentials are required.

`environment_freeze.txt` is the full tested shared analysis environment, including dependencies also used by other extensions. `requirements.txt` records this extension's direct dependencies. These pins and observed execution do not constitute an independently rebuilt environment from package wheels.

## Figures and audit trail

- `fig1_intervention_selectivity.png` and `.pdf`: protein effects with 95% Welch intervals and the within-sample neuronal/proteasome contrast. Dots in panel B are culture samples. The panel p value is exploratory.
- `fig2_synaptic_robustness.png` and `.pdf`: alternative protein aggregation and all complete SYT1 peptide effects. The failed SYP coverage criterion is displayed explicitly.
- `primary_contrasts.tsv`, `module_equivalence.tsv`, `validation_summary.json`: frozen endpoint and decision.
- `proteome_effects.tsv` and both `proteome_equal_backbones_*.tsv` files: all eligible proteins and multiplicity corrections.
- `sample_qc.tsv`, `gene_coverage.tsv`, `protein_sample_scores.tsv`, precursor/backbone tables: measurement and sample audit.
- `normalization_sensitivities.tsv`, `leave_one_run_out.tsv`, `synaptic_leave_one_run_out.tsv`, `neuronal_leave_gene_out.tsv`: every sensitivity, including weak and unevaluable cases.
- `shared8_*`, `screen_coverage.json`, `natb_existing_neuronal_survival.tsv`: original-lead coverage and already existing contextual results. None is a new independent confirmation.
- `docs/v2/natb_selectivity/analysis_plan.md`, `robustness_addendum.md`, `source_audit.md`: provenance, analysis order and alternative-source qualification.
- `SHA256SUMS.txt` and `docs/v2/natb_selectivity/verification.json`: release integrity and the exact isolated reproduction scope.

This standalone extension does not close the separate whole-project R 4.5.2 clean-clone execution gate.
