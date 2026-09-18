# Paper-focused overview of the analytical overhaul

Updated 18 September 2026. Read alongside [the reconstructed audit and new matched-perturbation comparison](reconstruction_20260918/README.md). This file is an evidence map for the student's own paper, not a claim that every result is independent or prospective.

The [subsequent specificity follow-up](specificity_followup_20260918/README.md) adds exact permutation-inverted human intervals, a CAMERA correlation curve and a descriptive dependency plot. All four human intervals remain below zero under the stated shift model. CAMERA's estimated mode remains nonsignificant and uses different denominator degrees of freedom from fixed-correlation mode. The global-disruption/realized-knockdown regression remains **blocked by missing broad expression inputs** and has not been performed. The dependency plot cannot substitute for that check.

## One scientific question

**Do environmental-toxicant transcriptomic signatures nominate protein-homeostasis processes whose genetic perturbation alters neuronal synaptic programs, and how specific are those responses?**

The starting project asked whether MPTP and maneb/paraquat altered network-prioritized PD genes. The earlier PPI work supplies the origin of that candidate panel; describe and cite it proportionately, distinguishing prior work from the present study. The overhaul replaces weak directional-concordance claims with formal tests, establishes the original panel's limits, examines broader pathway convergence and evaluates a biologically motivated perturbation bridge. Its argument is nomination, generalizability and experimental-data triangulation. It should not imply that each exploratory decision was planned at the project's beginning.

## Changes that matter to the paper

| Analytical change | Main result | Scientific role |
|---|---|---|
| Corrected exposure/sample contrasts, fixed probe/gene mappings, formal limma models, multiple-testing correction and matched-background tests for the original network panel | The panel did not establish robust external generalization | Establishes the limit of network-based nomination and corrects the earlier interpretation of concordance |
| Reactome analyses from formal transcriptomic statistics across toxicant and human postmortem models | Proteasome assembly has cross-model meta-FDR 0.00203; eight shared leading-edge genes | Nominates a protein-homeostasis process; neither toxin dataset individually passes pathway FDR |
| Frozen rotenone and held-out toxicant testing | Primary thresholds were not met; ciliary/Hedgehog validation remained NOT_CONFIRMED; aggrephagy gate stayed closed | Prevents selective retention of favorable follow-ups |
| Donor-level human substantia nigra analysis | One cohort met the prespecified one-sided threshold (p=0.0352); independent cohort was null (p=0.5095) | Human generalization was not independently replicated; significance differs across cohorts without proving formal heterogeneity |
| Published neuronal-survival and rescue-related follow-ups | Eight-gene survival primary p=0.14558; other unsupported/unevaluable classifications retained | Provides boundary tests rather than additional confirmation |
| Human NAA25-knockdown proteomic reanalysis | SYP and SYT1 decrease; relative synaptic scores remain lower than general and matched backgrounds | Shows a synaptic protein response to a direct NatB perturbation, with important competitive-enrichment and cell-health limits |
| Mouse genetic-perturbation reanalysis | Naa20, Pomp, Psmb4 and Psmc5 reduce the frozen synaptic RNA score relative to non-targeting neurons | Links perturbations of NatB and selected proteasome machinery to a common transcriptional readout using animals as replicates |
| New locked comparison against matched molecular perturbations | Naa20, Pomp and Psmc5 retain larger decreases; Psmb4 is unevaluable | Addresses an important generic-perturbation alternative within defined comparator panels |
| Reproducibility and reporting safeguards | Fixed plans, source/selection checksums, complete outcomes, animal/culture-level tables, deterministic figures and branch history | Makes the student's analytical contribution inspectable and prevents null findings disappearing from the narrative |

## The bridge between the eight leads and NatB

The recurring leads are **POMP, PSMA4, PSMB4, PSMC5, PSMD1, PSMD2, PSMD4 and PSMG1**. They nominate proteasome-related biology; they are not themselves evidence that NatB regulates all eight.

NatB became relevant through published experimental work connecting N-terminal acetylation to α-synuclein stability and clearance. [Santhosh Kumar et al.](https://www.science.org/doi/10.1126/sciadv.adj4767) provide that biological rationale and the human NAA25-knockdown experiments. Credit their mechanistic findings to them. The project then asks what additional neuronal protein changes accompany that perturbation, and whether NatB/proteasome perturbations share a synaptic response in another experimental source.

The mouse atlas supplies Naa20, Pomp, Psmb4 and Psmc5. It does not supply every original lead or the same human NAA25 perturbation. Thus, three of eight leads have the reconstructed non-targeting comparison; two have evaluable matched-perturbation specificity evidence. Naa20 and NAA25 are different subunits of the same NatB complex, studied in different species and readouts. This is triangulation across perturbations, not a demonstrated eight-gene-to-NatB regulatory chain.

## New specificity result to feature, with its limits

The locked mouse comparison finds additional mean synaptic-score differences of -0.01219 for Naa20 (38 mice, BH4=0.000540), -0.03517 for Pomp (19 mice, BH4=4.34e-6), and -0.03953 for Psmc5 (25 mice, BH4=1.57e-6) relative to their matched comparator perturbations. All three retain significance after litter averaging. Psmb4's single qualifying comparator is below the fixed minimum, so its specificity is not established.

Controls match external neuronal-dependency percentile, non-targeting baseline expression and mouse coverage. They do not match every consequence of essential-gene disruption or directly measure cell health. Naa20's effect is small in mean log1p score units. Always report its magnitude alongside significance.

In human proteomics, the four synaptic-minus-background comparisons pass exact-permutation BH4=0.0159. However, 994/4,511 proteins change at FDR<0.05, multiple functional modules decrease, and estimated-correlation CAMERA does not establish competitive synaptic enrichment (FDR=0.627). The supported claim is a reproducible relative synaptic response, not exclusive synaptic selectivity.

## Suggested paper argument

1. **Nomination and its limits:** explain why environmental transcriptomics can test the earlier network candidates; present the panel's formal failure and the broader proteasome convergence without preserving the old overclaim.
2. **Generalizability:** report the frozen toxicant and human tests together, including failed replication. This defines how far the nomination travels.
3. **Perturbation evidence:** explain the published NatB rationale before presenting the human protein and mouse RNA analyses. Distinguish source authors' experimental discoveries from the student's analytical tests.
4. **Specificity and unresolved mechanism:** present relative backgrounds and matched perturbations with their null/unevaluable components. State the remaining functional and causal gap directly.

At the start of Methods, distinguish retrospective discovery, outcome-blind choices made within an already selected source, frozen held-out tests, and the September reconstruction. Internal pre-fit commits improve auditability; they do not make all stages prospectively preregistered. The reconstructed analyses had known outcomes. The new comparator manifest was locked before comparator synaptic fitting.

## Findings that must remain visible

- Human PD replication failed in the independent cohort. A one-sided primary threshold with a two-sided interval crossing zero should not be described as broadly replicated disease validation.
- Proteasome RNA increases under all four mouse perturbations, whereas the motivating toxin-associated pathway signal decreased. Compensation is plausible but untested.
- Naa20 perturbation increases Snca RNA in the mouse data; this is not replication of the human α-synuclein protein decrease.
- Neither relative protein scores nor matched RNA responses measure proteasome activity, acetylation, synaptic physiology or rescue. No direct toxin-to-NatB mediation or therapy is established.
- The mouse atlas is an unreviewed preprint and a young-mouse brain experiment. General proteasome/synaptic observations already appear in the source work; the new matched comparison must not be credited as the original experiment.

## Figure plan

Keep the earlier five V2 figures as the historical evidence package. The new `specificity_comparison` figure supplies an audited two-panel summary of human relative scores and mouse matched controls. For final paper layout, use an evidence map that marks which relationships are measured, a nomination/validation figure including nulls, a human protein figure with culture points and both CAMERA settings, and a mouse figure with animal/litter points, comparator balance and missing Psmb4 support. Include a complete outcome matrix and source-selection/provenance table in the supplement. Do not draw an established causal arrow from toxins through NatB to synaptic failure.

The remaining highest-value work is presentation, independent methodological review and completion of the repository-wide reproduction gate. Further open-ended dataset searching is not needed to strengthen the present claim; direct functional or rescue evidence would be required for a substantially stronger causal claim.
