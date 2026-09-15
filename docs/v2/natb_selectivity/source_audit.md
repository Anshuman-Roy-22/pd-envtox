# Causal-source audit and route selection

This audit responds to the request to execute alternatives autonomously. Source selection used experimental design, relation to the prior proteasome lead and availability, rather than new computed p values.

| Route | What was verified | Disposition |
|---|---|---|
| Bi et al. 2021, DOI 10.1016/j.redox.2021.102167 | POMP/NRF2 and immunoproteasome interventions, PLK2/alpha-synuclein and dopaminergic outcomes; full text and supplementary archive acquired | Causal prior directly relevant to original lead POMP. Public supplement is a document containing figures, not an identified replicate-level outcome table. Do not invent new sample-level estimates. |
| Giandomenico et al. 2025, DOI 10.1101/2025.04.25.650603 | Full preprint recovered from University of Basel; POMP-versus-PSMB perturbations distinguish proteasome impairment and a nucleolar response | Attractive candidate-specific contrast, but no accession or downloadable replicate-level count matrix identified in the inspected manuscript/repository page. It remains a preprint. No raw-data analysis claimed. |
| Santhosh Kumar et al. 2024, DOI 10.1126/sciadv.adj4767 | Complete S1-S5 archive recovered from Europe PMC; ten explicitly named proteomic sample columns and peptide Q values; published dependency interventions | Selected for a new intervention-selectivity analysis. One-peptide SNCA coverage and missing POMP/UBE2W quantities are recorded before fitting. |
| Woo et al. 2025, DOI 10.1016/j.cell.2025.05.029 | GEO GSE279706 has sustained interferon exposure; GSE279707 has sorted nuclei in EAE | These verified transcriptome designs do not themselves provide the required perturbation-by-dependency contrast. The publication also has 2025 and 2026 corrections that would need reconciliation. No expression outcomes examined or claimed. |
| PSMG1 sex-specific genetics | Search surfaced a recent title and candidate mention, without a verified complete analysis source in this audit | Not counted as established genetic evidence or used to select an analysis result. |

The original eight leads are POMP, PSMA4, PSMB4, PSMC5, PSMD1, PSMD2, PSMD4 and PSMG1. The new intervention targets NAA25, not those eight. This must remain clear: it tests a possible substrate-clearance strategy related to the proteasome, not genetic rescue of the original eight-gene signature.

Retrieval issues were resolved through public repositories where possible. PMC and publisher direct supplement endpoints returned HTML/403; Europe PMC's public supplementary-files endpoint supplied the complete NatB archive. The preprint server returned 429; the University of Basel's public repository supplied the deposited PDF. No private datasets or unpublished results were accessed, and no authors were contacted.

Sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC8577461/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10857481/ ; https://edoc.unibas.ch/entities/publication/be0e3311-7de4-4f8a-a7be-83ae9e15c9a7/full ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279706 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279707 .

## Follow-up qualification after the selectivity result

The synaptic-protein signal prompted a bounded check for independent perturbation evidence. These searches followed inspection of the NatB result and are not an outcome-blind validation search.

- The existing Tian 2021 neuronal survival tables contain NAA25, NAA20 and METAP2. Their screen phenotypes and percentile ranks can provide orthogonal context, but they do not measure synaptic protein abundance or activity. The `source_rows` and `unique_tss` columns are counts, not p values.
- Replogle et al. 2022 includes NAA25 Perturb-seq in K562 cells. That cell type cannot independently test loss of a neuronal synaptic program, so it was not substituted for a neuronal validation experiment.
- Stephan et al. 2012 (DOI 10.1523/JNEUROSCI.3116-12.2012) concerns fly Psidin, neuronal survival and axon targeting. It is relevant prior biology, but no new replicate-level validation was performed. It is not evidence that human SYP/SYT1 loss has the same cause.
- A 2026 preprint (DOI 10.64898/2026.01.16.699530; PMC12871359) reports dopamine-neuron CRISPR effects under alpha-synuclein stress. Its linked public repository, `Arinze-BioX/daichi_etal_parkinsons`, contains guide-count files and alignment scripts. However, its inspected tree has no sample manifest identifying which trial/amplicon files belong to the alpha-synuclein versus untreated conditions. The preprint's text and Figure 6 legend also report differing animal counts and histology statistics. No sample mapping was inferred from A/D filenames, no count outcomes were fitted, and no independent replication was claimed. The article has preprint status in the retrieved record.
- UniProt canonical SYP P08247 and SYT1 P21579 sequences were retrieved on 2026-09-15. S5 lists those exact accessions. Their first two residues are ML and MV, respectively. These fall outside MD/ME/MN/MQ; the annotation does not measure acetylation, identify all possible isoforms, or prove an indirect response. The extracted source snapshot and response hashes are stored in `sequence_annotation.tsv`.

No newly trained virtual cell or protein-structure model was substituted for missing intervention outcomes. The implemented alternative is the fully reported experimental-proteomics robustness analysis in `robustness_addendum.md`.
