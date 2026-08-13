# Limma-confirmed Reactome convergence

Reactome enrichment was rerun from limma-moderated GSE46798 and GSE17542 t
statistics and a directional Stouffer meta-z built from limma-moderated human
postmortem results.

Five pathways satisfy the strict rule of a common NES direction and nominal
p < 0.05 in all three model classes. Proteasome assembly is the cleanest
PD-relevant result:

| Model | NES | nominal p |
| --- | ---: | ---: |
| Human corrected DA neurons, paraquat/maneb | -1.496 | 0.0237 |
| Mouse SN DA neurons, 10-day MPTP | -1.422 | 0.0383 |
| Human PD postmortem SN meta-analysis | -3.155 | <0.001 |

Cross-model meta-FDR is 0.00203. Eight leading-edge genes recur in all three:
`POMP`, `PSMA4`, `PSMB4`, `PSMC5`, `PSMD1`, `PSMD2`, `PSMD4`, and `PSMG1`.

Neither toxin dataset has any Reactome pathway passing within-dataset FDR. This
supports convergent pathway evidence across models, not independently confirmed
pesticide-induced proteasome suppression.
