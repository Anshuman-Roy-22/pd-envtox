# Signed Reactome convergence

## Primary conclusion

The full Reactome meta-analysis identifies broad cross-model negative pathway
shifts, but 328 significant overlapping pathways are too redundant to interpret
individually. A stricter recurrence rule required the same NES direction and
nominal p < 0.05 in all three model classes. Fourteen pathways pass.

The strongest PD-relevant result is **Proteasome assembly**:

| Model | NES | nominal p |
| --- | ---: | ---: |
| Human corrected DA neurons, paraquat/maneb | -1.546 | 0.0126 |
| Mouse SN DA neurons, 10-day MPTP | -1.475 | 0.0305 |
| Human PD postmortem SN meta-analysis | -3.231 | <0.001 |

Cross-model signed Stouffer FDR is 0.000512. Nine leading-edge genes recur in
all three rankings: `POMP`, `PSMA4`, `PSMB4`, `PSMC5`, `PSMD1`, `PSMD2`,
`PSMD4`, `PSME3`, and `PSMG1`.

Other strict recurring processes include kinesin transport, aggrephagy, and RAB
GEF regulation. Their shared leading edges are less pathway-specific and often
contain cytoskeletal genes, so proteasome assembly is the cleanest mechanistic
lead.

## Interpretation boundary

- GSE46798 has zero Reactome pathways at within-dataset FDR <= 0.05.
- GSE17542 has two; human postmortem has 389.
- Therefore this is evidence of **cross-model pathway convergence**, not proof
  that pesticide exposure independently establishes proteasome suppression.
- Reactome pathways overlap extensively; the 14 strict results are related, not
  14 independent mechanisms.
- Final release should rerun enrichment from limma-moderated rankings and test
  the proteasome result in an additional replicated toxin dataset if feasible.

This pathway result is independent of the failed 21-gene validation and supports
a revised systems-level project more strongly than an individual-candidate claim.
