# GSE46798 validation result

The prespecified design is in `docs/v2/gse46798_frozen_analysis_plan.md`.

## Main result

- All 21 fixed genes were evaluable.
- Neither pesticide exposure contrast nor the genotype-by-exposure interaction
  produced an FDR-significant fixed-panel gene.
- No fixed module or the full panel was enriched for response magnitude or a
  coordinated signed response after comparison with expression-matched genes.
- The primary corrected-cell exposure contrast therefore does not independently
  validate the network panel as a pesticide-responsive set.

## Secondary genotype observation

Under vehicle conditions, three panel genes were higher in SNCA A53T neurons at
genome-wide FDR <= 0.05:

| Gene | log2FC | p | genome-wide FDR |
| --- | ---: | ---: | ---: |
| SNX2 | 1.424 | 0.000169 | 0.00726 |
| WLS | 1.054 | 0.000414 | 0.0117 |
| VPS26A | 0.919 | 0.00303 | 0.0376 |

This was a prespecified context contrast, but it is not evidence of pesticide
responsiveness. WLS is directionally consistent across all five unambiguous
probes. VPS26A is positive across both probes, with unequal strength. SNX2 has
one unambiguous annotated probe. The retromer module was not enriched against
the matched background, so these genes must be reported individually rather
than as a confirmed module-level effect.

## Interpretation boundary

These ordinary factorial-model results require empirical-Bayes limma confirmation
before competition release. The balanced design makes coefficient estimates
stable, but moderated uncertainty and final genome-wide FDR may change. The
large background genotype signal also means the A53T observations should be
treated as independent leads, not proof that the network method broadly worked.
