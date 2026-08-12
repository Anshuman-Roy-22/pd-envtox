# Frozen Reactome convergence plan

Frozen before pathway enrichment was calculated.

## Gene-set source

Official Reactome `ReactomePathways.gmt`, downloaded August 12, 2026 from the
Reactome current-release endpoint. Analyze pathways with 10–500 measured genes.

## Primary signed rankings

1. GSE46798 pesticide effect in genetically corrected human DA neurons: t.
2. GSE17542 10-day MPTP effect in SN DA neurons: t.
3. Human PD postmortem directional meta-analysis: meta-z.

Run preranked GSEA with 1,000 phenotype-independent gene-set permutations and a
fixed seed. Benjamini-Hochberg FDR is calculated within each dataset.

Cross-model signed Stouffer analysis combines nominal pathway p-values and NES
directions across the three independent model classes. Replace permutation
p = 0 by 1/1001. Apply FDR across all pathways common to all three results.
Cross-model significance requires meta-FDR <= 0.05; recurrence and leading-edge
genes must also be shown. Opposite NES directions cancel rather than being
recast as convergence.

## Prespecified targeted family

- Aerobic respiration and respiratory electron transport
- Respiratory electron transport
- Complex I biogenesis
- Mitochondrial biogenesis
- Cellular response to mitochondrial stress
- Mitophagy
- PINK1-PRKN Mediated Mitophagy
- Autophagy
- Macroautophagy
- Chaperone Mediated Autophagy
- Lysosome Vesicle Biogenesis
- Endosomal Sorting Complex Required For Transport (ESCRT)
- Vesicle-mediated transport
- Golgi Associated Vesicle Biogenesis
- Protein ubiquitination
- Proteasome assembly
- Neurotransmitter release cycle
- Dopamine Neurotransmitter Release Cycle
- Innate Immune System
- Cytokine Signaling in Immune system

The complete Reactome screen remains primary for correction. The targeted table
is an interpretation aid, not a separately corrected discovery universe.
