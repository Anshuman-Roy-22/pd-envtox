# Exploratory follow-up after the frozen selectivity result

Written 2026-09-15 after the primary results were known. This addendum is not a new preregistration and does not change the frozen `SELECTIVITY_NOT_ESTABLISHED` classification.

The original analysis found lower SYP and SYT1 abundance at proteome-wide BH FDR <0.05 and lower abundance of the fixed neuronal panel. The question is now whether this signal survives defensible changes in peptide coverage and aggregation, and what mechanisms it can distinguish.

## Complete analysis, including unsuccessful checks

1. Re-estimate every measurable protein using only precursors valid in all ten runs. Aggregate charge/modification features to peptide backbones first, then give each backbone equal weight. Require at least two backbones per protein; estimate loading offsets from this complete set, excluding NAA25 and SNCA. Correct the entire eligible proteome with BH. This avoids both changing peptide membership across samples and dominance by multiple charge states. SNCA's known single, incomplete peptide cannot pass this quality criterion.
2. Also report an equal-backbone analysis retaining the original >=8/10 feature rule. This separates completeness from aggregation changes.
3. For SYP and SYT1, report every eligible backbone's effect and count how many point down. Peptides are technical features, not independent replicates. Do not calculate a binomial significance test over peptides.
4. Report all single-run omissions for SYP, SYT1, the neuronal panel, the proteasome panel and their within-run difference, using the primary fixed protein scores. NAA25_07 has substantially fewer detected precursors; report its omission explicitly among all ten omissions, without deleting it from the primary result.
5. Use both median and mean versions of the original, fixed modules and calculate their neuronal-minus-proteasome contrast. Give exact two-sided sample-label permutation p values; these are exploratory and are not separately confirmatory.
6. Drop each neuronal-panel gene in turn, and drop SYP and SYT1 together. Report effect sizes and intervals to determine how much the panel finding depends on the two highlighted proteins. Do not promote a favorable sensitivity over an unfavorable one.
7. Record SYP and SYT1 canonical N termini from UniProt. The MD/ME/MN/MQ rule is a motif annotation, not a measurement of acetylation or proof that a protein is an indirect target. Match canonical accessions and record the isoform ambiguity in the source MS protein groups.
8. Extract NAA25, NAA20 and METAP2 from the already-analyzed Tian neuronal survival data, keeping the phenotype, condition and percentile labels. These are orthogonal context, not replication of a synaptic phenotype. Percentiles are not p values, and a lack of depletion does not establish neuronal function or safety.

## Interpretation boundary

A reproducible SYP/SYT1 abundance decrease supports a synaptic-protein liability hypothesis in this NAA25 intervention experiment. It cannot establish defective neurotransmitter release, neuron death, altered maturation, direct substrate degradation, SNCA-mediated effects, or the same effect after pharmacological METAP2 inhibition. The experiment does not separate those possibilities. The previously published proteasome/UBE2W dependencies for alpha-synuclein belong to the source authors.

No causal mediation analysis is appropriate with these ten samples and unmeasured competing pathways. No gene-set or effect threshold will be optimized to achieve a positive label.
