# Continuous network-predictive test

The fixed 21-gene cutoff could have failed even if network proximity retained a
weak genome-wide relationship with transcriptional response. This secondary,
post-null test evaluated that possibility across all eligible non-seed genes.

Unadjusted network closeness was unrelated to absolute response in all four
prespecified outcomes. After rank adjustment for STRING degree and GTEx SN TPM,
only corrected-neuron pesticide response in GSE46798 was statistically nonzero:

- partial Spearman r = 0.03075
- p = 0.00106
- FDR across four outcomes = 0.00425

The effect explains less than 0.1% of rank variation and did not replicate in
SN MPTP response, SN-versus-VTA interaction, or human postmortem PD. It is not
evidence of practically useful prediction and should not be used to rescue the
fixed panel. The defensible conclusion is that this minimum-distance network
score has negligible and non-generalizing transcriptomic predictive value.
