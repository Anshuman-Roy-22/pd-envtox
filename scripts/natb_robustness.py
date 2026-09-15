#!/usr/bin/env python3
"""Post-hoc, fully reported checks of the neuronal-protein liability signal."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from natb_selectivity import ROOT, DOCS, RUNS, GROUP, acquire, adjust, difference, dump, prepare, save


def two_sided_exact(values):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('Module permutation requires complete run-level values')
    observed = values[GROUP].mean() - values[~GROUP].mean()
    null = []
    for ix in itertools.combinations(range(10), 5):
        mask = np.zeros(10, dtype=bool); mask[list(ix)] = True
        null.append(values[mask].mean() - values[~mask].mean())
    return float(np.mean(np.abs(null) >= abs(observed) - 1e-12))


def balanced_proteome(features, complete):
    f = features[features[RUNS].notna().all(axis=1)].copy() if complete else features.copy()
    peptide = f.groupby(['gene', 'peptide_backbone'])[RUNS].median()
    gene = peptide.groupby('gene').median()
    counts = peptide.groupby('gene').size()
    gene = gene.loc[counts[counts >= 2].index]
    gene = gene.sub(gene.median(axis=1), axis=0)
    reference = gene.index[gene.notna().all(axis=1)].difference(['NAA25', 'SNCA'])
    offset = gene.loc[reference].median(axis=0)
    norm = gene.sub(offset, axis=1)
    rows = []
    for g, v in norm.iterrows():
        if min(v.iloc[:5].notna().sum(), v.iloc[5:].notna().sum()) < 4:
            continue
        rows.append(dict(gene=g, peptide_backbones=int(counts[g]), **difference(v)))
    table = pd.DataFrame(rows)
    table['FDR_BH'] = adjust(table.p_two_sided, 'BH')
    return table, norm, peptide, offset


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, default=ROOT / 'data_raw/natb_selectivity')
    ap.add_argument('--output', type=Path, default=ROOT / 'results/v2/natb_selectivity')
    ap.add_argument('--fetch', action='store_true')
    a = ap.parse_args(); o = a.output; o.mkdir(parents=True, exist_ok=True)
    acquire(a.cache, a.fetch)
    raw, norm, coverage, loading, features = prepare(a.cache / 'adj4767_Table_S5.csv', o, return_features=True)
    sets = json.loads((DOCS / 'gene_sets.json').read_text())
    membership = json.loads((o / 'module_membership.json').read_text())
    targets = list(dict.fromkeys(['SNCA', 'NAA25', 'SNCG'] + sets['neuronal_identity']))
    selected = []; backbone_rows = []; qc = []
    for complete, label in [(False, 'equal_backbones_at_least_8_runs'), (True, 'equal_backbones_complete_10_runs')]:
        table, balanced, peptides, offset = balanced_proteome(features, complete)
        save(table, o / f'proteome_{label}.tsv')
        subset = table[table.gene.isin(targets)].copy(); subset.insert(0, 'analysis', label)
        selected.append(subset)
        qc.append(dict(analysis=label, tested_genes=len(table), complete=complete,
                       requires_two_backbones=True, SNCA_evaluable='SNCA' in table.gene.values))
        if complete:
            save(balanced.loc[[g for g in targets if g in balanced.index]].rename_axis('gene').reset_index(),
                 o / 'complete_backbone_sample_scores.tsv')
        for gene in ['SYP', 'SYT1']:
            if gene not in peptides.index.get_level_values('gene'):
                continue
            for peptide, v in peptides.loc[gene].iterrows():
                z = v - offset
                if min(z.iloc[:5].notna().sum(), z.iloc[5:].notna().sum()) < 2:
                    continue
                backbone_rows.append(dict(analysis=label, gene=gene, peptide_backbone=peptide,
                                          **difference(z)))
    save(pd.concat(selected, ignore_index=True), o / 'backbone_target_effects.tsv')
    backbone = pd.DataFrame(backbone_rows); save(backbone, o / 'synaptic_backbone_effects.tsv')
    save(features[features.gene.isin(['SYP', 'SYT1'])], o / 'synaptic_precursor_scores.tsv')
    modules = {}
    module_rows = []
    for aggregation in ['median', 'mean']:
        s = {}
        for name in ['neuronal_identity', 'proteasome_assembly']:
            s[name] = getattr(norm.loc[membership[name]], aggregation)(axis=0).to_numpy()
        s['neuronal_minus_proteasome'] = s['neuronal_identity'] - s['proteasome_assembly']
        for name, v in s.items():
            module_rows.append(dict(aggregation=aggregation, contrast=name, **difference(v),
                                    p_exact_two_sided=two_sided_exact(v)))
        if aggregation == 'median': modules = s
    save(pd.DataFrame(module_rows), o / 'module_robustness.tsv')
    omission = []
    vectors = {g: norm.loc[g].to_numpy() for g in ['SYP', 'SYT1']}
    vectors.update(modules)
    for i, run in enumerate(RUNS):
        keep = np.arange(10) != i
        for name, v in vectors.items():
            omission.append(dict(excluded_run=run, contrast=name, **difference(v[keep], GROUP[keep])))
    save(pd.DataFrame(omission), o / 'synaptic_leave_one_run_out.tsv')
    gene_omission = []
    drop_sets = [[g] for g in membership['neuronal_identity']] + [['SYP', 'SYT1']]
    for drops in drop_sets:
        left = [g for g in membership['neuronal_identity'] if g not in drops]
        v = norm.loc[left].median(axis=0).to_numpy()
        gene_omission.append(dict(excluded_genes=','.join(drops), retained_genes=','.join(left), **difference(v)))
    save(pd.DataFrame(gene_omission), o / 'neuronal_leave_gene_out.tsv')
    survival = []
    for mode in ['CRISPRi', 'CRISPRa']:
        d = pd.read_csv(ROOT / f'results/v2/neuronal_survival/{mode}_gene_level.tsv', sep='\t')
        survival.append(d[d.gene.isin(['NAA25', 'NAA20', 'METAP2'])])
    save(pd.concat(survival, ignore_index=True), o / 'natb_existing_neuronal_survival.tsv')
    # These canonical sequence annotations are a stored source snapshot, not predicted acetylation.
    seq = pd.read_csv(DOCS / 'sequence_annotation.tsv', sep='\t', dtype=str)
    seq['canonical_NatB_motif'] = seq.N_terminus.str[:2].isin(['MD', 'ME', 'MN', 'MQ'])
    save(seq, o / 'synaptic_sequence_annotation.tsv')
    summary = dict(status='EXPLORATORY_FOLLOWUP', proteome_sensitivity_coverage=qc,
                   backbone_direction_counts=backbone.assign(down=backbone.log2_difference < 0)
                   .groupby(['analysis', 'gene']).agg(backbones=('down', 'size'), down=('down', 'sum'))
                   .reset_index().to_dict('records'),
                   independent_synaptic_function_replication=False,
                   sample_level_neurotransmission_measured=False,
                   primary_classification_unchanged=True)
    dump(o / 'robustness_summary.json', summary)
    print(pd.concat(selected, ignore_index=True).to_string(index=False))
    print(pd.DataFrame(module_rows).to_string(index=False))
    print(pd.DataFrame(gene_omission).to_string(index=False))


if __name__ == '__main__':
    main()
