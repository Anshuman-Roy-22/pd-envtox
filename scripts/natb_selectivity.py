#!/usr/bin/env python3
"""Frozen culture-level proteomic selectivity reanalysis of public Table S5."""
import argparse
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import re
from urllib.request import Request, urlopen
import zipfile

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs/v2/natb_selectivity'
RUNS = [f'NTC_{i:02d}' for i in range(1, 6)] + [f'NAA25_{i:02d}' for i in range(6, 11)]
GROUP = np.array([False] * 5 + [True] * 5)
MARGIN = float(np.log2(1.25))


def dump(path, data):
    def clean(x):
        if isinstance(x, dict): return {k: clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [clean(v) for v in x]
        if isinstance(x, np.generic): return clean(x.item())
        if isinstance(x, float) and not math.isfinite(x): return None
        return x
    path.write_text(json.dumps(clean(data), indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def save(frame, path):
    frame.to_csv(path, sep='\t', index=False, float_format='%.12g', na_rep='NA', lineterminator='\n')


def acquire(cache, fetch):
    m = json.loads((DOCS / 'input_manifest.json').read_text())
    missing = [x for x in m['files'] if not (cache / x['path']).exists()]
    if missing:
        if not fetch: raise FileNotFoundError('Missing inputs; run with --fetch')
        raw = urlopen(Request(m['archive_url'], headers={'User-Agent': 'pd-envtox-reanalysis'}), timeout=60).read()
        outer = zipfile.ZipFile(io.BytesIO(raw))
        inner = zipfile.ZipFile(io.BytesIO(outer.read(m['inner_archive'])))
        cache.mkdir(parents=True, exist_ok=True)
        for item in missing:
            content = inner.read(item['path'])
            if hashlib.sha256(content).hexdigest() != item['sha256']:
                raise ValueError('Downloaded source differs: ' + item['path'])
            (cache / item['path']).write_bytes(content)
    for item in m['files']:
        raw = (cache / item['path']).read_bytes()
        if len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']:
            raise ValueError('Source mismatch: ' + item['path'])


def adjust(p, method='holm'):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind='stable'); n = len(p)
    if method == 'holm': values = np.maximum.accumulate(p[order] * np.arange(n, 0, -1))
    else: values = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(values, 1)
    return out


def difference(values, group=GROUP, confidence=.95, permutation=False):
    values = np.asarray(values, dtype=float); valid = np.isfinite(values)
    values = values[valid]; labels = np.asarray(group)[valid]
    c = values[~labels]; k = values[labels]
    if min(len(c), len(k)) < 2: raise ValueError('Insufficient samples for contrast')
    delta = float(k.mean() - c.mean())
    vc, vk = c.var(ddof=1) / len(c), k.var(ddof=1) / len(k)
    se = float(np.sqrt(vc + vk))
    df = float((vc + vk) ** 2 / (vc ** 2 / (len(c) - 1) + vk ** 2 / (len(k) - 1)))
    critical = stats.t.ppf((1 + confidence) / 2, df)
    record = dict(n_control=len(c), n_knockdown=len(k), log2_difference=delta,
                  ratio=2 ** delta, standard_error=se, df=df,
                  ci_low=delta - critical * se, ci_high=delta + critical * se,
                  p_two_sided=float(2 * stats.t.sf(abs(delta / se), df)))
    if permutation:
        null = []
        for selected in itertools.combinations(range(len(values)), len(k)):
            mask = np.zeros(len(values), dtype=bool); mask[list(selected)] = True
            null.append(float(values[mask].mean() - values[~mask].mean()))
        record['p_exact_decrease'] = float(np.mean(np.asarray(null) <= delta + 1e-12))
        record['assignments'] = len(null)
    return record


def prepare(path, output, return_features=False):
    d = pd.read_csv(path, dtype=str)
    original_rows = len(d)
    def single(s):
        if not isinstance(s, str): return None
        labels = {x for x in s.split(';') if x}
        return next(iter(labels)) if len(labels) == 1 else None
    d['gene'] = d.Genes.map(single)
    d = d[d.gene.notna()].copy()
    quantity_cols = [f'{run}.raw.EG.TotalQuantity (Settings)' for run in RUNS]
    qcols = [f'{run}.raw.EG.Qvalue' for run in RUNS]
    keys = ['gene', 'EG.PrecursorId']
    repeated = d[d.duplicated(keys, keep=False)]
    for _, frame in repeated.groupby(keys):
        if len(frame[qcols + quantity_cols].drop_duplicates()) > 1:
            raise ValueError('Conflicting duplicate gene/precursor rows')
    d = d.drop_duplicates(keys)
    q = d[qcols].apply(pd.to_numeric, errors='coerce').to_numpy()
    quant = d[quantity_cols].apply(pd.to_numeric, errors='coerce').to_numpy()
    valid = np.isfinite(quant) & (quant > 0) & np.isfinite(q) & (q <= .01)
    qc = pd.DataFrame({'run': RUNS, 'condition': ['NTC'] * 5 + ['NAA25_KD'] * 5,
                       'valid_precursors_before_filter': valid.sum(axis=0)})
    eligible = valid.sum(axis=1) >= 8
    d = d.loc[eligible].copy(); quant = quant[eligible]; valid = valid[eligible]
    log = np.full_like(quant, np.nan)
    log[valid] = np.log2(quant[valid])
    centered = log - np.nanmedian(log, axis=1, keepdims=True)
    features = pd.DataFrame(centered, columns=RUNS)
    features['gene'] = d.gene.to_numpy()
    features['precursor'] = d['EG.PrecursorId'].to_numpy()
    features['peptide_backbone'] = features.precursor.map(
        lambda x: re.sub(r'\[[^\]]+\]', '', x.rsplit('.', 1)[0]).replace('_', ''))
    gene = features.groupby('gene')[RUNS].median().sort_index()
    gene = gene.sub(gene.median(axis=1), axis=0)
    coverage = features.groupby('gene').agg(precursor_features=('precursor', 'nunique'),
                                           peptide_backbones=('peptide_backbone', 'nunique'))
    coverage['runs_observed'] = gene.notna().sum(axis=1)
    coverage['control_runs'] = gene.iloc[:, :5].notna().sum(axis=1)
    coverage['knockdown_runs'] = gene.iloc[:, 5:].notna().sum(axis=1)
    stable = coverage.index[(coverage.peptide_backbones >= 2) & (coverage.runs_observed == 10)]
    stable = stable.difference(['SNCA', 'NAA25'])
    if len(stable) < 100: raise ValueError('Too few loading-reference proteins')
    offset = gene.loc[stable].median(axis=0)
    normalized = gene.sub(offset, axis=1)
    qc['loading_offset_log2'] = offset.to_numpy()
    qc['valid_precursors_after_filter'] = valid.sum(axis=0)
    save(qc, output / 'sample_qc.tsv')
    save(coverage.rename_axis('gene').reset_index(), output / 'gene_coverage.tsv')
    save(normalized.rename_axis('gene').reset_index(), output / 'protein_sample_scores.tsv')
    save(features[features.gene.isin(['SNCA', 'NAA25', 'PSMG1', 'SNCG'])], output / 'target_precursor_scores.tsv')
    dump(output / 'preprocessing.json', dict(source_rows=original_rows,
        retained_precursors=len(features), quantified_genes=len(gene),
        loading_reference_genes=len(stable), missing_values_imputed=False))
    result = (gene, normalized, coverage, offset)
    return result + (features,) if return_features else result


def modules(raw, norm, sets):
    membership = {}; scores = {}; raw_scores = {}
    for name, minimum in [('proteasome_assembly', 20), ('neuronal_identity', 6)]:
        measured = [g for g in sets[name] if g in norm.index and norm.loc[g].notna().all()]
        if len(measured) < minimum: raise ValueError(f'Insufficient complete genes: {name} ({len(measured)})')
        membership[name] = measured
        scores[name] = norm.loc[measured].median(axis=0).to_numpy()
        raw_scores[name] = raw.loc[measured].median(axis=0).to_numpy()
    return membership, scores, raw_scores


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache', type=Path, default=ROOT / 'data_raw/natb_selectivity')
    p.add_argument('--output', type=Path, default=ROOT / 'results/v2/natb_selectivity')
    p.add_argument('--fetch', action='store_true')
    args = p.parse_args(); o = args.output; o.mkdir(parents=True, exist_ok=True)
    acquire(args.cache, args.fetch)
    sets = json.loads((DOCS / 'gene_sets.json').read_text())
    raw, norm, coverage, offset = prepare(args.cache / 'adj4767_Table_S5.csv', o)
    membership, score, raw_score = modules(raw, norm, sets)
    dump(o / 'module_membership.json', membership)
    snca = norm.loc['SNCA'].to_numpy(); naa = norm.loc['NAA25'].to_numpy()
    vectors = {'NAA25_manipulation': naa,
               'SNCA_relative_to_proteasome': snca - score['proteasome_assembly'],
               'SNCA_relative_to_neuronal_identity': snca - score['neuronal_identity']}
    records = []
    for name, values in vectors.items():
        if min(np.isfinite(values[:5]).sum(), np.isfinite(values[5:]).sum()) < 4:
            raise ValueError('Primary contrast has fewer than four samples per group: ' + name)
        records.append(dict(contrast=name, **difference(values, permutation=True)))
    primary = pd.DataFrame(records); primary['p_holm'] = adjust(primary.p_exact_decrease)
    save(primary, o / 'primary_contrasts.tsv')
    eq = []
    for name, values in score.items():
        result = difference(values, confidence=.90)
        result['equivalent_within_margin'] = result['ci_low'] > -MARGIN and result['ci_high'] < MARGIN
        eq.append(dict(module=name, genes=len(membership[name]), margin_log2=MARGIN, **result))
    equivalence = pd.DataFrame(eq); save(equivalence, o / 'module_equivalence.tsv')
    sample = pd.DataFrame({'run': RUNS, 'condition': ['NTC'] * 5 + ['NAA25_KD'] * 5,
                           'SNCA': snca, 'NAA25': naa, **score, **vectors})
    save(sample, o / 'primary_sample_scores.tsv')
    all_effects = []
    for gene in norm.index:
        values = norm.loc[gene].to_numpy()
        if min(np.isfinite(values[:5]).sum(), np.isfinite(values[5:]).sum()) < 4: continue
        all_effects.append(dict(gene=gene, **difference(values)))
    proteins = pd.DataFrame(all_effects); proteins['FDR_BH'] = adjust(proteins.p_two_sided, 'BH')
    proteins = proteins.merge(coverage.reset_index(), on='gene')
    save(proteins, o / 'proteome_effects.tsv')
    # Supporting specificity comparator was fixed in the protocol.
    save(proteins[proteins.gene.isin(sets['shared8'] + ['SNCA', 'NAA25', 'SNCG'])], o / 'candidate_effects.tsv')
    comparisons = []
    raw_snca = raw.loc['SNCA'].to_numpy()
    for name in ['proteasome_assembly', 'neuronal_identity']:
        a = snca - score[name]; b = raw_snca - raw_score[name]
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12, equal_nan=True)
        comparisons.append(dict(sensitivity='no_extra_loading_normalization', contrast='SNCA_relative_to_' + name,
                                **difference(b)))
        mean_score = norm.loc[membership[name]].mean(axis=0).to_numpy()
        comparisons.append(dict(sensitivity='mean_module_score', contrast='SNCA_relative_to_' + name,
                                **difference(snca - mean_score)))
    save(pd.DataFrame(comparisons), o / 'normalization_sensitivities.tsv')
    loo = []
    for i, run in enumerate(RUNS):
        keep = np.arange(10) != i
        for name, values in vectors.items():
            loo.append(dict(excluded_run=run, contrast=name, **difference(values[keep], GROUP[keep])))
    save(pd.DataFrame(loo), o / 'leave_one_run_out.tsv')
    s1 = pd.read_csv(args.cache / 'adj4767_Table_S1.csv')
    s1 = s1[s1.gene_name.isin(sets['shared8'])].copy()
    save(s1, o / 'shared8_published_melanoma_screen.tsv')
    screen_coverage = {}
    for number, col in [(1, 'gene_name'), (2, 'gene_name'), (3, 'gene')]:
        frame = pd.read_csv(args.cache / f'adj4767_Table_S{number}.csv')
        screen_coverage[f'S{number}'] = {key: sorted(set(frame[col]) & set(sets[key]))
                                        for key in ['shared8', 'proteasome_assembly']}
    dump(o / 'screen_coverage.json', screen_coverage)
    human = pd.read_csv(ROOT / 'results/v2/human_sn_proteasome/gene_level_results.tsv', sep='\t')
    human = human[human.gene.isin(sets['shared8'])]
    save(human, o / 'shared8_existing_human_results.tsv')
    supported = bool(primary.p_holm.le(.05).all() and
                     primary.log2_difference.iloc[1:].le(-MARGIN).all() and
                     equivalence.equivalent_within_margin.all())
    summary = dict(classification='SELECTIVITY_SUPPORTED_SINGLE_PEPTIDE' if supported else 'SELECTIVITY_NOT_ESTABLISHED',
                   freeze_commit='a19b720', primary_holm_all_pass=bool(primary.p_holm.le(.05).all()),
                   both_selectivity_effects_pass=bool(primary.log2_difference.iloc[1:].le(-MARGIN).all()),
                   both_modules_equivalent=bool(equivalence.equivalent_within_margin.all()),
                   SNCA_peptide_backbones=int(coverage.loc['SNCA', 'peptide_backbones']),
                   two_peptide_SNCA_sensitivity='NOT_EVALUABLE',
                   sample_unit='culture sample; not independent donor',
                   generalized_neuronal_survival_tested=False,
                   independent_dependency_experiment_reanalyzed=False)
    dump(o / 'validation_summary.json', summary)
    dump(o / 'mathematical_checks.json', {'relative_contrast_loading_invariance': 'PASS',
         'permutation_group_counts': [int(x) for x in primary.assignments],
         'sample_columns_matched_by_name': True, 'missing_values_imputed': False})
    dump(o / 'environment.json', dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__))
    print(primary.to_string(index=False)); print(equivalence.to_string(index=False)); print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
