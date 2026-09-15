#!/usr/bin/env python3
"""Render source-linked molecular selectivity and robustness figures."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TEAL, ORANGE, GRAY = '#197C80', '#B54C25', '#5C6370'


def finish(fig, path):
    fig.savefig(path.with_suffix('.png'), dpi=300, facecolor='white')
    fig.savefig(path.with_suffix('.pdf'), facecolor='white',
                metadata={'Creator': 'pd-envtox reproducible analysis', 'CreationDate': None, 'ModDate': None})
    plt.close(fig)


def decorate(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.axvline(0, color='#8A8F98', lw=1, zorder=0)
    ax.grid(axis='x', color='#E7E9ED', lw=.7, zorder=0)
    ax.set_axisbelow(True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=ROOT / 'results/v2/natb_selectivity')
    a = p.parse_args(); o = a.output
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42,
                         'axes.titleweight': 'bold', 'axes.labelcolor': '#20242A',
                         'text.color': '#20242A', 'savefig.dpi': 300})
    effects = pd.read_csv(o / 'proteome_effects.tsv', sep='\t').set_index('gene')
    modules = pd.read_csv(o / 'module_robustness.tsv', sep='\t')
    modules = modules[modules.aggregation == 'median'].set_index('contrast')
    samples = pd.read_csv(o / 'primary_sample_scores.tsv', sep='\t')

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 6.7), gridspec_kw={'width_ratios': [1.55, 1]})
    fig.subplots_adjust(left=.19, right=.965, top=.78, bottom=.24, wspace=.48)
    fig.text(.05, .945, 'NAA25 knockdown reveals a synaptic-protein selectivity concern',
             fontsize=16, weight='bold')
    fig.text(.05, .89, 'Human iNeuron proteomics · 5 control and 5 knockdown culture samples · retrospective reanalysis', fontsize=10.5)
    labels = ['NAA25', 'SNCA', 'SNCG', 'SYP', 'SYT1', 'Proteasome panel (39)', 'Neuronal panel (7)']
    for i, name in enumerate(labels):
        if i < 5:
            row = effects.loc[name]
            color = ORANGE if name in ['SYP', 'SYT1'] else TEAL
        else:
            row = modules.loc['proteasome_assembly' if i == 5 else 'neuronal_identity']
            color = GRAY if i == 5 else ORANGE
        ax.errorbar(row.log2_difference, i, xerr=[[row.log2_difference-row.ci_low],
                    [row.ci_high-row.log2_difference]], fmt='o', color=color, capsize=3, lw=1.5)
    ax.set_yticks(range(7), labels); ax.invert_yaxis(); ax.set_ylim(6.6, -.6)
    ax.set_xlim(-2.35, 1.05); ax.set_xticks([-2, -1, 0, 1])
    ax.set_xlabel('Knockdown − control (log₂ abundance)')
    ax.set_title('A  Protein and panel effects', loc='left', pad=16, fontsize=12)
    decorate(ax)
    values = samples.neuronal_identity - samples.proteasome_assembly
    jitter = np.array([-.11, -.05, 0, .05, .11])
    for start, x, color in [(0, 0, TEAL), (5, 1, ORANGE)]:
        v = values.iloc[start:start+5].to_numpy()
        bx.scatter(x+jitter, v, color=color, s=48, edgecolor='white', linewidth=.6, zorder=3)
        bx.plot([x-.19,x+.19], [v.mean(),v.mean()], color=color, lw=2)
    bx.set_xticks([0,1], ['Control', 'NAA25 KD']); bx.set_xlim(-.42,1.42)
    bx.set_ylabel('Neuronal − proteasome panel score (log₂)')
    bx.spines[['top','right']].set_visible(False)
    bx.axhline(0, color='#8A8F98', lw=1, zorder=0)
    bx.set_title('B  Within-sample panel contrast', loc='left', pad=16, fontsize=12)
    bx.text(.5, -.26, 'Exact two-sided p = 0.00794\nExploratory contrast; dots are culture samples',
            ha='center', va='top', transform=bx.transAxes, fontsize=9.5)
    fig.text(.05, .105, 'Bars: Welch 95% confidence intervals. Protein FDR: SYP 0.024; SYT1 0.034 (4,511 proteins tested).', fontsize=10)
    fig.text(.05, .068, 'SNCA: one peptide, 5 vs 4 observed samples. Panel preservation concerns relative abundance, not enzyme or synaptic activity.', fontsize=9)
    fig.text(.05, .035, 'Source: Santhosh Kumar et al., Science Advances (2024), Table S5. The frozen selectivity gate did not pass.', fontsize=9)
    finish(fig, o / 'fig1_intervention_selectivity')

    robust = pd.read_csv(o / 'backbone_target_effects.tsv', sep='\t')
    backbones = pd.read_csv(o / 'synaptic_backbone_effects.tsv', sep='\t')
    fig, (ax,bx) = plt.subplots(1,2,figsize=(12.7,7.4), gridspec_kw={'width_ratios':[1.35,1]})
    fig.subplots_adjust(left=.21,right=.96,top=.79,bottom=.23,wspace=.6)
    fig.text(.05,.945,'The SYT1 decrease survives complete-peptide reanalysis',fontsize=16,weight='bold')
    fig.text(.05,.89,'Sensitivity analyses of peptide coverage, charge-state weighting and sample influence.',fontsize=10.5)
    yt=[]
    for i,(gene,variant,label) in enumerate([
        ('SYP',None,'SYP · original'),
        ('SYP','equal_backbones_at_least_8_runs','SYP · equal backbones'),
        ('SYP','equal_backbones_complete_10_runs','SYP · complete, ≥2 backbones'),
        ('SYT1',None,'SYT1 · original'),
        ('SYT1','equal_backbones_at_least_8_runs','SYT1 · equal backbones'),
        ('SYT1','equal_backbones_complete_10_runs','SYT1 · complete, ≥2 backbones')]):
        yt.append(label)
        if variant is None: row=effects.loc[gene]
        else:
            z=robust[(robust.gene==gene)&(robust.analysis==variant)]
            if z.empty:
                ax.text(-1.95,i,'Not evaluable: only 1 complete backbone',fontsize=9,va='center',color=GRAY)
                continue
            row=z.iloc[0]
        ax.errorbar(row.log2_difference,i,xerr=[[row.log2_difference-row.ci_low],[row.ci_high-row.log2_difference]],
                    fmt='o',color=ORANGE if gene=='SYT1' else TEAL,capsize=3,lw=1.5)
    ax.set_yticks(range(6),yt); ax.set_ylim(5.7,-.7); ax.set_xlim(-2.3,.15)
    ax.set_xlabel('Knockdown − control (log₂ abundance)')
    ax.set_title('A  Protein estimates with 95% intervals',loc='left',pad=16,fontsize=11.5)
    decorate(ax)
    z=backbones[(backbones.gene=='SYT1')&(backbones.analysis=='equal_backbones_complete_10_runs')].sort_values('peptide_backbone')
    bx.scatter(z.log2_difference,np.arange(1,len(z)+1),s=40,color=ORANGE,zorder=3)
    bx.set_yticks([1,4,7,10,14]);bx.set_ylim(14.8,.2);bx.set_xlim(min(-3.1,z.log2_difference.min()-.2),.2)
    bx.set_ylabel('SYT1 peptide backbone index')
    bx.set_xlabel('Change in log₂ abundance')
    bx.set_title('B  Complete SYT1 peptides\nAll 14 point downward',loc='left',pad=16,fontsize=11.5)
    decorate(bx)
    fig.text(.05,.135,'SYT1: 63% lower in the complete-peptide analysis; FDR = 0.038 across 2,183 eligible proteins.',fontsize=10.5)
    fig.text(.05,.095,'Both proteins and the neuronal-panel effect remain negative after every single-sample omission.',fontsize=10)
    fig.text(.05,.058,'Peptides are repeated measurements of the same samples, not 14 independent experiments. SYP fails the strict coverage gate.',fontsize=9)
    fig.text(.05,.026,'These checks are post hoc within-study robustness analyses. They do not establish independent replication or impaired neurotransmission.',fontsize=9)
    finish(fig,o/'fig2_synaptic_robustness')


if __name__ == '__main__':
    main()
