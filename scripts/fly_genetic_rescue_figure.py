#!/usr/bin/env python3
"""Draw the fixed benchmark from its generated result tables."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, default=ROOT / 'results/v2/fly_genetic_rescue')
    args = p.parse_args()
    d = args.results
    effects = pd.read_csv(d / 'genotype_effects.tsv', sep='\t')
    sensitivities = pd.read_csv(d / 'sensitivity_summary.tsv', sep='\t')
    summary = json.loads((d / 'validation_summary.json').read_text())
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'svg.hashsalt': 'fly-genetic-rescue-v1'})
    fig = plt.figure(figsize=(13.3, 8.3), facecolor='white')
    grid = fig.add_gridspec(2, 2, left=.09, right=.97, bottom=.16, top=.82,
                           wspace=.56, hspace=.62, width_ratios=[1.03, 1])
    ax = fig.add_subplot(grid[:, 0])
    y = np.arange(len(effects))
    ax.axvline(0, color='#ADB6BF', lw=1, zorder=0)
    ax.scatter(effects.Q10_adjusted_log2_response, y + .12, color='#246A9A', s=33, label='Q10')
    ax.scatter(effects.R55_adjusted_log2_response, y - .12, color='#D17923', s=33, label='R55')
    ax.set_yticks(y, effects.locus)
    ax.invert_yaxis()
    ax.set_xlim(-.54, 1.02)
    ax.set_xlabel('Control-adjusted response (log2 ratio)')
    ax.set_title('A. Innervation responses', loc='left', fontweight='bold', pad=14)
    ax.legend(loc='lower right', bbox_to_anchor=(1, 1.005), ncol=2,
              borderaxespad=0, frameon=False, fontsize=9)
    ax.grid(axis='x', color='#E8ECEF', lw=.7)
    ax.set_axisbelow(True)
    ax2 = fig.add_subplot(grid[0, 1])
    labels = ['Average response', 'Baseline severity', 'Genetic interactions']
    vals = [summary['MAE_mean'], summary['MAE_severity'], summary['MAE_GI']]
    colors = ['#A5ADB7', '#798590', '#943F4E']
    ax2.barh(np.arange(3), vals, color=colors, height=.52)
    for i, value in enumerate(vals):
        ax2.text(value + .012, i, f'{value:.3f}', va='center', fontsize=10)
    ax2.set_yticks(np.arange(3), labels)
    ax2.invert_yaxis()
    ax2.set_xlim(0, .64)
    ax2.set_xlabel('Mean absolute prediction error; lower is better')
    ax2.set_title('B. Primary genotype-withheld prediction', loc='left', fontweight='bold', pad=14)
    ax2.text(0, -.34, f"NOT_SUPPORTED  |  randomization p = {summary['randomization_p']:.3f}",
             transform=ax2.transAxes, fontsize=10, color='#943F4E')
    ax3 = fig.add_subplot(grid[1, 1])
    names = ['Primary: 3 neighbors', '1 neighbor', '5 neighbors', 'Median arm values', 'Exclude younger ages']
    gains = sensitivities.improvement_over_mean.to_numpy() * 100
    ax3.barh(np.arange(5), gains, color='#943F4E', height=.52)
    for i, value in enumerate(gains):
        ax3.text(value - .6, i, f'{value:.1f}%', ha='right', va='center', fontsize=9)
    ax3.axvline(0, color='#505A64', lw=1)
    ax3.set_yticks(np.arange(5), names)
    ax3.invert_yaxis()
    ax3.set_xlim(-36, 12)
    ax3.set_xlabel('Error improvement over average-response baseline (%)')
    ax3.set_title('C. Fixed sensitivity analyses', loc='left', fontweight='bold', pad=14)
    fig.text(.06, .95, 'Interaction-based model did not improve response prediction',
             fontsize=19, fontweight='bold', color='#253A4D')
    fig.text(.06, .908, 'Retrospective reanalysis of published fly experiments | 18 genotypes | one study',
             fontsize=11, color='#52606D')
    fig.text(.06, .077, 'Responses are pooled source means of normalized TH-positive innervation area. Missing vial/batch identifiers prevent animal-level inference.',
             fontsize=9, color='#52606D')
    fig.text(.06, .051, 'Drug preference is not proof of benefit. Prediction failure does not invalidate the authors\' intervention results. Source: Kaempf et al., 2026.',
             fontsize=9, color='#52606D')
    fig.savefig(d / 'genetic_interaction_prediction.png', dpi=300,
                metadata={'Software': 'pd-envtox frozen fly benchmark'})
    fig.savefig(d / 'genetic_interaction_prediction.pdf',
                metadata={'CreationDate': None, 'ModDate': None, 'Creator': 'pd-envtox'})
    plt.close(fig)


if __name__ == '__main__':
    main()
