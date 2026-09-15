#!/usr/bin/env python3
"""Descriptive diagnostics and static figures; never refit or relabel outcomes."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

from functional_rescue_analysis import ROOT, write_table

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'svg.hashsalt':'functional-rescue-v1'})
BLUE, RED, GRAY = '#22748D', '#BD593D', '#66717A'


def save(fig, out, name):
    fig.savefig(out/(name+'.png'),dpi=300,facecolor='white')
    fig.savefig(out/(name+'.pdf'),facecolor='white',
                metadata={'CreationDate':None,'ModDate':None,'Creator':'pd-envtox functional rescue'})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=ROOT/'results/v2/functional_rescue')
    args = parser.parse_args()
    out = args.results
    matched = pd.read_csv(out/'matched_benchmark.tsv',sep='\t')
    summary = json.loads((out/'analysis_summary.json').read_text())
    controls, plates, per_well = [], [], {}
    genotype_names = {'Mut - DMSO':'Mutant vehicle','WT - DMSO':'Corrected vehicle',
                      'Mut - PRO':'Prostratin','Mut;DMSO':'Mutant vehicle',
                      'WT;DMSO':'Corrected vehicle','Mut;PRO':'Prostratin'}
    for dataset in ['control','screen','followup']:
        frame = pd.read_csv(out/(dataset+'_wells.tsv.gz'),sep='\t')
        per_well[dataset] = frame
        # Use deposited control tags, never absence of a catalogue as a proxy.
        frame = frame[frame.tags.isin(genotype_names)].copy()
        frame['control_group'] = frame.tags.map(genotype_names)
        for group, g in frame.groupby('control_group',sort=True):
            row = {'dataset':dataset,'control_group':group,'n_wells':len(g)}
            for col in ['alpha_score','multivariate_score','ratio4','z1','z2','z3','z5']:
                row['mean_'+col] = g[col].mean()
            for col in ['alpha_response','integrity_0.8','viability_0.8','joint_0.8']:
                row['passing_'+col] = int(g[col].sum())
                row['fraction_'+col] = g[col].mean()
            controls.append(row)
        for (plate, group),g in frame.groupby(['Plate','control_group'],sort=True):
            plates.append({'dataset':dataset,'Plate':plate,'control_group':group,'n_wells':len(g),
                           'mean_alpha_score':g.alpha_score.mean(),
                           'mean_multivariate_score':g.multivariate_score.mean()})
    control = pd.DataFrame(controls)
    plate = pd.DataFrame(plates)
    write_table(control,out/'control_transfer_diagnostics.tsv')
    write_table(plate,out/'control_plate_scores.tsv')

    # Component pass rates are diagnostic. They do not replace joint labels.
    component_rows = []
    for _, row in matched.iterrows():
        g = per_well['followup']
        g = g[(g.SC_cat==row.SC_cat)&np.isclose(g.dose_M,row.followup_dose_M,rtol=1e-8,atol=0)&g.valid]
        component_rows.append({'SC_cat':row.SC_cat,'compound_name':row.compound_name,'n_wells':len(g),
                               'alpha':g.alpha_response.mean(),'MAP2_neurites':(g.z1>=-.5).mean(),
                               'TH_neurites':(g.z2>=-.5).mean(),'TH_intensity':(g.z3>=-.5).mean(),
                               'living_fraction':(g.z5>=-.5).mean(),'living_count':(g.ratio4>=.8).mean(),
                               'joint':g['joint_0.8'].mean(),
                               'marker_and_integrity':(g.alpha_response&g['integrity_0.8']).mean()})
    comp = pd.DataFrame(component_rows)
    write_table(comp,out/'matched_component_pass_rates.tsv')
    responders = matched[matched.followup_alpha_response].copy()
    cols = ['SC_cat','compound_name','followup_dose_M','screen_alpha_score','followup_alpha_score',
            'followup_Nuclei_Number_Living_ratio','followup_integrity_fraction_0.8',
            'followup_viability_fraction_0.8','followup_joint_fraction_0.8']
    write_table(responders[cols],out/'marker_responder_diagnostics.tsv')

    fig, axes = plt.subplots(2,2,figsize=(11.5,8.5))
    fig.subplots_adjust(left=.08,right=.97,top=.87,bottom=.20,wspace=.35,hspace=.66)
    fig.suptitle('Functional-rescue reanalysis: the combined benchmark is unevaluable',
                 x=.08,ha='left',y=.97,fontsize=15,fontweight='bold')
    fig.text(.08,.925,'31 matched compounds; published experiments; fixed analysis with a documented source-scale amendment',fontsize=10,color=GRAY)

    ax=axes[0,0]
    ax.set_title('A  Corrected-control scores shift across experiments',loc='left')
    corrected=control[control.control_group=='Corrected vehicle'].set_index('dataset')
    for offset,field,label,color in [(-.12,'mean_alpha_score','α-synuclein',BLUE),(.12,'mean_multivariate_score','Multivariate',RED)]:
        for x,dataset in enumerate(['control','screen','followup']):
            p=plate[(plate.dataset==dataset)&(plate.control_group=='Corrected vehicle')][field].to_numpy()
            jitter=np.linspace(-.045,.045,len(p))
            ax.scatter(x+offset+jitter,p,s=12,alpha=.45,color=color)
            mean=corrected.loc[dataset,field]
            ax.plot([x+offset-.07,x+offset+.07],[mean,mean],lw=3,color=color,label=label if x==0 else None)
    ax.axhline(.5,color=GRAY,lw=1,ls='--')
    ax.set_xticks([0,1,2],['Training','Screen','Follow-up'])
    ax.set_ylabel('Fixed-model score (training corrected mean = 1)')
    ax.legend(frameon=False,fontsize=9,loc='upper right')
    ax.text(.01,-.35,'Dots: plate means. Lines: pooled well means.\nNo biological error bars.',transform=ax.transAxes,fontsize=8,color=GRAY)

    ax=axes[0,1]
    ax.set_title('B  Marker ranking is repeatable within this panel',loc='left')
    ax.scatter(matched.screen_alpha_score,matched.followup_alpha_score,c=BLUE,s=28,alpha=.8)
    for _,r in responders.iterrows():
        ax.scatter(r.screen_alpha_score,r.followup_alpha_score,c=RED,s=35,zorder=3)
    lo=min(matched.screen_alpha_score.min(),matched.followup_alpha_score.min())-.15
    hi=max(matched.screen_alpha_score.max(),matched.followup_alpha_score.max())+.15
    ax.plot([lo,hi],[lo,hi],color=GRAY,lw=1,ls=':')
    ax.axhline(.5,color=GRAY,lw=1,ls='--')
    ax.set_xlabel('Screen α-synuclein score')
    ax.set_ylabel('Follow-up α-synuclein score')
    ax.text(.04,.96,f"Spearman ρ = {summary['primary']['alpha_repeatability_spearman']:.3f}\nn = 31 selected compounds",va='top',transform=ax.transAxes,fontsize=10)
    ax.text(.01,-.35,'Orange: three marker responders.\nRank correlation does not establish rescue.',transform=ax.transAxes,fontsize=8,color=GRAY)

    ax=axes[1,0]
    ax.set_title('C  Marker responders have lower living-cell counts',loc='left')
    for i,(_,r) in enumerate(responders.iterrows()):
        g=per_well['followup']
        g=g[(g.SC_cat==r.SC_cat)&np.isclose(g.dose_M,r.followup_dose_M,rtol=1e-8,atol=0)]
        median=r.followup_Nuclei_Number_Living_ratio
        ax.bar(i,median,width=.55,color=RED,alpha=.75)
        ax.scatter(i+np.linspace(-.12,.12,len(g)),g.ratio4,s=21,facecolor='white',edgecolor=RED,zorder=4)
        ax.text(i,.04,f'{median:.1%}',ha='center',color='white',fontweight='bold')
    ax.axhline(1,color=GRAY,lw=1,ls=':')
    ax.axhline(.8,color=GRAY,lw=1,ls='--')
    ax.set_xticks(range(len(responders)),responders.compound_name,fontsize=9)
    ax.set_ylim(0,1.12)
    ax.set_ylabel('Living nuclei / plate mutant-vehicle median')
    ax.text(.01,-.35,'Bars: median of four wells at 5.13 μM.\nDots: individual wells.',transform=ax.transAxes,fontsize=8,color=GRAY)

    ax=axes[1,1]
    ax.set_title('D  No positives under the fixed combined rule',loc='left')
    values=[len(matched),int(matched.followup_alpha_response.sum()),int((comp.marker_and_integrity>=.75).sum()),int(matched['followup_joint_response_0.8'].sum())]
    labels=['Matched','Marker','Marker +\nneuronal features','All criteria']
    ax.bar(range(4),values,color=[GRAY,BLUE,RED,RED],width=.6)
    for i,v in enumerate(values):ax.text(i,v+.65,str(v),ha='center',fontweight='bold')
    ax.set_xticks(range(4),labels,fontsize=9)
    ax.set_ylim(0,36)
    ax.set_ylabel('Compounds (≥75% passing wells)')
    ax.text(.01,-.35,'Corrected follow-up controls also have\n0/48 jointly passing wells.',transform=ax.transAxes,fontsize=8,color=RED)
    fig.text(.08,.03,'Calibration and endpoint limitations prevent a rescue verdict.\nThis reanalysis does not establish a molecular mechanism.',fontsize=10,color=GRAY,linespacing=1.5)
    save(fig,out,'fig1_benchmark_and_calibration')

    order=matched.sort_values(['followup_alpha_score','SC_cat'],ascending=[False,True]).SC_cat
    c=comp.set_index('SC_cat').loc[order]
    fields=['alpha','MAP2_neurites','TH_neurites','TH_intensity','living_fraction','living_count','joint']
    fig,ax=plt.subplots(figsize=(10,12))
    fig.subplots_adjust(left=.31,right=.9,top=.88,bottom=.13)
    cmap=LinearSegmentedColormap.from_list('passfraction',['#F1DAD2','#FAF8F3','#27788B'])
    im=ax.imshow(c[fields].to_numpy(),aspect='auto',vmin=0,vmax=1,cmap=cmap)
    ax.set_yticks(range(len(c)),c.compound_name,fontsize=9)
    ax.set_xticks(range(len(fields)),['α-synuclein','MAP2\nneurites','TH\nneurites','TH\nintensity','Living\nfraction','Living\ncount','Joint'],fontsize=9)
    ax.tick_params(axis='both',length=0)
    for i in range(len(c)):
        for j,v in enumerate(c[fields].iloc[i]):
            ax.text(j,i,f'{v:.0%}',ha='center',va='center',fontsize=8,color='white' if v>=.75 else '#31383D')
    ax.set_xticks(np.arange(-.5,len(fields),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(c),1),minor=True)
    ax.grid(which='minor',color='white',lw=1)
    ax.tick_params(which='minor',length=0)
    cb=fig.colorbar(im,ax=ax,location='right',fraction=.03,pad=.04)
    cb.set_ticks([0,.25,.5,.75,1],labels=['0%','25%','50%','75%','100%'])
    cb.set_label('Fraction of four wells meeting each criterion')
    fig.suptitle('Follow-up response components for all 31 matched compounds',x=.06,ha='left',y=.965,fontsize=15,fontweight='bold')
    fig.text(.06,.925,'Ordered by follow-up marker score. Component values are diagnostic and do not change the fixed endpoint.',fontsize=10,color=GRAY)
    fig.text(.06,.058,'Marker score ≥0.5; normalized neuronal features and living fraction ≥−0.5 control SD; living count ≥80% of vehicle.\nA compound needs ≥75% jointly passing wells. These are operational cutoffs, not validated toxicity or safety thresholds.\nThe corrected controls also fail the joint rule, so a red cell cannot by itself be interpreted as biological harm.',fontsize=9,color=GRAY,linespacing=1.5)
    save(fig,out,'fig2_all_matched_response_components')
    print('Wrote diagnostic tables and two PNG/PDF figures from unchanged result tables.')


if __name__=='__main__':
    main()
