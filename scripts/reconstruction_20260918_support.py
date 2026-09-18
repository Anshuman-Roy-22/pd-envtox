#!/usr/bin/env python3
"""Descriptive balance, litter sensitivity, discrepancy checks and scientific figure."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from reconstruct_mouse_20260918 import ROOT, DOC, OUT, FOCAL, save, dump
from reconstruct_mouse_fit_20260918 import test

def main():
    e = pd.read_csv(OUT/'all_mouse_effects.tsv', sep='\t')
    locked = pd.read_csv(DOC/'locked_comparators.tsv', sep='\t')
    pairs = pd.read_csv(OUT/'specificity_mouse_pairs.tsv', sep='\t')
    # Litter identities already cross-referenced to the author Animal table by the fit.
    litter = pairs[['mouse','litter']].drop_duplicates().set_index('mouse').litter
    # All focal mice, including the excluded specificity mice, require the full author map.
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);a=p.parse_args()
    animals=pd.read_excel(a.cache/'media-3.xlsx',header=1)[['Animal','Litter #']].dropna()
    litter=dict(zip('mouse'+animals.Animal.astype(int).astype(str),animals['Litter #']))
    rows=[];lr=[]
    for (target,readout),q in e[e.primary&e.target.isin(FOCAL)&e.readout.isin(['synapse','proteasome'])].groupby(['target','readout']):
        v=q.assign(litter=q.mouse.map(litter)).groupby('litter').difference.mean()
        assert q.mouse.map(litter).notna().all()
        rows.append(dict(target=target,readout=readout,**test(v)))
        lr.extend(dict(target=target,readout=readout,litter=k,effect=x) for k,x in v.items())
    t=pd.DataFrame(rows);t['BH_eight_two']=stats.false_discovery_control(t.p_two.to_numpy())
    save(t,OUT/'focal_litter_tests.tsv');save(pd.DataFrame(lr),OUT/'focal_litter_effects.tsv')
    balance=[];refrows=[];challenge=[]
    for target in FOCAL:
        q=locked[locked.focal==target]
        row=dict(target=target,comparators=len(q))
        for col in ['dependency_delta','coverage_log2_ratio','expression_log2_ratio']:
            row[col+'_min']=q[col].min();row[col+'_max']=q[col].max();row[col+'_mean']=q[col].mean()
        balance.append(row)
        mice=set(pairs[pairs.target==target].mouse)
        for gene in q.mouse:
            vals=e[e.primary&(e.target==gene)&(e.readout=='synapse')&e.mouse.isin(mice)].difference
            refrows.append(dict(focal=target,comparator=gene,overlapping_mice=len(vals),mean=float(vals.mean()) if len(vals) else np.nan))
        vals=pd.read_csv(OUT/f'{target}_essential_challenge.tsv',sep='\t').paired_difference
        challenge.append(dict(target=target,interpretation='Descriptive, dependency-bottom-quintile shortlist; not expression-matched',**test(vals)))
    save(pd.DataFrame(balance),OUT/'comparator_balance.tsv')
    ref=pd.DataFrame(refrows);ref['within_panel_descriptive_rank']=ref.groupby('focal')['mean'].rank(method='min');save(ref,OUT/'comparator_reference_means.tsv')
    save(pd.DataFrame(challenge),OUT/'essential_challenge_summary.tsv')
    # These known values are discrepancy checks, never optimizer targets.
    focal=pd.read_csv(OUT/'focal_tests.tsv',sep='\t');f=focal[(focal.analysis=='primary')&(focal.readout=='synapse')].set_index('target')
    expected={'Naa20':-.0130293359097,'Pomp':-.0325933585595,'Psmb4':-.0770265357707}
    delta={g:float(f.loc[g,'mean']-x) for g,x in expected.items()}
    assert max(map(abs,delta.values()))<1e-10
    cal=json.loads((OUT/'calibration_summary.json').read_text());assert cal['evaluable_groups']==99 and cal['nominal_down_positives']==3
    h=ROOT/'results/v2/reconstruction_20260918/human';proteins=pd.read_csv(h/'limma_primary.tsv',sep='\t')
    assert len(proteins)==4511 and (proteins['adj.P.Val']<.05).sum()==994
    relative=pd.read_csv(h/'relative_tests.tsv',sep='\t')
    assert np.allclose(relative.BH_four,4/252,atol=1e-10)
    cam=pd.read_csv(h/'focus_camera.tsv',sep='\t');q=cam[(cam.pathway=='R-HSA-112315')&(cam.model=='primary')&(cam.setting=='estimated')].FDR.iloc[0]
    assert abs(q-.626542587395541)<1e-10
    dump(OUT/'reconstruction_checks.json',{'focal_difference_from_known_effect':delta,'calibration_reproduced':True,'human_994_of_4511_reproduced':True,'human_relative_four_tests_reproduced':True,'estimated_CAMERA_FDR':float(q),'scope':'Numeric discrepancy checks against previously reported results; not independent confirmation'})
    figure(h, focal)

def figure(h,focal):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(12,5.8),gridspec_kw={'width_ratios':[1,1.2]})
    r=pd.read_csv(h/'relative_tests.tsv',sep='\t')
    labels=['Chemical synapse\nvs all background','Chemical synapse\nvs matched background','Neurotransmitter release\nvs all background','Neurotransmitter release\nvs matched background']
    for y,(_,row) in zip(range(4),r.iterrows()):
        ax[0].errorbar(row.log2_difference,y,xerr=[[row.log2_difference-row.ci_low],[row.ci_high-row.log2_difference]],fmt='o',color='#196d82',capsize=3)
    ax[0].set(yticks=range(4),yticklabels=labels,ylim=(3.6,-.6),xlabel='Relative protein score difference (log2 units)',title='A  Human NAA25 knockdown')
    ax[0].text(.02,-.23,'5 cultures/group; exact permutation BH4 = 0.0159\nEstimated-correlation CAMERA FDR = 0.627',transform=ax[0].transAxes,fontsize=9,va='top')
    m=pd.read_csv(OUT/'specificity_tests.tsv',sep='\t')
    for y,target in enumerate(FOCAL):
        z=focal[(focal.analysis=='primary')&(focal.target==target)&(focal.readout=='synapse')].iloc[0]
        ax[1].errorbar(z['mean'],y-.12,xerr=[[z['mean']-z.ci_low],[z.ci_high-z['mean']]],fmt='o',color='#9aa4ac',capsize=3,label='Versus non-targeting' if y==0 else None)
        z=m[m.target==target].iloc[0]
        if z.status=='EVALUABLE':
            ax[1].errorbar(z['mean'],y+.12,xerr=[[z['mean']-z.ci_low],[z.ci_high-z['mean']]],fmt='s',color='#b45325',capsize=3,label='Versus matched perturbations' if y==0 else None)
        else:ax[1].text(-.006,y+.12,'Unevaluable',ha='right',va='center',fontsize=8,color='#b45325')
    ax[1].set(yticks=range(4),yticklabels=FOCAL,ylim=(3.6,-.6),xlabel='Synaptic RNA score difference (mean log1p units)',title='B  Mouse perturbation comparison')
    ax[1].legend(loc='lower left',bbox_to_anchor=(-.08,-.32),frameon=False,fontsize=9)
    for a in ax:a.axvline(0,color='#333333',linewidth=.8,linestyle='--');a.grid(axis='x',alpha=.15)
    fig.suptitle('Relative synaptic responses with explicit specificity limits',fontsize=15,y=.99)
    fig.subplots_adjust(left=.18,right=.98,top=.86,bottom=.27,wspace=.65)
    fig.savefig(OUT.parent/'specificity_comparison.png',dpi=300)
    fig.savefig(OUT.parent/'specificity_comparison.pdf',metadata={'CreationDate':None,'ModDate':None})
    plt.close(fig)

if __name__=='__main__':main()
