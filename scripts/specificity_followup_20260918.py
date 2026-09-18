#!/usr/bin/env python3
"""Exact shift-test inversion and a descriptive dependency plot from frozen tables."""
from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/v2/specificity_followup_20260918'
OLD=ROOT/'results/v2/reconstruction_20260918'
DOC=ROOT/'docs/v2/reconstruction_20260918'
def save(d,name):d.to_csv(OUT/name,sep='\t',index=False,float_format='%.14g',na_rep='NA',lineterminator='\n')
def invert(y,treated):
    masks=np.array([[i in ix for i in range(10)] for ix in itertools.combinations(range(10),5)],dtype=float)
    weights=(2*masks-1)/5; a=weights@y; b=weights@treated
    delta=float(y[treated.astype(bool)].mean()-y[~treated.astype(bool)].mean())
    def p(tau):return np.mean(np.abs(a-b*tau)>=abs(delta-tau)-1e-12)
    critical=[]
    for sign in [-1,1]:
        denominator=b-sign; numerator=a-sign*delta;ok=abs(denominator)>1e-12
        critical.extend((numerator[ok]/denominator[ok]).tolist())
    critical=np.unique(np.array(critical))
    bound=2*max(1,float(np.max(abs(critical))))+1
    knots=np.r_[-bound,critical,bound];points=[]
    for i,x in enumerate(knots):
        points.append((x,p(x),'breakpoint' if 0<i<len(knots)-1 else 'tail_check'))
        if i+1<len(knots):
            mid=(x+knots[i+1])/2;points.append((mid,p(mid),'interior'))
    frame=pd.DataFrame(points,columns=['shift','p_exact','location']);frame['accepted_95']=frame.p_exact>.05
    accepted=np.flatnonzero(frame.accepted_95)
    assert len(accepted) and accepted[0]>0 and accepted[-1]<len(frame)-1
    assert np.all(np.diff(accepted)==1),'Disconnected confidence set; stop rather than report a hull'
    lo=float(frame.iloc[accepted[0]]['shift']);hi=float(frame.iloc[accepted[-1]]['shift'])
    assert p(lo)>.05 and p(hi)>.05 and p(delta)==1
    # Independent direct permutations at the estimate, null and both boundaries.
    for tau in [0,delta,lo,hi]:
        shifted=y-tau*treated;null=weights@shifted;observed=shifted[treated.astype(bool)].mean()-shifted[~treated.astype(bool)].mean()
        direct=np.mean(abs(null)>=abs(observed)-1e-12);assert direct==p(tau)
    return delta,lo,hi,float(p(0)),frame
def human():
    d=pd.read_csv(OLD/'human/relative_culture_scores.tsv',sep='\t');old=pd.read_csv(OLD/'human/relative_tests.tsv',sep='\t');rows=[];detail=[]
    for (pathway,comparison),q in d.groupby(['pathway','comparison'],sort=False):
        q=q.sort_values('run');z=q.run.str.startswith('NAA25').to_numpy(dtype=float);assert z.sum()==5 and len(z)==10
        delta,lo,hi,p,grid=invert(q.relative.to_numpy(),z);ref=old[(old.pathway==pathway)&(old.comparison==comparison)].iloc[0]
        assert abs(delta-ref.log2_difference)<1e-10 and abs(p-ref.p_exact)<1e-10
        rows.append(dict(pathway=pathway,comparison=comparison,effect=delta,permutation_CI_low=lo,permutation_CI_high=hi,p_exact=p,BH_four=ref.BH_four,Welch_CI_low=ref.ci_low,Welch_CI_high=ref.ci_high,assignments=252,two_sided_minimum_p=2/252))
        detail.append(grid.assign(pathway=pathway,comparison=comparison))
    result=pd.DataFrame(rows);save(result,'human_permutation_intervals.tsv');save(pd.concat(detail,ignore_index=True),'human_inversion_breakpoints.tsv');print(result.to_string(index=False),flush=True)
    return result
def dependency():
    e=pd.read_csv(OLD/'mouse/all_mouse_effects.tsv',sep='\t');m=pd.read_csv(DOC/'locked_comparators.tsv',sep='\t');cov=pd.read_csv(DOC/'all_candidate_covariates.tsv',sep='\t')
    comparators=sorted(set(m.mouse));focal=['Naa20','Pomp','Psmb4','Psmc5']
    q=e[e.primary&(e.readout=='synapse')&e.target.isin(comparators+focal)].groupby('target').difference.agg(['mean','count']).reset_index()
    q=q.merge(cov[['mouse','percentile_plusAO']],left_on='target',right_on='mouse',validate='one_to_one').drop(columns='mouse');q['role']=np.where(q.target.isin(focal),'focal','comparator');save(q,'dependency_gene_means.tsv')
    c=q[q.role=='comparator'];pr=stats.pearsonr(c.percentile_plusAO,c['mean']);sp=stats.spearmanr(c.percentile_plusAO,c['mean'])
    save(pd.DataFrame([dict(n_comparator_genes=len(c),pearson_r=pr.statistic,pearson_descriptive_p=pr.pvalue,spearman_rho=sp.statistic,spearman_descriptive_p=sp.pvalue)]),'dependency_associations.tsv')
    return q
def figures(h,d):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(12,5.7));labels=['Chemical synapse / all','Chemical synapse / matched','Release / all','Release / matched']
    for i,r in h.iterrows():
        ax[0].errorbar(r.effect,i,xerr=[[r.effect-r.permutation_CI_low],[r.permutation_CI_high-r.effect]],fmt='o',color='#196d82',capsize=3)
    ax[0].set(yticks=range(4),yticklabels=labels,ylim=(3.6,-.6),xlabel='Relative protein-score difference (log2)',title='A  Exact permutation-inverted 95% intervals');ax[0].axvline(0,color='gray',ls='--',lw=1)
    c=d[d.role=='comparator'];ax[1].scatter(c.percentile_plusAO,c['mean'],s=28,alpha=.75,color='#7d8d98',label='Matched comparator genes')
    for i,g in enumerate(['Naa20','Pomp','Psmb4','Psmc5']):
        r=d[d.target==g].iloc[0];ax[1].scatter(r.percentile_plusAO,r['mean'],marker='D' if g!='Psmb4' else 'x',s=50,color='#b45325')
        ax[1].annotate(g,(r.percentile_plusAO,r['mean']),xytext=(6,-12 if g=='Naa20' else 6),textcoords='offset points',fontsize=9)
    ax[1].axhline(0,color='gray',lw=1,ls='--');ax[1].set(xlabel='External neuronal-dependency percentile',ylabel='Mean synaptic RNA difference (log1p)',title='B  Descriptive dependency comparison');ax[1].legend(frameon=False,fontsize=9,loc='lower left')
    fig.subplots_adjust(left=.20,right=.98,bottom=.22,top=.88,wspace=.5)
    fig.text(.02,.035,'A: constant additive-shift/exchangeability model; pointwise intervals. B: one point per gene; mouse averages.\nPsmb4 is shown descriptively; its matched specificity test remains unevaluable. Global-disruption regression is pending.',fontsize=9)
    fig.savefig(OUT/'intervals_and_dependency.png',dpi=300);fig.savefig(OUT/'intervals_and_dependency.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    cam=pd.read_csv(OUT/'camera_grid_all_pathways.tsv',sep='\t');focus=['R-HSA-112315','R-HSA-112310'];estimated=pd.read_csv(OUT/'camera_estimated_all_pathways.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(8,5.4));colors=['#196d82','#b45325'];names=['Chemical synapse','Neurotransmitter release']
    for pid,color,name in zip(focus,colors,names):
        z=cam[cam.pathway==pid].sort_values('fixed_correlation');ax.plot(z.fixed_correlation,z.FDR,'o-',color=color,label=name,ms=4)
    ax.axhline(.05,color='#333333',ls='--',lw=1,label='FDR 0.05');ax.set(xlabel='Common fixed inter-gene correlation',ylabel='CAMERA FDR (862-pathway family)',yscale='log',title='Common-correlation CAMERA sensitivity');ax.legend(frameon=False,fontsize=9,loc='lower right')
    notes=[]
    for pid,name in zip(focus,names):
        z=estimated[estimated.pathway==pid].iloc[0];notes.append(f'{name}: estimated correlation {z.Correlation:.4f}, FDR {z.FDR:.3f}')
    fig.subplots_adjust(bottom=.28,top=.90,left=.12,right=.97);fig.text(.08,.05,'\n'.join(notes)+'\nEstimated mode uses pathway-specific correlations and a different denominator df.',fontsize=9)
    fig.savefig(OUT/'camera_correlation_sensitivity.png',dpi=300);fig.savefig(OUT/'camera_correlation_sensitivity.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
def main():
    OUT.mkdir(parents=True,exist_ok=True);h=human();d=dependency();figures(h,d)
if __name__=='__main__':main()
