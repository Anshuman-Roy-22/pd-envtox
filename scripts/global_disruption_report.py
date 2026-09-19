"""Report every model outcome without upgrading mechanistic claims."""
import json
import numpy as np
import pandas as pd
from global_disruption_metrics import ROOT,OUT
from global_disruption_io import dump

def main():
    status=json.loads((OUT/'model_status.json').read_text())
    text=['# Global-disruption sensitivity analysis','',f"Model status: **{status['status']}**.",'','Original matched-perturbation findings are preserved. This is a descriptive adjustment of measured transcriptomic disruption and realized target expression. It does not measure viability, proteasome activity, or synaptic function, and cannot exclude residual confounding.','']
    if status['status']=='NOT_EVALUABLE':
        text+=['The prespecified model could not be evaluated: '+status['reason'], 'No replacement model was fitted.']
    else:
        q=pd.read_csv(OUT/'focal_residual_summary.tsv',sep='\t')
        text+=['Focal residual = observed synaptic contrast minus the comparator-model population prediction. Negative residuals mean a more negative synaptic response than predicted from the measured covariates.','',f"Comparator genes: {status['comparator_genes']}; observations: {status['training_observations']}; successful litter refits: {status['bootstrap_successful']}/1000.",'','| Target | Scope | Residual | Descriptive 95% interval | Mice | Extrapolated mice |','|---|---|---:|---|---:|---:|']
        for r in q.itertuples():
            interval=f'[{r.ci_low:.5g}, {r.ci_high:.5g}]' if pd.notna(r.ci_low) else 'NOT_EVALUABLE'
            text.append(f'| {r.target} | {r.scope} | {r.mean_residual:.5g} | {interval} | {r.mice} | {r.extrapolated_mice} |')
        text+=['','Intervals use whole-litter resampling and are descriptive. Range overlap is a marginal check on two biological predictors, not proof of joint covariate support. Out-of-range results are extrapolations. Realized-knockdown measurement and overlap exclusions may restrict generalizability.','', 'Plots of raw gene means are descriptive; their slopes are not the multivariable adjusted model. No favorable direction or interval constitutes a new mechanistic success criterion.']
        figure(q)
    text+=['','Source: public mouse CRISPR experiments; observational secondary analysis of experimental perturbations. Source selection and this follow-up are retrospective. Psmb4 remains descriptive and is not restored to the primary matched panel.','', 'Review global_metrics.tsv, comparator_measurement_eligibility.tsv, focal_predictions.tsv, coefficients.tsv, source_reconciliation.json, and every bootstrap failure before using these results in the paper.']
    (OUT/'RESULTS.md').write_text('\n'.join(text)+'\n',encoding='utf-8')

def figure(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=pd.read_csv(OUT/'global_metrics.tsv',sep='\t')
    genes=d.groupby(['target','role']).agg(global_rms=('global_rms','mean'),synaptic_effect=('difference','mean')).reset_index()
    dep=pd.read_csv(ROOT/'results/v2/specificity_followup_20260918/dependency_gene_means.tsv',sep='\t')
    genes=genes.merge(dep[['target','percentile_plusAO']],on='target',validate='one_to_one')
    colors={'Naa20':'#167b93','Pomp':'#d17b22','Psmc5':'#754a99','Psmb4':'#69747b'}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,3,figsize=(14,5.5),layout='constrained')
    for ax,x,label in [(axes[0],'global_rms','Global transcriptomic RMS'),(axes[1],'percentile_plusAO','Independent dependency percentile')]:
        c=genes[genes.role=='comparator'];ax.scatter(c[x],c.synaptic_effect,color='#a7afb6',s=23,label='Fixed comparators')
        for target,col in colors.items():
            p=genes[genes.target==target]
            ax.scatter(p[x],p.synaptic_effect,color=col,s=42,label=target)
            if len(p):ax.annotate(target,(p[x].iloc[0],p.synaptic_effect.iloc[0]),xytext=(4,4),textcoords='offset points',fontsize=8)
        ax.axhline(0,color='black',lw=.6);ax.set(xlabel=label,ylabel='Mean synaptic effect (log1p units)')
    axes[0].set_title('A  Equal-mouse gene means')
    axes[1].set_title('B  Dependency comparison')
    ax=axes[2]
    for i,t in enumerate(['Naa20','Pomp','Psmc5']):
        for off,scope,marker in [(-.12,'all_available','o'),(.12,'within_comparator_ranges','s')]:
            r=summary[(summary.target==t)&(summary.scope==scope)].iloc[0]
            if pd.notna(r.mean_residual):
                ax.plot(r.mean_residual,i+off,marker,color=colors[t])
                if pd.notna(r.ci_low):ax.hlines(i+off,r.ci_low,r.ci_high,color=colors[t])
    ax.axvline(0,color='black',lw=.6,ls='--');ax.set(yticks=range(3),yticklabels=['Naa20','Pomp','Psmc5'],xlabel='Litter-weighted adjusted residual',title='C  Descriptive 95% intervals',ylim=(2.6,-.6))
    ax.text(0,-.21,'Circle: all available; square: within ranges',transform=ax.transAxes,fontsize=8)
    fig.suptitle('Global-disruption adjustment: descriptive specificity check',fontsize=14)
    fig.savefig(OUT/'global_disruption.png',dpi=300);fig.savefig(OUT/'global_disruption.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)

if __name__=='__main__':main()
