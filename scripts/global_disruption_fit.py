"""Frozen comparator-only REML model and whole-litter descriptive bootstrap."""
import hashlib
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg, optimize, stats
from global_disruption_io import dump
from global_disruption_metrics import OUT, ROOT, save

PRED=['global_rms','realized_target_difference','log_target_cells']
FOCAL=['Naa20','Pomp','Psmc5']

class Unevaluable(RuntimeError): pass

def design(d,levels=None,center=None,scale=None):
    if levels is None:levels=sorted(d.mouse.unique())
    continuous=d[PRED].to_numpy(dtype=float)
    if center is None:center=continuous.mean(0)
    if scale is None:scale=continuous.std(0)
    if not np.isfinite(continuous).all() or (scale<=0).any():raise Unevaluable('Non-finite or constant predictor')
    X=np.column_stack([np.ones(len(d)),(continuous-center)/scale]+[(d.mouse==m).to_numpy(dtype=float) for m in levels[1:]])
    return X,levels,center,scale

def reml(y,X,groups):
    y=np.asarray(y,dtype=float);X=np.asarray(X,dtype=float)
    if not np.isfinite(X).all() or not np.isfinite(y).all():raise Unevaluable('Non-finite model inputs')
    n,p=X.shape
    if n<=p or np.linalg.matrix_rank(X)!=p:raise Unevaluable('Fixed design is rank deficient')
    names,codes=np.unique(groups,return_inverse=True);k=len(names);ng=np.bincount(codes)
    sx=np.zeros((k,p));sy=np.zeros(k);np.add.at(sx,codes,X);np.add.at(sy,codes,y)
    xx=X.T@X;xy=X.T@y;yy=y@y
    def calc(lam,return_fit=False):
        w=lam/(1+lam*ng);vx=xx-sx.T@(w[:,None]*sx);vy=xy-sx.T@(w*sy)
        try:
            cf=linalg.cho_factor(vx,lower=True,check_finite=False)
            beta=linalg.cho_solve(cf,vy,check_finite=False)
        except linalg.LinAlgError:
            if return_fit:raise Unevaluable('Singular REML crossproduct')
            return np.inf
        rss=yy-np.dot(w,sy*sy)-vy@beta
        if rss<=0 or not np.isfinite(rss):
            if return_fit:raise Unevaluable('Nonpositive REML residual variance')
            return np.inf
        var=rss/(n-p)
        objective=(n-p)*np.log(var)+np.log1p(lam*ng).sum()+2*np.log(np.diag(cf[0])).sum()
        if not return_fit:return objective
        return {'beta':beta,'residual_variance':float(var),'gene_variance':float(lam*var),'lambda':float(lam),'objective':float(objective),'se':np.sqrt(np.diag(linalg.cho_solve(cf,np.eye(p),check_finite=False))*var)}
    opt=optimize.minimize_scalar(lambda z:calc(np.exp(z)),bounds=(-20,20),method='bounded',options={'xatol':1e-9})
    if not opt.success or opt.x>19.9:raise Unevaluable('REML optimizer failed or reached upper boundary')
    lam=np.exp(opt.x) if opt.fun<calc(0) else 0.
    return calc(lam,True)

def fit(train):
    if train.target.nunique()<20:raise Unevaluable('Fewer than 20 represented comparator genes')
    X,levels,center,scale=design(train)
    model=reml(train.difference.to_numpy(),X,train.target.to_numpy())
    model.update(levels=levels,center=center,scale=scale,n=len(train),genes=train.target.nunique())
    return model

def validate_statsmodels(train,model):
    import statsmodels.api as sm
    X,*_=design(train,model['levels'],model['center'],model['scale'])
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter('always')
        result=sm.MixedLM(train.difference.to_numpy(),X,groups=train.target.to_numpy()).fit(reml=True,method=['lbfgs','powell'],disp=False)
    if not result.converged:raise Unevaluable('Independent statsmodels validation did not converge')
    err=float(np.max(np.abs(X@model['beta']-X@result.fe_params)))
    gv=float(np.asarray(result.cov_re)[0,0]);rv=float(result.scale)
    if err>1e-6 or abs(rv-model['residual_variance'])>max(1e-10,rv*.001) or abs(gv-model['gene_variance'])>max(1e-8,gv*.001):
        raise Unevaluable(f'Independent REML validation mismatch: prediction error {err}')
    return {'maximum_fixed_prediction_difference':err,'statsmodels_gene_variance':gv,'statsmodels_residual_variance':rv,'warnings':[str(x.message) for x in recorded]}

def predict(d,train,model):
    q=d.copy();q['prediction']=np.nan;q['residual']=np.nan
    q['known_mouse']=q.mouse.isin(model['levels'])
    q['evaluable_prediction']=q.measurable&q.known_mouse
    q['extrapolation']=False
    for col in PRED[:2]:
        q[col+'_outside_range']=(q[col]<train[col].min())|(q[col]>train[col].max())
        q['extrapolation']|=q[col+'_outside_range']
    valid=q.evaluable_prediction
    if valid.any():
        X,*_=design(q[valid],model['levels'],model['center'],model['scale'])
        q.loc[valid,'prediction']=X@model['beta']
        q.loc[valid,'residual']=q.loc[valid,'difference']-q.loc[valid,'prediction']
    return q

def summaries(p):
    rows=[]
    for target in FOCAL:
        q=p[(p.target==target)&p.evaluable_prediction]
        for scope,sub in [('all_available',q),('within_comparator_ranges',q[~q.extrapolation])]:
            v=sub.groupby(['mouse','litter']).residual.mean().reset_index().groupby('litter').residual.mean()
            rows.append({'target':target,'scope':scope,'mean_residual':float(v.mean()) if len(v) else None,'mice':sub.mouse.nunique(),'litters':len(v),'extrapolated_mice':sub.loc[sub.extrapolation,'mouse'].nunique()})
    return rows

def resample(d,draw):
    pieces=[]
    for i,lit in enumerate(draw):
        p=d[d.litter==lit].copy()
        p['mouse']=p.mouse.astype(str)+'__copy'+str(i)
        p['litter']=str(lit)+'__copy'+str(i)
        pieces.append(p)
    return pd.concat(pieces,ignore_index=True)

def main():
    d=pd.read_csv(OUT/'global_metrics.tsv',sep='\t',dtype={'mouse':str,'litter':str})
    reconciliation=json.loads((OUT/'source_reconciliation.json').read_text());assert reconciliation['status']=='PASS'
    controls=d[(d.role=='comparator')&d.measurable]
    counts=controls.groupby('target').mouse.nunique();qualified=sorted(counts[counts>=3].index)
    training=controls[controls.target.isin(qualified)].copy()
    save(counts.rename('measurable_mice').reset_index().assign(qualified=lambda q:q.target.isin(qualified)),OUT/'comparator_measurement_eligibility.tsv')
    try:
        model=fit(training)
        validation=validate_statsmodels(training,model)
    except Unevaluable as error:
        dump(OUT/'model_status.json',{'status':'NOT_EVALUABLE','reason':str(error),'original_matched_result':'Unchanged; no replacement model'})
        return
    dump(OUT/'observed_model.json',{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in model.items()})
    dump(OUT/'observed_model_validation.json',validation)
    labels=['intercept']+[c+'_per_training_SD' for c in PRED]+['mouse_'+m for m in model['levels'][1:]]
    save(pd.DataFrame({'term':labels,'coefficient':model['beta'],'model_SE':model['se']}),OUT/'coefficients.tsv')
    focal=d[d.target.isin(FOCAL)].copy();pred=predict(focal,training,model);save(pred,OUT/'focal_predictions.tsv')
    # Comparator and focal population predictions consistently use gene intercept zero.
    comp=predict(training,training,model);save(comp,OUT/'comparator_predictions.tsv')
    obs=summaries(pred)
    bootdata=d[d.target.isin(qualified+FOCAL)].copy();lits=sorted(bootdata.litter.unique())
    digest=hashlib.sha256((OUT/'global_metrics.tsv').read_bytes()+Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    checkpoint=OUT/'bootstrap_checkpoint.json';records=[]
    if checkpoint.exists():
        old=json.loads(checkpoint.read_text())
        if old['input_and_code_sha256']!=digest:raise RuntimeError('Bootstrap checkpoint belongs to different inputs/code')
        records=old['records']
        assert [x['draw'] for x in records]==list(range(len(records)))
    rng=np.random.default_rng(20260918)
    for i in range(1000):
        draw=rng.choice(lits,size=len(lits),replace=True)
        if i<len(records):continue
        b=resample(bootdata,draw);tr=b[(b.role=='comparator')&b.measurable]
        try:
            m=fit(tr);p=predict(b[b.target.isin(FOCAL)],tr,m)
            row={'draw':i,'status':'FIT','summaries':summaries(p)}
        except (Unevaluable,np.linalg.LinAlgError) as error:
            row={'draw':i,'status':'FAILED','reason':str(error),'summaries':[]}
        records.append(row)
        if (i+1)%10==0:
            dump(checkpoint,{'input_and_code_sha256':digest,'records':records})
            print('Litter bootstrap',i+1,'/1000; successful',sum(r['status']=='FIT' for r in records),flush=True)
    valid=sum(r['status']=='FIT' for r in records)
    allrows=[]
    for r in records:
        if r['summaries']:
            allrows.extend(dict(draw=r['draw'],fit_status=r['status'],**v) for v in r['summaries'])
        else:allrows.append({'draw':r['draw'],'fit_status':r['status'],'failure':r['reason']})
    save(pd.DataFrame(allrows),OUT/'bootstrap_draws.tsv')
    for row in obs:
        values=[v['mean_residual'] for r in records for v in r['summaries'] if v['target']==row['target'] and v['scope']==row['scope'] and v['mean_residual'] is not None]
        row['available_bootstrap_summaries']=len(values)
        row['ci_low']=row['ci_high']=None
        row['interval_status']='NOT_EVALUABLE'
        if valid>=900 and len(values)>=900 and row['mean_residual'] is not None:
            row['ci_low'],row['ci_high']=map(float,np.quantile(values,[.025,.975]))
            row['interval_status']='DESCRIPTIVE_PERCENTILE'
    save(pd.DataFrame(obs),OUT/'focal_residual_summary.tsv')
    dump(OUT/'model_status.json',{'status':'FITTED_DESCRIPTIVE','comparator_genes':len(qualified),'training_observations':len(training),'bootstrap_successful':valid,'bootstrap_draws':1000,'interval_requirement':900,'scope':'Sensitivity analysis; does not exclude global-disruption confounding, prove viability or establish causal mediation.'})

if __name__=='__main__':main()
