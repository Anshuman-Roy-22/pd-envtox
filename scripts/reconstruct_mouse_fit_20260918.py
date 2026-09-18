#!/usr/bin/env python3
"""Fit only after the comparator manifest has been committed. No result-driven comparator replacement."""
import argparse,json,gzip,hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd,pyarrow.parquet as pq
from scipy import stats
from reconstruct_mouse_20260918 import ROOT,DOC,OUT,FOCAL,save,dump,vector,sha

def test(x):
 x=np.asarray(x,dtype=float);n=len(x)
 if n<2:return dict(n=n,mean=float(x.mean()) if n else None,se=None,ci_low=None,ci_high=None,p_two=None,p_down=None,lomo_low=None,lomo_high=None)
 m=x.mean();se=stats.sem(x);v=m/se;p=2*stats.t.sf(abs(v),n-1);half=stats.t.ppf(.975,n-1)*se;lo=(x.sum()-x)/(n-1)
 return dict(n=n,mean=float(m),se=float(se),ci_low=float(m-half),ci_high=float(m+half),p_two=float(p),p_down=float(stats.t.cdf(v,n-1)),lomo_low=float(lo.min()),lomo_high=float(lo.max()))

def scale_check(cache):
 pf=pq.ParquetFile(cache/'first_raw.parquet');raw=pf.read_row_group(0,columns=['cell_id','genes','expressions','num_rna_umi']).to_pandas()
 with gzip.open(cache/'barcodes.tsv.gz','rt') as f:barcodes=pd.Index([line.strip() for line in f])
 ix=barcodes.get_indexer(raw.cell_id);raw=raw[ix>=0];ix=ix[ix>=0]
 token=int(pq.read_table(cache/'gene_metadata.parquet').to_pandas().query('gene_name=="Xkr4"').gene_token_id.iloc[0])
 counts=np.array([v[np.flatnonzero(g==token)[0]] if (g==token).any() else 0 for g,v in zip(raw.genes,raw.expressions)])
 expected=np.log1p(10000*counts/raw.num_rna_umi.to_numpy());actual=vector(cache,'Xkr4')[ix];error=float(np.max(abs(expected-actual)));assert error<1e-5
 dump(OUT/'scale_validation.json',dict(matched_cells=len(raw),nonzero_cells=int((counts>0).sum()),maximum_absolute_error=error,scale='natural log(1 + 10000 * raw gene UMI / total RNA UMI)',raw_sha256=sha(cache/'first_raw.parquet')))

def scores(cache,targets):
 d=pd.read_pickle(cache/'eligible_metadata.pkl.gz');d=d[d.gene_target.isin(targets+['Non_target'])].reset_index(drop=True);pos=d.position.to_numpy();panels=pd.read_csv(DOC/'panel_mapping.tsv',sep='\t');panels=panels[panels.included]
 genes=sorted(set(panels.mouse)|{'Syp','Syt1','Snca','Naa20','Naa25','Pomp','Psmb4','Psmc5'})
 def acquire(g):vector(cache,g);return g
 for i,g in enumerate(ThreadPoolExecutor(5).map(acquire,genes)):
  if i%25==0:print('Expression vectors',i+1,'/',len(genes),flush=True)
 arrays={}
 for g in genes:arrays[g]=vector(cache,g)[pos].astype(float)
 for name,q in panels.groupby('panel'):d[name]=np.mean([arrays[g] for g in q.mouse],axis=0)
 remaining=[g for g in panels[panels.panel=='synapse'].mouse if g not in ['Syp','Syt1']];d['synapse_without_Syt1']=np.mean([arrays[g] for g in remaining],axis=0)
 for g in ['Syp','Syt1','Snca','Naa20','Naa25','Pomp','Psmb4','Psmc5']:d['gene_'+g]=arrays[g]
 dump(DOC/'response_vector_checksums.json',[dict(gene=g,bytes=(cache/'vectors'/(g+'.bin')).stat().st_size,sha256=sha(cache/'vectors'/(g+'.bin'))) for g in genes+['Xkr4']])
 return d

def effects(d,targets):
 readouts=['synapse','proteasome','synapse_without_Syt1']+['gene_'+g for g in ['Syp','Syt1','Snca','Naa20','Naa25','Pomp','Psmb4','Psmc5']]
 control=d[d.gene_target=='Non_target'].groupby(['source','predicted_group'],observed=True)[readouts].mean()
 elig=pd.read_csv(OUT/'mouse_eligibility.tsv',sep='\t');rows=[];guide=[]
 for target in targets:
  q=d[(d.gene_target==target)&d.retained];group=q.groupby(['source','predicted_group'],observed=True)
  mean=group[readouts].mean();count=group.size();delta=mean-control.reindex(mean.index);weighted=delta.mul(count,axis=0);eff=weighted.groupby(level=0,observed=True).sum().div(count.groupby(level=0,observed=True).sum(),axis=0)
  for mouse,v in eff.iterrows():
   e=elig[(elig.gene_target==target)&(elig.source==mouse)].iloc[0]
   for readout,value in v.items():rows.append(dict(target=target,mouse=str(mouse),readout=readout,difference=float(value),primary=bool(e.primary),min10=bool(e.min10),target_cells=int(e.retained_cells)))
  for (source,g),n in q.groupby(['source','guide_call'],observed=True).size().items():guide.append(dict(target=target,mouse=str(source),guide=str(g),cells=int(n)))
 e=pd.DataFrame(rows);save(e,OUT/'all_mouse_effects.tsv');save(pd.DataFrame(guide),OUT/'guide_coverage.tsv');return e

def calibration(d,cache):
 ntc=d[d.gene_target=='Non_target'].copy();guides=sorted(ntc.guide_call.astype(str).unique());rng=np.random.default_rng(20260916);shuffled=rng.permutation(guides);n=len(shuffled)//4;mapping={g:i//4 for i,g in enumerate(shuffled[:n*4])};ntc['pseudo']=ntc.guide_call.astype(str).map(mapping)
 save(pd.DataFrame([dict(guide=g,pseudo=i) for g,i in mapping.items()]),OUT/'calibration_groups.tsv')
 orig=pd.read_parquet(cache/'original_guides.parquet');orig=orig[orig.gene_target=='Non_target'].copy();orig['pseudo']=orig.guide_call.astype(str).map(mapping);before=orig.dropna(subset=['pseudo']).groupby(['pseudo','source'],observed=True).cells.sum()
 totals=ntc.groupby(['source','predicted_group'],observed=True).synapse.agg(['sum','count']);g=ntc.dropna(subset=['pseudo']).groupby(['pseudo','source','predicted_group'],observed=True).synapse.agg(['sum','count']).reset_index();g=g.join(totals,on=['source','predicted_group'],rsuffix='_all')
 g['control_n']=g['count_all']-g['count'];g=g[g.control_n>=20].copy();g['effect']=g['sum']/g['count']-(g['sum_all']-g['sum'])/g.control_n;g['weighted']=g.effect*g['count'];m=g.groupby(['pseudo','source'],observed=True).agg(weighted=('weighted','sum'),cells=('count','sum'));m=m.join(before.rename('original'));m['retention']=m.cells/m.original;m['effect']=m.weighted/m.cells;m=m[(m.cells>=20)&(m.retention>=.8)].reset_index();save(m,OUT/'calibration_mouse_effects.tsv')
 tests=[]
 for i in range(n):
  x=m[m.pseudo==i].effect;z=test(x);z['pseudo']=i;z['evaluable']=len(x)>=8;tests.append(z)
 t=pd.DataFrame(tests);save(t,OUT/'calibration_tests.tsv');valid=t[t.evaluable];positive=int((valid.p_down<.05).sum());bp=stats.binomtest(positive,len(valid),.05,alternative='greater').pvalue
 result=dict(groups=n,evaluable_groups=len(valid),nominal_down_positives=positive,positive_fraction=positive/len(valid),binomial_excess_p=bp,ignored_incomplete_guides=len(guides)%4,status='CALIBRATED' if len(valid)>=50 and bp>=.05 else 'NOT_CALIBRATED',caveat='Pseudo-tests share controls and are dependent; diagnostic only');dump(OUT/'calibration_summary.json',result);return result

def summarize(e):
 rows=[]
 for name in ['primary','min10']:
  for (target,readout),d in e[(e.target.isin(FOCAL))&e[name]].groupby(['target','readout']):rows.append(dict(analysis=name,target=target,readout=readout,**test(d.difference)))
 t=pd.DataFrame(rows);t['BH_four_down']=np.nan;t['BH_eight_two']=np.nan
 for name in ['primary','min10']:
  q=(t.analysis==name)&(t.readout=='synapse');t.loc[q,'BH_four_down']=stats.false_discovery_control(t.loc[q,'p_down'].to_numpy())
  q=(t.analysis==name)&t.readout.isin(['synapse','proteasome']);t.loc[q,'BH_eight_two']=stats.false_discovery_control(t.loc[q,'p_two'].to_numpy())
 save(t,OUT/'focal_tests.tsv');return t

def specificity(e,cache):
 manifest=pd.read_csv(DOC/'locked_comparators.tsv',sep='\t');allrows=[];tests=[];loo=[]
 animals=pd.read_excel(cache/'media-3.xlsx',header=1)[['Animal','Litter #']].dropna();litter=dict(zip('mouse'+animals.Animal.astype(int).astype(str),animals['Litter #']))
 littertests=[];litterrows=[]
 for focal in FOCAL:
  q=manifest[manifest.focal==focal];genes=q.mouse.to_list();f=e[(e.target==focal)&(e.readout=='synapse')&e.primary].set_index('mouse').difference
  ref=e[e.target.isin(genes)&(e.readout=='synapse')&e.primary].pivot(index='mouse',columns='target',values='difference').reindex(f.index)
  keep=ref.notna().sum(1)>=5;z=f[keep]-ref[keep].mean(1);qualified=len(genes)>=10 and len(z)>=8
  row=dict(target=focal,comparators=len(genes),status='EVALUABLE' if qualified else 'NOT_EVALUABLE',**test(z if qualified else []));tests.append(row)
  for mouse,value in z.items():allrows.append(dict(target=focal,mouse=mouse,focal_effect=float(f[mouse]),comparator_mean=float(ref.loc[mouse].mean()),comparators_in_mouse=int(ref.loc[mouse].notna().sum()),paired_difference=float(value),litter=litter[mouse],primary_evaluable=qualified))
  if qualified:
   l=pd.DataFrame({'effect':z,'litter':[litter[x] for x in z.index]}).groupby('litter').effect.mean();littertests.append(dict(target=focal,**test(l)))
   for lit,v in l.items():litterrows.append(dict(target=focal,litter=lit,effect=v))
   for g in genes:loo.append(dict(target=focal,omitted=g,mean=float((f[keep]-ref.loc[keep].drop(columns=g,errors='ignore').mean(1)).mean())))
  # Prespecified broad essential-like challenge, descriptive because support may be poor.
  shortlist=pd.read_csv(DOC/'baseline_shortlist.tsv',sep='\t');ess=sorted(set(shortlist.loc[shortlist.percentile_plusAO<=.2,'mouse']))
  r=e[e.target.isin(ess)&(e.readout=='synapse')&e.primary].groupby('mouse').difference.agg(['mean','count']).reindex(f.index);v=f[r['count']>=5]-r.loc[r['count']>=5,'mean']
  save(pd.DataFrame({'mouse':v.index,'paired_difference':v.values}),OUT/f'{focal}_essential_challenge.tsv')
 t=pd.DataFrame(tests);t['BH_four_two']=stats.false_discovery_control(t.p_two.fillna(1).to_numpy());save(t,OUT/'specificity_tests.tsv');save(pd.DataFrame(allrows),OUT/'specificity_mouse_pairs.tsv');save(pd.DataFrame(loo),OUT/'specificity_leave_one_comparator_out.tsv');save(pd.DataFrame(litterrows),OUT/'specificity_litter_effects.tsv')
 l=pd.DataFrame(littertests)
 if len(l):l['BH_four_two']=stats.false_discovery_control(np.array([float(l[l.target==g].p_two.iloc[0]) if g in set(l.target) else 1 for g in FOCAL]))[[FOCAL.index(g) for g in l.target]]
 save(l,OUT/'specificity_litter_tests.tsv');return t

def main():
 p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--focal-only',action='store_true');a=p.parse_args()
 # The lock must exist before comparator expression is read.
 manifest=pd.read_csv(DOC/'locked_comparators.tsv',sep='\t');short=pd.read_csv(DOC/'baseline_shortlist.tsv',sep='\t');targets=FOCAL if a.focal_only else sorted(set(FOCAL)|set(manifest.mouse)|set(short.loc[short.percentile_plusAO<=.2,'mouse']))
 scale_check(a.cache);d=scores(a.cache,targets);e=effects(d,targets);cal=calibration(d,a.cache);t=summarize(e)
 print(t[(t.analysis=='primary')&t.readout.isin(['synapse','proteasome'])].to_string(index=False),flush=True);print(cal,flush=True)
 if not a.focal_only:print(specificity(e,a.cache).to_string(index=False),flush=True)
if __name__=='__main__':main()
