#!/usr/bin/env python3
"""Reconstruct the known NatB audit and relative-specificity tests from preserved inputs."""
import argparse,itertools,json,subprocess,os
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from natb_selectivity import prepare,acquire,RUNS,save,difference
from natb_robustness import balanced_proteome,two_sided_exact
ROOT=Path(__file__).resolve().parents[1]
FOCUS=['R-HSA-112315','R-HSA-112310','R-HSA-9907900','R-HSA-611105','R-HSA-72766','R-HSA-9612973','R-HSA-397014']
def main():
 a=argparse.ArgumentParser();a.add_argument('--output',type=Path,default=ROOT/'results/v2/reconstruction_20260918/human');args=a.parse_args();o=args.output;o.mkdir(parents=True,exist_ok=True)
 acquire(ROOT/'data_raw/natb_selectivity',False)
 raw,norm,cov,offset,features=prepare(ROOT/'data_raw/natb_selectivity/adj4767_Table_S5.csv',o,return_features=True)
 old=pd.read_csv(ROOT/'results/v2/natb_selectivity/protein_sample_scores.tsv',sep='\t').set_index('gene')
 assert np.allclose(norm,old,equal_nan=True,atol=1e-10)
 _,strict,_,_=balanced_proteome(features,True);save(strict.rename_axis('gene').reset_index(),o/'strict_matrix.tsv')
 complete=norm.dropna();sets={};names={}
 for line in (ROOT/'data_raw/pathways/ReactomePathways.gmt').read_text().splitlines():
  name,pid,*g=line.split('\t');sets[pid]=set(g);names[pid]=name
 rows=[];members=[];samples=[]
 for pid in FOCUS:
  genes=sorted(set(complete.index)&sets[pid]);vals=complete.loc[genes].mean().to_numpy();rows.append(dict(pathway=pid,n_genes=len(genes),**difference(vals),p_exact=two_sided_exact(vals)))
  members.extend(dict(pathway=pid,gene=g) for g in genes)
  samples.extend(dict(pathway=pid,run=r,score=v) for r,v in zip(RUNS,vals))
 t=pd.DataFrame(rows);t['BH_seven']=stats.false_discovery_control(t.p_exact.to_numpy());save(t,o/'pathway_score_tests.tsv');save(pd.DataFrame(members),o/'pathway_membership.tsv');save(pd.DataFrame(samples),o/'pathway_culture_scores.tsv')
 # Label-blind detectability covariates, keeping original precursor filters.
 d=pd.read_csv(ROOT/'data_raw/natb_selectivity/adj4767_Table_S5.csv',dtype=str)
 def single(v):
  if not isinstance(v,str):return None
  g={s for s in v.split(';') if s};return next(iter(g)) if len(g)==1 else None
 d['gene']=d.Genes.map(single);d=d[d.gene.notna()].drop_duplicates(['gene','EG.PrecursorId'])
 q=d[[f'{r}.raw.EG.Qvalue' for r in RUNS]].apply(pd.to_numeric,errors='coerce').to_numpy();x=d[[f'{r}.raw.EG.TotalQuantity (Settings)' for r in RUNS]].apply(pd.to_numeric,errors='coerce').to_numpy()
 valid=np.isfinite(q)&(q<=.01)&np.isfinite(x)&(x>0);keep=valid.sum(1)>=8;log=np.full_like(x[keep],np.nan);log[valid[keep]]=np.log2(x[keep][valid[keep]])
 ab=pd.DataFrame({'gene':d.loc[keep,'gene'],'abundance':np.nanmedian(log,axis=1)}).groupby('gene').abundance.median()
 c=pd.DataFrame({'pooled_abundance':ab,'backbones':cov.peptide_backbones}).loc[complete.index];c['abundance_rank']=c.pooled_abundance.rank(pct=True);c['backbone_rank']=np.log2(c.backbones).rank(pct=True)
 save(c.rename_axis('gene').reset_index(),o/'matching_covariates.tsv')
 syn=sets[FOCUS[0]]|sets[FOCUS[1]];background=sorted(set(complete.index)-syn-{'NAA25','SNCA'});matches=[];tests=[];culture=[]
 for pid in FOCUS[:2]:
  genes=sorted(set(complete.index)&sets[pid]);panel=complete.loc[genes].mean().to_numpy();ref=[]
  for g in genes:
   dist=((c.loc[background,['abundance_rank','backbone_rank']]-c.loc[g,['abundance_rank','backbone_rank']])**2).sum(1)
   near=pd.DataFrame({'gene':background,'distance':dist.to_numpy()}).sort_values(['distance','gene']).head(5)
   for z in near.itertuples():matches.append(dict(pathway=pid,gene=g,control=z.gene,distance=z.distance))
   ref.append(complete.loc[near.gene].mean().to_numpy())
  for kind,reference in [('all_nonsynaptic',complete.loc[background].mean().to_numpy()),('matched_nonsynaptic',np.mean(ref,axis=0))]:
   val=panel-reference;tests.append(dict(pathway=pid,comparison=kind,**difference(val),p_exact=two_sided_exact(val)))
   for r,p,b,v in zip(RUNS,panel,reference,val):culture.append(dict(pathway=pid,comparison=kind,run=r,panel=p,reference=b,relative=v))
   shifts=np.linspace(-3,3,10);assert np.max(abs((panel+shifts)-(reference+shifts)-val))<1e-12
 t=pd.DataFrame(tests);t['BH_four']=stats.false_discovery_control(t.p_exact.to_numpy());save(t,o/'relative_tests.tsv');save(pd.DataFrame(matches),o/'protein_matches.tsv');save(pd.DataFrame(culture),o/'relative_culture_scores.tsv')
 print(t[['pathway','comparison','log2_difference','ci_low','ci_high','p_exact','BH_four']].to_string(index=False))
if __name__=='__main__':main()
