#!/usr/bin/env python3
"""Reconstruction stages: metadata, baseline matching, focal fit, then locked specificity."""
import argparse,gzip,hashlib,json,struct,zlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import Request,urlopen
import numpy as np,pandas as pd,pyarrow.parquet as pq
from scipy import stats
ROOT=Path(__file__).resolve().parents[1];DOC=ROOT/'docs/v2/reconstruction_20260918';OUT=ROOT/'results/v2/reconstruction_20260918/mouse'
FOCAL=['Naa20','Pomp','Psmb4','Psmc5'];UCSC='https://cells.ucsc.edu/whole-brain-perturb/combined/'
def save(d,p):d.to_csv(p,sep='\t',index=False,float_format='%.12g',na_rep='NA',lineterminator='\n')
def dump(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def orthology(cache):
 b=pd.read_csv(cache/'biomart.csv').dropna().drop_duplicates();b.columns=['human_id','human','mouse_id','mouse']
 b=b[(b.groupby('human_id').mouse_id.transform('nunique')==1)&(b.groupby('mouse_id').human_id.transform('nunique')==1)]
 b=b[(b.groupby('human').mouse.transform('nunique')==1)&(b.groupby('mouse').human.transform('nunique')==1)]
 return b.drop_duplicates(['human','mouse']).sort_values(['human','mouse'])
def metadata(cache):
 OUT.mkdir(parents=True,exist_ok=True)
 columns=['cell_id','source','gene_target','guide_call','passes_qc','num_guides','neuron_type','predicted_group']
 parts=[]
 pf=pq.ParquetFile(cache/'all_obs.parquet')
 for i in range(pf.num_row_groups):
  q=pf.read_row_group(i,columns=columns).to_pandas();q=q[q.passes_qc&(q.num_guides==1)&(q.neuron_type!='Non-Neuron')];parts.append(q.drop(columns=['passes_qc','num_guides','neuron_type']))
 d=pd.concat(parts,ignore_index=True);del parts
 # Safe-target controls have collapsed browser labels and are not molecular comparator genes.
 d=d[~d.gene_target.str.startswith('Safe_target_')].copy()
 assert d.cell_id.is_unique
 for c in ['source','gene_target','guide_call','predicted_group']:d[c]=d[c].astype('category')
 original=d.groupby(['gene_target','source'],observed=True).size().rename('original_cells').reset_index()
 origguide=d.groupby(['gene_target','source','guide_call'],observed=True).size().rename('cells').reset_index();origguide.to_parquet(cache/'original_guides.parquet',index=False)
 with gzip.open(cache/'barcodes.tsv.gz','rt') as f:barcodes=pd.Index([line.strip() for line in f])
 assert barcodes.is_unique
 pos=barcodes.get_indexer(d.cell_id);d['position']=pos
 missing=d[d.position<0];save(missing.astype({'gene_target':str,'source':str,'guide_call':str,'predicted_group':str}),OUT/'missing_browser_cells.tsv');d=d[d.position>=0].copy()
 config=json.loads((cache/'dataset.json').read_text())
 for c in ['source','gene_target','predicted_group']:
  field=next(x for x in config['metaFields'] if x['name']==c);dtype={'Uint8':'u1','Uint16':'<u2'}[field['arrType']];values=np.frombuffer(gzip.decompress((cache/('browser_'+c+'.bin.gz')).read_bytes()),dtype=dtype)
  labels=np.array([x[0] for x in field['valCounts']],dtype=object);mapped=labels[values[d.position.to_numpy()]];assert (mapped==d[c].astype(str).to_numpy()).all(),c
 cnt=d[d.gene_target=='Non_target'].groupby(['source','predicted_group'],observed=True).size().rename('control_cells').reset_index()
 d=d.merge(cnt,on=['source','predicted_group'],how='left');d['retained']=d.control_cells.fillna(0)>=20
 kept=d[d.retained].groupby(['gene_target','source'],observed=True).size().rename('retained_cells').reset_index()
 elig=original.merge(kept,on=['gene_target','source'],how='left').fillna({'retained_cells':0});elig['retention']=elig.retained_cells/elig.original_cells;elig['primary']=(elig.retained_cells>=20)&(elig.retention>=.8);elig['min10']=(elig.retained_cells>=10)&(elig.retention>=.8)
 save(elig,OUT/'mouse_eligibility.tsv')
 d.to_pickle(cache/'eligible_metadata.pkl.gz',compression='gzip')
 assert len(pd.read_pickle(cache/'eligible_metadata.pkl.gz'))==len(d)
 b=orthology(cache);save(b,DOC/'orthology.tsv')
 index=json.loads((cache/'exprMatrix.json').read_text());human=pd.read_csv(ROOT/'results/v2/reconstruction_20260918/human/pathway_membership.tsv',sep='\t')
 panels=[]
 for pid,label in [('R-HSA-112315','synapse'),('R-HSA-9907900','proteasome')]:
  s=human[human.pathway==pid].merge(b,on='human',how='left') if 'human'in human.columns else human[human.pathway==pid].rename(columns={'gene':'human'}).merge(b,on='human',how='left')
  s['panel']=label;s['included']=s.mouse.isin(index);panels.append(s)
 p=pd.concat(panels,ignore_index=True);save(p,DOC/'panel_mapping.tsv')
 for label,g in p.groupby('panel'):assert g.included.sum()/len(g)>=.7
 dep=pd.read_csv(ROOT/'results/v2/neuronal_survival/CRISPRi_gene_level.tsv',sep='\t').rename(columns={'gene':'human'})
 dep=b.merge(dep[['human','phenotype_plusAO','percentile_plusAO']],on='human')
 excluded={'NAA20','NAA25','POMP','PSMA4','PSMB4','PSMC5','PSMD1','PSMD2','PSMD4','PSMG1'}
 for line in (ROOT/'data_raw/pathways/ReactomePathways.gmt').read_text().splitlines():
  name,pid,*genes=line.split('\t')
  if 'proteasom' in name.lower():excluded.update(genes)
 qualified=elig[elig.primary].groupby('gene_target',observed=True).size().rename('n_mice')
 guides=d[d.retained].groupby('gene_target',observed=True).guide_call.nunique().rename('n_guides')
 candidates=dep.merge(qualified,left_on='mouse',right_index=True).merge(guides,left_on='mouse',right_index=True)
 candidates=candidates[(candidates.n_mice>=8)&(candidates.n_guides==4)&candidates.mouse.isin(index)].copy();candidates['excluded_family']=candidates.human.isin(excluded)
 save(candidates,DOC/'all_candidate_covariates.tsv')
 short=[]
 for focal in FOCAL:
  row=candidates[candidates.mouse==focal].iloc[0];g=candidates[~candidates.excluded_family].copy();g['dependency_delta']=g.percentile_plusAO-row.percentile_plusAO;g['coverage_log2_ratio']=np.log2(g.n_mice/row.n_mice)
  g=g[(abs(g.dependency_delta)<=.10)&(abs(g.coverage_log2_ratio)<=1)];g['distance']=(g.dependency_delta/.10)**2+g.coverage_log2_ratio**2;g=g.sort_values(['distance','mouse']).head(30);g['focal']=focal;short.append(g)
 s=pd.concat(short,ignore_index=True);save(s,DOC/'baseline_shortlist.tsv')
 dump(OUT/'metadata_reconstruction.json',{'qualified_all_targets':int(original.original_cells.sum()),'matched_all_targets':len(d),'missing_all_targets':len(missing),'focal_control_original':int(original[original.gene_target.isin(FOCAL+['Non_target'])].original_cells.sum()),'focal_control_matched':int(d.gene_target.isin(FOCAL+['Non_target']).sum()),'metadata_disagreements':0,'mapped_panels':{str(k):int(g.included.sum()) for k,g in p.groupby('panel')},'shortlist_counts':{str(k):len(g) for k,g in s.groupby('focal')}})
 print('Focal eligibility',candidates[candidates.mouse.isin(FOCAL)][['mouse','n_mice','n_guides','percentile_plusAO']].to_string(index=False),flush=True);print(s.groupby('focal').size().to_string(),flush=True)
def vector(cache,g):
 idx=json.loads((cache/'exprMatrix.json').read_text());offset,size=idx[g];folder=cache/'vectors';folder.mkdir(exist_ok=True);path=folder/(g+'.bin')
 if not path.exists():
  req=Request(UCSC+'exprMatrix.bin?'+g,headers={'Range':f'bytes={offset}-{offset+size-1}'})
  for attempt in range(3):
   try:
    with urlopen(req,timeout=90) as f:
     assert f.status==206,'Range not honored';data=f.read();assert len(data)==size
    path.write_bytes(data);break
   except Exception:
    if attempt==2:raise
 raw=zlib.decompress(path.read_bytes());n=struct.unpack('<H',raw[:2])[0];desc=raw[2:2+n].decode();assert g in desc,(g,desc)
 x=np.frombuffer(raw,offset=2+n,dtype='<f4');assert len(x)==6348631
 return x

def baseline(cache):
 d=pd.read_pickle(cache/'eligible_metadata.pkl.gz');ntc=d[(d.gene_target=='Non_target')&d.retained].copy();short=pd.read_csv(DOC/'baseline_shortlist.tsv',sep='\t');genes=sorted(set(short.mouse)|set(FOCAL));results=[]
 def one(g):
  x=vector(cache,g)[ntc.position.to_numpy()];means=pd.Series(x).groupby(ntc.source.astype(str).to_numpy()).mean();return {'mouse':g,'baseline_expression':float(means.mean()),'control_mice':len(means)}
 for i,r in enumerate(ThreadPoolExecutor(5).map(one,genes)):
  results.append(r)
  if i%20==0:print('Baseline targets',i+1,'/',len(genes),flush=True)
 b=pd.DataFrame(results);save(b,DOC/'baseline_expression.tsv');joined=short.merge(b,on='mouse');bm=b.set_index('mouse').baseline_expression;matched=[];decisions=[]
 for focal in FOCAL:
  q=joined[joined.focal==focal].copy();q=q[q.baseline_expression>0];q['expression_log2_ratio']=np.log2(q.baseline_expression/bm[focal]);q=q[abs(q.expression_log2_ratio)<=1];q['matching_distance']=(q.dependency_delta/.1)**2+q.coverage_log2_ratio**2+q.expression_log2_ratio**2;q=q.sort_values(['matching_distance','mouse']).head(20);q['support']=len(q)>=10;matched.append(q);decisions.append({'focal':focal,'matched_genes':len(q),'baseline_expression':float(bm[focal]),'status':'QUALIFIED' if len(q)>=10 else 'NOT_EVALUABLE'})
 save(pd.concat(matched,ignore_index=True),DOC/'locked_comparators.tsv');dump(DOC/'comparator_qualification.json',decisions)
 manifest=[dict(gene=g,bytes=(cache/'vectors'/(g+'.bin')).stat().st_size,sha256=sha(cache/'vectors'/(g+'.bin'))) for g in genes];dump(DOC/'baseline_vector_checksums.json',manifest);print(decisions,flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['metadata','baseline']);p.add_argument('--cache',type=Path,required=True);a=p.parse_args();globals()[a.stage](a.cache)
if __name__=='__main__':main()
