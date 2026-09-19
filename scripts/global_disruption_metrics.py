"""Original-weight whole-transcriptome distances; no outcome-driven gene selection."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import h5py
from scipy import sparse
from global_disruption_io import download, dump, sha

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/'docs/v2/global_disruption_20260919'
OLD=ROOT/'docs/v2/reconstruction_20260918'
OLDOUT=ROOT/'results/v2/reconstruction_20260918/mouse'
OUT=ROOT/'results/v2/global_disruption_20260919'
FOCAL=['Naa20','Pomp','Psmb4','Psmc5']

def save(d,path):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(path,sep='\t',index=False,float_format='%.12g',lineterminator='\n',na_rep='NA')

def text_hash(path):
    return hashlib.sha256(Path(path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def read_col(group,key):
    obj=group[key]
    if isinstance(obj,h5py.Group):
        cats=read_col(obj,'categories');codes=obj['codes'][:]
        if (codes<0).any(): raise RuntimeError('Missing categorical values: '+key)
        return cats[codes]
    if obj.dtype.kind in 'OSU': return obj.asstr()[:]
    return obj[:]

def prepare(cache):
    OUT.mkdir(parents=True,exist_ok=True);cache.mkdir(parents=True,exist_ok=True)
    manifests=json.loads((OLD/'public_source_manifest.json').read_text())
    exclude={'first_raw.parquet','wholebrain_supp.zip','biomart.csv','gene_metadata.parquet'}
    for spec in manifests:
        if spec['file'] not in exclude:
            print('Verifying/acquiring',spec['file'],flush=True);download(spec,cache)
    wb=ROOT/'data_raw/global_disruption/media-3.xlsx'
    litter_spec=json.loads((DOC/'litter_source.json').read_text())
    assert sha(wb)==litter_spec['workbook_sha256'],'Author workbook checksum'
    a=pd.read_excel(wb,header=1)[['Animal','Litter #']].dropna()
    assert a.Animal.is_unique
    litter=dict(zip('mouse'+a.Animal.astype(int).astype(str),a['Litter #'].astype(str)))
    for file in ['specificity_mouse_pairs.tsv']:
        prev=pd.read_csv(OLDOUT/file,sep='\t')
        assert np.allclose(prev.mouse.map(litter).astype(float),prev.litter),'Author litter mapping changed'
    locked=pd.read_csv(OLD/'locked_comparators.tsv',sep='\t')
    comparators=sorted(set(locked.mouse));targets=sorted(set(comparators)|set(FOCAL))
    assert len(comparators)==51 and len(targets)==55
    index=json.loads((cache/'exprMatrix.json').read_text())
    genes=sorted(g for g in index if g!='_range')
    panel=pd.read_csv(OLD/'panel_mapping.tsv',sep='\t')
    syn=panel.loc[(panel.panel=='synapse')&panel.included,'mouse'].tolist()
    assert len(syn)==85 and len(set(syn))==85
    background=sorted(set(genes)-set(syn)-set(targets))
    assert len(genes)==19070 and len(background)==18931
    save(pd.DataFrame({'gene':genes,'background':[g in set(background) for g in genes],'synapse':[g in syn for g in genes]}),OUT/'gene_universe.tsv')
    cols=['cell_id','source','gene_target','passes_qc','num_guides','neuron_type','predicted_group','num_rna_umi']
    parts=[];pf=pq.ParquetFile(cache/'all_obs.parquet')
    for i in range(pf.num_row_groups):
        q=pf.read_row_group(i,columns=cols).to_pandas()
        q=q[q.passes_qc&(q.num_guides==1)&(q.neuron_type!='Non-Neuron')&q.gene_target.isin(targets+['Non_target'])]
        parts.append(q[['cell_id','source','gene_target','predicted_group','num_rna_umi']])
    d=pd.concat(parts,ignore_index=True);del parts
    for c in ['cell_id','source','gene_target','predicted_group']:d[c]=d[c].astype(str)
    assert d.cell_id.is_unique
    original=d.groupby(['gene_target','source']).size()
    with gzip.open(cache/'barcodes.tsv.gz','rt') as f:barcodes=pd.Index(line.strip() for line in f)
    assert barcodes.is_unique
    d['position']=barcodes.get_indexer(d.cell_id);d=d[d.position>=0].copy()
    cfg=json.loads((cache/'dataset.json').read_text())
    for col in ['source','gene_target','predicted_group']:
        field=next(x for x in cfg['metaFields'] if x['name']==col)
        dt={'Uint8':'u1','Uint16':'<u2'}[field['arrType']]
        values=np.frombuffer(gzip.decompress((cache/f'browser_{col}.bin.gz').read_bytes()),dtype=dt)
        labels=np.array([x[0] for x in field['valCounts']],dtype=object)
        assert np.array_equal(labels[values[d.position.to_numpy()]],d[col].to_numpy()),col
    ntc=d[d.gene_target=='Non_target'].groupby(['source','predicted_group']).size().rename('control_cells')
    d=d.join(ntc,on=['source','predicted_group']);d=d[d.control_cells>=20].copy()
    retained=d.groupby(['gene_target','source']).size()
    eligibility=original.rename('original_cells').to_frame().join(retained.rename('retained_cells')).fillna(0)
    eligibility['primary']=(eligibility.retained_cells>=20)&(eligibility.retained_cells/eligibility.original_cells>=.8)
    eligible=eligibility.reset_index().query('primary and gene_target != "Non_target"')
    prev=pd.read_csv(OLDOUT/'mouse_eligibility.tsv',sep='\t')
    merged=eligible.merge(prev,on=['gene_target','source'],suffixes=('','_previous'),validate='one_to_one')
    assert len(merged)==len(eligible) and merged.primary_previous.all()
    assert (merged.retained_cells==merged.retained_cells_previous).all()
    assert (merged.original_cells==merged.original_cells_previous).all()
    expected=prev[prev.primary&prev.gene_target.isin(targets)]
    assert len(expected)==len(eligible),'Original eligibility set changed'
    keys=set(zip(eligible.gene_target,eligible.source))
    d=d[(d.gene_target=='Non_target')|pd.MultiIndex.from_frame(d[['gene_target','source']]).isin(keys)].copy()
    # Controls not sharing a stratum with any eligible target carry zero weight; omit safely.
    used_strata=pd.MultiIndex.from_frame(d[d.gene_target!='Non_target'][['source','predicted_group']].drop_duplicates())
    d=d[(d.gene_target!='Non_target')|pd.MultiIndex.from_frame(d[['source','predicted_group']]).isin(used_strata)].copy()
    d=d.sort_values('cell_id').reset_index(drop=True)
    groups=d.groupby(['gene_target','source','predicted_group'],sort=True).size().rename('n').reset_index()
    groups['group']=np.arange(len(groups))
    d=d.merge(groups,on=['gene_target','source','predicted_group'],validate='many_to_one').sort_values('cell_id').reset_index(drop=True)
    obs=eligible.rename(columns={'gene_target':'target','source':'mouse'}).sort_values(['target','mouse']).reset_index(drop=True)
    obs['litter']=obs.mouse.map(litter);assert obs.litter.notna().all()
    obs['role']=np.where(obs.target.isin(FOCAL),'focal','comparator')
    old=pd.read_csv(OLDOUT/'all_mouse_effects.tsv',sep='\t')
    old=old[old.primary&(old.readout=='synapse')][['target','mouse','difference','target_cells']]
    obs=obs.merge(old,on=['target','mouse'],validate='one_to_one')
    assert np.array_equal(obs.retained_cells.to_numpy(),obs.target_cells.to_numpy())
    lookup={(r.gene_target,r.source,r.predicted_group):(r.group,r.n) for r in groups.itertuples()}
    rows=[];columns=[];vals=[];cr=[];cc=[];cv=[]
    target_group=groups[groups.gene_target!='Non_target'].groupby(['gene_target','source'])
    for i,r in enumerate(obs.itertuples()):
        q=target_group.get_group((r.target,r.mouse));total=q.n.sum();assert total==r.target_cells
        for s in q.itertuples():
            ctl,nc=lookup[('Non_target',r.mouse,s.predicted_group)]
            rows.extend([i,i]);columns.extend([s.group,ctl]);vals.extend([1/total,-s.n/total/nc])
            cr.append(i);cc.append(ctl);cv.append(s.n/total/nc)
    A=sparse.csr_matrix((vals,(rows,columns)),shape=(len(obs),len(groups)))
    C=sparse.csr_matrix((cv,(cr,cc)),shape=A.shape)
    assert np.allclose(A@groups.n.to_numpy(),0) and np.allclose(C@groups.n.to_numpy(),1)
    d.to_parquet(cache/'selected_cells.parquet',index=False)
    groups.to_parquet(cache/'groups.parquet',index=False)
    sparse.save_npz(cache/'A.npz',A);sparse.save_npz(cache/'C.npz',C)
    save(obs,OUT/'eligible_observations.tsv');save(groups,OUT/'aggregation_groups.tsv')
    h=hashlib.sha256()
    for c in ['cell_id','group','num_rna_umi','position']:h.update(d[c].to_csv(index=False,lineterminator='\n').encode())
    for mat in [A,C]:
        for arr in [mat.data,mat.indices,mat.indptr]:h.update(arr.tobytes())
    h.update(text_hash(OUT/'gene_universe.tsv').encode())
    dump(cache/'preparation.json',{'fingerprint':h.hexdigest(),'cells':len(d),'observations':len(obs),'groups':len(groups),'genes':genes,'synapse':syn,'background':background,'targets':targets})
    dump(OUT/'preparation.json',{'fingerprint':h.hexdigest(),'cells':len(d),'observations':len(obs),'groups':len(groups),'genes':len(genes),'background':len(background),'comparators':len(comparators)})
    print('Preparation complete:',len(d),'cells;',len(obs),'target-mouse observations',flush=True)

def batch(cache,spec,delete_raw=True):
    prep=json.loads((cache/'preparation.json').read_text());name=Path(spec['file']).stem
    cp=cache/'checkpoints';cp.mkdir(exist_ok=True)
    dest=cp/(name+'.npz');marker=cp/(name+'.json')
    if dest.exists() and marker.exists():
        m=json.loads(marker.read_text())
        if m['source_sha256']==spec['sha256'] and m['fingerprint']==prep['fingerprint'] and m['reduction_code_sha256']==text_hash(__file__) and sha(dest)==m['checkpoint_sha256']:
            print('Reuse',name,flush=True);return
        raise RuntimeError('Checkpoint mismatch: '+str(dest))
    path=download(spec,cache/'raw')
    d=pd.read_parquet(cache/'selected_cells.parquet');A=sparse.load_npz(cache/'A.npz');C=sparse.load_npz(cache/'C.npz')
    genes=prep['genes'];ng=len(genes);targetcols=np.array([genes.index(g) for g in prep['targets']])
    obs=pd.read_csv(OUT/'eligible_observations.tsv',sep='\t');own=np.array([prep['targets'].index(g) for g in obs.target])
    with h5py.File(path,'r') as f:
        ids=read_col(f['obs'],f['obs'].attrs['_index']);assert len(set(ids))==len(ids)
        pos=pd.Index(ids).get_indexer(d.cell_id);selected=np.flatnonzero(pos>=0);order=np.argsort(pos[selected]);selected=selected[order];rawrows=pos[selected]
        var=f['var'];symbols=read_col(var,'gene_name') if 'gene_name' in var else read_col(var,var.attrs['_index'])
        # Gene features must uniquely match the complete biological browser universe.
        mapping={}; gene_set=set(genes)
        for j,s in enumerate(symbols):
            if s in gene_set:
                if s in mapping:raise RuntimeError('Duplicate raw symbol: '+s)
                mapping[s]=j
        assert set(mapping)==set(genes),'Raw/browser gene universe mismatch'
        feature_types=read_col(var,'feature_types')
        assert all(feature_types[mapping[g]]=='Gene Expression' for g in genes),'Non-RNA feature in gene universe'
        gi=np.array([mapping[g] for g in genes]);raw_to_gene=np.full(len(symbols),-1,dtype=int);raw_to_gene[gi]=np.arange(ng)
        umi=read_col(f['obs'],'num_rna_umi')
        assert np.array_equal(umi[rawrows],d.iloc[selected].num_rna_umi.to_numpy()),'RNA denominators changed'
        assert (umi[rawrows]>0).all()
        X=f['X'];assert X.attrs['encoding-type']=='csr_matrix'
        ptr=X['indptr'][:];active=np.unique(d.iloc[selected].group.to_numpy());local=np.full(A.shape[1],-1,dtype=int);local[active]=np.arange(len(active))
        sums=np.zeros((len(active),ng),dtype=np.float64);detect=np.zeros((len(active),len(targetcols)),dtype=np.float64)
        for start in range(0,len(selected),512):
            ix=selected[start:start+512];rr=rawrows[start:start+512]
            values=[];indices=[];ip=[0]
            for r in rr:
                val=X['data'][ptr[r]:ptr[r+1]];col=X['indices'][ptr[r]:ptr[r+1]]
                assert np.isfinite(val).all() and (val>=0).all() and np.equal(val,np.floor(val)).all(),'X is not raw integer counts'
                c=raw_to_gene[col];keep=c>=0
                values.append(np.log1p(10000.*val[keep].astype(float)/float(umi[r])).astype(np.float32));indices.append(c[keep]);ip.append(ip[-1]+int(keep.sum()))
            block=sparse.csr_matrix((np.concatenate(values),np.concatenate(indices),ip),shape=(len(ix),ng))
            groupids=local[d.iloc[ix].group.to_numpy()]
            np.add.at(sums,groupids,block.toarray())
            np.add.at(detect,groupids,(block[:,targetcols].toarray()>0).astype(float))
            if start%10240==0:print(name,'selected cells',start,'/',len(selected),flush=True)
        delta=np.asarray(A[:,active]@sums)
        control=np.asarray(C[:,active]@sums[:,targetcols])[np.arange(len(obs)),own]
        detection=np.asarray(C[:,active]@detect)[np.arange(len(obs)),own]
    tmp=dest.with_name(dest.name+'.tmp')
    with tmp.open('wb') as out:np.savez_compressed(out,delta=delta,control=control,detection=detection,selected=selected)
    os.replace(tmp,dest)
    dump(marker,{'source_sha256':spec['sha256'],'fingerprint':prep['fingerprint'],'checkpoint_sha256':sha(dest),'reduction_code_sha256':text_hash(__file__),'selected_cells':len(selected)})
    if delete_raw:path.unlink()
    print('Verified reduction',name,flush=True)

def finalize(cache):
    prep=json.loads((cache/'preparation.json').read_text());obs=pd.read_csv(OUT/'eligible_observations.tsv',sep='\t')
    specs=json.loads((DOC/'batch_sources.json').read_text());delta=np.zeros((len(obs),len(prep['genes'])))
    control=np.zeros(len(obs));detection=np.zeros(len(obs));seen=np.zeros(prep['cells'],dtype=np.int16)
    sources=[]
    for spec in specs:
        name=Path(spec['file']).stem;dest=cache/'checkpoints'/(name+'.npz');m=json.loads(dest.with_suffix('.json').read_text())
        assert m['source_sha256']==spec['sha256'] and m['fingerprint']==prep['fingerprint'] and m['reduction_code_sha256']==text_hash(__file__) and sha(dest)==m['checkpoint_sha256']
        with np.load(dest) as z:
            delta+=z['delta'];control+=z['control'];detection+=z['detection'];seen[z['selected']]+=1
        sources.append(m)
    assert np.all(seen==1),'Selected cells missing or duplicated across batches'
    gi={g:i for i,g in enumerate(prep['genes'])}
    syn=delta[:,[gi[g] for g in prep['synapse']]].mean(1)
    err=np.abs(syn-obs.difference.to_numpy());assert err.max()<=1e-6,'Original synaptic effects not reproduced'
    obs['reconstructed_synaptic_effect']=syn
    obs['global_rms']=np.sqrt(np.mean(delta[:,[gi[g] for g in prep['background']]]**2,axis=1))
    obs['realized_target_difference']=delta[np.arange(len(obs)),[gi[g] for g in obs.target]]
    obs['control_target_mean']=control;obs['control_target_detection']=detection
    obs['measurable']=(control>=.1)&(detection>=.1);obs['log_target_cells']=np.log(obs.target_cells)
    assert np.isfinite(obs[['global_rms','realized_target_difference','control_target_mean','control_target_detection']].to_numpy()).all()
    assert ((detection>=-1e-10)&(detection<=1+1e-10)).all()
    save(obs,OUT/'global_metrics.tsv');dump(OUT/'source_reconciliation.json',{'status':'PASS','maximum_synaptic_discrepancy':float(err.max()),'all_selected_cells_once':True,'sources':sources,'fingerprint':prep['fingerprint']})

def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','batches','finalize']);p.add_argument('--cache',required=True,type=Path);a=p.parse_args()
    if a.stage=='prepare':prepare(a.cache)
    elif a.stage=='finalize':finalize(a.cache)
    else:
        from global_disruption_validate import first_batch_check
        specs=json.loads((DOC/'batch_sources.json').read_text())
        prep=json.loads((a.cache/'preparation.json').read_text());marker=a.cache/'first_validation.json'
        already=False
        if marker.exists():
            check=json.loads(marker.read_text());already=check['fingerprint']==prep['fingerprint'] and check['status']=='PASS'
        if not already:
            first_batch_check(a.cache)
            batch(a.cache,specs[0],delete_raw=False)
            check=first_batch_check(a.cache);check['fingerprint']=prep['fingerprint'];dump(marker,check)
            (a.cache/'raw'/specs[0]['file']).unlink()
        else:
            batch(a.cache,specs[0]);dump(OUT/'first_batch_validation.json',check)
        for spec in specs[1:]:batch(a.cache,spec)

if __name__=='__main__':main()
