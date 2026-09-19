"""Independent source-scale check; stops execution on browser/raw disagreement."""
import gzip
import hashlib
import json
import struct
import zlib
from pathlib import Path
from urllib.request import Request
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from global_disruption_metrics import ROOT,DOC,OLD,OUT,read_col
from global_disruption_io import dump,download,sha,verified_urlopen

def first_batch_check(cache):
    spec=json.loads((DOC/'batch_sources.json').read_text())[0]
    raw=download(spec,cache/'raw')
    index=json.loads((cache/'exprMatrix.json').read_text());offset,length=index['Xkr4']
    vectors=cache/'vectors';vectors.mkdir(exist_ok=True);path=vectors/'Xkr4.bin'
    expected=next(x['sha256'] for x in json.loads((OLD/'response_vector_checksums.json').read_text()) if x['gene']=='Xkr4')
    if not path.exists():
        req=Request('https://cells.ucsc.edu/whole-brain-perturb/combined/exprMatrix.bin?Xkr4',headers={'Range':f'bytes={offset}-{offset+length-1}'})
        with verified_urlopen(req,timeout=90) as r:
            assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {offset}-{offset+length-1}/')
            data=r.read(length+1);assert len(data)==length
        assert hashlib.sha256(data).hexdigest()==expected
        path.write_bytes(data)
    assert sha(path)==expected
    data=zlib.decompress(path.read_bytes());n=struct.unpack('<H',data[:2])[0]
    assert 'Xkr4' in data[2:2+n].decode()
    x=np.frombuffer(data,offset=2+n,dtype='<f4');assert len(x)==6348631
    d=pd.read_parquet(cache/'selected_cells.parquet');A=sparse.load_npz(cache/'A.npz')
    with h5py.File(raw) as f:
        ids=read_col(f['obs'],f['obs'].attrs['_index']);pos=pd.Index(ids).get_indexer(d.cell_id)
        selected=np.flatnonzero(pos>=0);rows=pos[selected]
        syms=read_col(f['var'],f['var'].attrs['_index']);col=np.flatnonzero(syms=='Xkr4');assert len(col)==1
        umi=read_col(f['obs'],'num_rna_umi');ptr=f['X/indptr'][:];vals=[]
        for r in rows:
            ix=f['X/indices'][ptr[r]:ptr[r+1]];v=f['X/data'][ptr[r]:ptr[r+1]]
            count=v[ix==col[0]].sum();vals.append(np.float32(np.log1p(10000*float(count)/float(umi[r]))))
    actual=x[d.iloc[selected].position.to_numpy()]
    err=float(np.max(np.abs(np.asarray(vals)-actual)))
    assert err<1e-5,'Original browser expression scale disagrees with raw counts'
    groups=np.bincount(d.iloc[selected].group.to_numpy(),weights=actual,minlength=A.shape[1])
    expected_delta=A@groups
    # Compare to the saved batch reduction if it already exists.
    cp=cache/'checkpoints'/Path(spec['file']).with_suffix('.npz').name
    aggregation_error=None
    if cp.exists():
        prep=json.loads((cache/'preparation.json').read_text())
        with np.load(cp) as z:delta=z['delta'][:,prep['genes'].index('Xkr4')]
        aggregation_error=float(np.max(np.abs(delta-expected_delta)));assert aggregation_error<1e-6
    result={'status':'PASS','selected_cells':len(selected),'maximum_raw_browser_error':err,'maximum_batch_aggregation_error':aggregation_error,'raw_sha256':spec['sha256'],'browser_vector_sha256':expected}
    dump(OUT/'first_batch_validation.json',result)
    return result
