#!/usr/bin/env python3
"""Windows/Linux entry point; explicit status, checkpoint resume, small result export."""
import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
import zipfile
from global_disruption_io import dump,sha

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/'docs/v2/global_disruption_20260919'
OUT=ROOT/'results/v2/global_disruption_20260919'

def artifact_hash(path):
    data=Path(path).read_bytes()
    if Path(path).suffix in {'.py','.md','.txt','.json','.tsv'}:data=data.replace(b'\r\n',b'\n')
    return hashlib.sha256(data).hexdigest()

def verify():
    if sys.version_info[:2]!=(3,12):raise RuntimeError('Use Python 3.12 for this pinned runner.')
    if sys.flags.optimize:raise RuntimeError('Do not run with python -O: validation assertions are required.')
    for line in (DOC/'requirements.txt').read_text().splitlines():
        if not line.strip():continue
        name,expected=line.split('==');actual=importlib.metadata.version(name)
        if actual!=expected:raise RuntimeError(f'{name} {actual}; required {expected}. Install the requirements file.')
    locks=json.loads((DOC/'runner_checksums.json').read_text())
    for rel,digest in locks.items():
        path=ROOT/rel
        if not path.is_file() or artifact_hash(path)!=digest:raise RuntimeError('Frozen runner/input checksum mismatch: '+rel)
    print('Pinned environment and frozen artifacts verified.',flush=True)

def child(args,cache):
    env=os.environ.copy()
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
    env['PYTHONUNBUFFERED']='1'
    cmd=[sys.executable,*args]
    print('Stage:', ' '.join(args),flush=True)
    with (cache/'execution.log').open('a',encoding='utf-8') as log:
        log.write('\nCOMMAND '+repr(cmd)+'\n');log.flush()
        proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
        try:
            for line in proc.stdout:
                print(line,end='',flush=True);log.write(line);log.flush()
            rc=proc.wait()
        except BaseException:
            proc.terminate();proc.wait();raise
        if rc:raise RuntimeError(f'Stage failed with exit {rc}: {args[0]}')

def export(path,cache):
    path.parent.mkdir(parents=True,exist_ok=True)
    files=list(OUT.rglob('*'))+list(DOC.rglob('*'))+list((ROOT/'scripts').glob('*global_disruption*.py'))
    files.extend([ROOT/'docs/v2/specificity_followup_20260918/PLAN.md',ROOT/'docs/v2/reconstruction_20260918/locked_comparators.tsv'])
    tmp=path.with_name(path.name+'.tmp')
    with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for f in sorted(set(files)):
            if f.is_file():z.write(f,str(f.relative_to(ROOT)).replace('\\','/'))
        log=cache/'execution.log'
        if log.exists():z.write(log,'execution.log')
        cp=cache/'checkpoints'
        if cp.exists():
            for f in cp.glob('*.json'):z.write(f,'checkpoint_manifests/'+f.name)
    os.replace(tmp,path)
    print('Results/diagnostics ZIP:',path,flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--export',type=Path,required=True);p.add_argument('--verify-only',action='store_true');a=p.parse_args()
    if a.verify_only:
        verify();return
    cache=a.cache.resolve();destination=a.export.resolve()
    if cache==ROOT or ROOT in cache.parents:raise RuntimeError('Use a dedicated cache outside the project repository.')
    if cache==destination or cache in destination.parents:raise RuntimeError('Export ZIP must be outside the cache.')
    cache.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    sentinel=cache/'runner_cache.json'
    if sentinel.exists():
        if json.loads(sentinel.read_text()).get('purpose')!='pd-envtox item1 rebuilt 20260919':raise RuntimeError('Unrecognized cache directory')
    else:
        if any(cache.iterdir()):raise RuntimeError('Choose a new empty cache folder for the rebuilt runner.')
        dump(sentinel,{'purpose':'pd-envtox item1 rebuilt 20260919','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    status={'status':'INCOMPLETE','stage':'verification','scientific_results':'Pending'}
    dump(OUT/'execution_status.json',status)
    try:
        verify()
        for stage in ['prepare','batches','finalize']:
            status['stage']=stage;dump(OUT/'execution_status.json',status)
            child(['scripts/global_disruption_metrics.py',stage,'--cache',str(cache)],cache)
        status['stage']='model';dump(OUT/'execution_status.json',status)
        child(['scripts/global_disruption_fit.py'],cache)
        status['stage']='report';dump(OUT/'execution_status.json',status)
        child(['scripts/global_disruption_report.py'],cache)
        model=json.loads((OUT/'model_status.json').read_text())
        status.update(status='COMPLETE',stage='finished',scientific_results=model['status'],note='Computational execution complete; scientific audit still required. NOT_EVALUABLE is a valid preserved outcome.')
    except BaseException as error:
        status['error']=str(error);status['error_type']=type(error).__name__
        with (cache/'execution.log').open('a',encoding='utf-8') as f:traceback.print_exc(file=f)
        raise
    finally:
        dump(OUT/'execution_status.json',status)
        export(destination,cache)

if __name__=='__main__':main()
