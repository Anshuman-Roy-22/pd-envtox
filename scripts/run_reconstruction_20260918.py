#!/usr/bin/env python3
"""Verify frozen inputs, reproduce this addition, and verify its scientific outputs."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/'docs/v2/reconstruction_20260918'
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while b:=f.read(4*1024*1024):h.update(b)
    return h.hexdigest()
def verify(records):
    for path,expected in records.items():
        if digest(ROOT/path)!=expected:raise ValueError('Frozen artifact differs: '+path)
def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--r',default='R',help='R executable, version 4.3.3; not Rscript');p.add_argument('--offline',action='store_true');p.add_argument('--verify-only',action='store_true');a=p.parse_args()
    lock=json.loads((DOC/'artifact_checksums.json').read_text());verify(lock['inputs_and_code'])
    if not a.verify_only:
        cache=str(a.cache.resolve())
        commands=[['fetch_reconstruction_20260918.py','--cache',cache]+(['--offline'] if a.offline else []),['reconstruct_human_20260918.py']]
        for cmd in commands:subprocess.run([sys.executable,str(ROOT/'scripts'/cmd[0]),*cmd[1:]],cwd=ROOT,check=True)
        subprocess.run([a.r,'--vanilla','--slave','-f',str(ROOT/'scripts/reconstruct_human_20260918.R'),'--args',str(ROOT),str(ROOT/'results/v2/reconstruction_20260918/human')],cwd=ROOT,check=True)
        for cmd in [ ['reconstruct_mouse_20260918.py','metadata','--cache',cache],['reconstruct_mouse_20260918.py','baseline','--cache',cache] ]:
            subprocess.run([sys.executable,str(ROOT/'scripts'/cmd[0]),*cmd[1:]],cwd=ROOT,check=True)
        # Reproduction may regenerate selection, but it must equal the pre-outcome lock.
        verify(lock['inputs_and_code'])
        for name in ['reconstruct_mouse_fit_20260918.py','reconstruction_20260918_support.py']:
            subprocess.run([sys.executable,str(ROOT/'scripts'/name),'--cache',cache],cwd=ROOT,check=True)
    verify(lock['inputs_and_code']);verify(lock['scientific_outputs'])
    print('PASS: frozen selection and all scientific artifact checksums match.',flush=True)
if __name__=='__main__':main()
