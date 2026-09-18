#!/usr/bin/env python3
"""Acquire/verify the frozen reconstruction inputs without regenerating old analyses."""
import argparse, hashlib, json, zipfile
from pathlib import Path
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor
from natb_selectivity import acquire
from reconstruct_mouse_20260918 import ROOT, DOC, vector

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while block:=f.read(4*1024*1024):h.update(block)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--offline',action='store_true');a=p.parse_args();a.cache.mkdir(parents=True,exist_ok=True)
    acquire(ROOT/'data_raw/natb_selectivity',not a.offline)
    def get(item):
        dest=a.cache/item['file']
        if not dest.exists():
            if a.offline:raise FileNotFoundError(dest)
            tmp=dest.with_suffix(dest.suffix+'.part')
            with urlopen(item['url'],timeout=120) as src,tmp.open('wb') as out:
                while block:=src.read(4*1024*1024):out.write(block)
            if tmp.stat().st_size!=item['bytes'] or sha(tmp)!=item['sha256']:raise ValueError('Source changed: '+str(dest))
            tmp.replace(dest)
        if dest.stat().st_size!=item['bytes'] or sha(dest)!=item['sha256']:raise ValueError('Source changed: '+str(dest))
        return item['file']
    records=json.loads((DOC/'public_source_manifest.json').read_text())
    for name in ThreadPoolExecutor(4).map(get,records):print('Verified',name,flush=True)
    with zipfile.ZipFile(a.cache/'wholebrain_supp.zip') as z:
        for name in ['media-1.xlsx','media-3.xlsx']:(a.cache/name).write_bytes(z.read(name))
    records={}
    for name in ['baseline_vector_checksums.json','response_vector_checksums.json']:
        for row in json.loads((DOC/name).read_text()):
            if row['gene'] in records:assert row==records[row['gene']]
            records[row['gene']]=row
    def check(row):
        dest=a.cache/'vectors'/(row['gene']+'.bin')
        if a.offline and not dest.exists():raise FileNotFoundError(dest)
        vector(a.cache,row['gene'])
        if dest.stat().st_size!=row['bytes'] or sha(dest)!=row['sha256']:raise ValueError('Expression vector changed: '+row['gene'])
    list(ThreadPoolExecutor(4).map(check,records.values()))
    print('Verified',len(records),'expression vectors',flush=True)

if __name__=='__main__':main()
