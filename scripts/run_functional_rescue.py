#!/usr/bin/env python3
"""Acquire, reproduce, plot and verify the bounded functional-rescue extension."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch',action='store_true',help='Download missing checksum-locked author files')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    jobs=[['acquire_functional_rescue_metadata.py']+(['--fetch'] if args.fetch else []),
          ['functional_rescue_design_audit.py'],['functional_rescue_analysis.py'],
          ['functional_rescue_figures.py'],['verify_functional_rescue.py']]
    for job in jobs:
        print('Running '+job[0],flush=True)
        subprocess.run([sys.executable,str(root/'scripts'/job[0]),*job[1:]],cwd=root,check=True)


if __name__=='__main__':
    main()
