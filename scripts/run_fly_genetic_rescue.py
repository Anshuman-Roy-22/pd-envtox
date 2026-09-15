#!/usr/bin/env python3
"""Reproduce the standalone fly extension; earlier V2/R stages are separate."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch', action='store_true', help='Acquire missing checksum-locked inputs')
    args = parser.parse_args()
    command = [sys.executable, str(ROOT / 'scripts/fly_genetic_rescue.py')]
    if args.fetch:
        command.append('--fetch')
    subprocess.run(command, cwd=ROOT, check=True)
    for script in ['fly_genetic_rescue_figure.py', 'verify_fly_genetic_rescue.py']:
        subprocess.run([sys.executable, str(ROOT / 'scripts' / script)], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
