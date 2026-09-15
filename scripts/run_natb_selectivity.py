#!/usr/bin/env python3
"""Reproduce and verify the complete NatB extension without running older analyses."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def verify_manifest(path):
    entries = path.read_text().splitlines()
    for line in entries:
        digest, relative = line.split('  ', 1)
        source = ROOT / relative
        if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != digest:
            raise RuntimeError('Checksum mismatch: ' + relative)
    return len(entries)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fetch', action='store_true', help='Acquire absent source files and check frozen inner-file hashes')
    p.add_argument('--verify-only', action='store_true', help='Check tracked dependencies, code, result bytes and labels')
    a = p.parse_args()
    inputs = verify_manifest(ROOT / 'docs/v2/natb_selectivity/SHA256SUMS.txt')
    if not a.verify_only:
        for script in ['natb_selectivity.py', 'natb_robustness.py', 'plot_natb_selectivity.py']:
            cmd = [sys.executable, str(ROOT / 'scripts' / script)]
            if a.fetch and script == 'natb_selectivity.py': cmd.append('--fetch')
            subprocess.run(cmd, cwd=ROOT, check=True)
    outputs = verify_manifest(ROOT / 'results/v2/natb_selectivity/SHA256SUMS.txt')
    verdict = json.loads((ROOT / 'results/v2/natb_selectivity/validation_summary.json').read_text())
    if verdict['classification'] != 'SELECTIVITY_NOT_ESTABLISHED':
        raise RuntimeError('Frozen result label changed; audit before making a new release')
    print(json.dumps({'release_verification': 'PASS', 'input_code_files': inputs, 'result_files': outputs,
                      'primary_classification': verdict['classification'],
                      'scope': 'NatB extension only; does not close whole-project R execution gate'}, indent=2))


if __name__ == '__main__':
    main()
