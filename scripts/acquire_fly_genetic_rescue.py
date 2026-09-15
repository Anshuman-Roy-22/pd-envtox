#!/usr/bin/env python3
"""Acquire checksum-locked workbooks; inventory labels without effect estimates."""
import argparse
import hashlib
import json
import os
from pathlib import Path
from urllib.request import Request, urlopen

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs/v2/fly_genetic_rescue'
SHEETS = ['Fig3d_GI_cluster', 'Fig4b_Q10', 'Fig4d_R55',
          'Suppl Fig12b_Q10_SING', 'Suppl Fig12c_R55_SING']


def acquire(cache, fetch=False):
    manifest = json.loads((DOCS / 'input_manifest.json').read_text())
    inventory = {'scope': 'headers, dimensions, nonempty counts; no effect estimates', 'files': {}}
    for item in manifest['files']:
        name = Path(item['path'])
        if name.is_absolute() or '..' in name.parts:
            raise ValueError('Unsafe manifest path')
        path = cache / name
        if not path.exists():
            if not fetch:
                raise FileNotFoundError(f'{path}: run with --fetch')
            request = Request(item['url'], headers={'User-Agent': 'pd-envtox-source-audit'})
            with urlopen(request, timeout=60) as response:
                raw = response.read()
            if hashlib.sha256(raw).hexdigest() != item['sha256']:
                raise ValueError(f'Unexpected source checksum: {name}')
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(path.name + '.partial')
            temporary.write_bytes(raw)
            os.replace(temporary, path)
        raw = path.read_bytes()
        if len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']:
            raise ValueError(f'Cached source differs: {name}')
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        sheets = SHEETS if name.name == 'source_data.xlsx' else wb.sheetnames
        record = {'bytes': len(raw), 'sha256': item['sha256'], 'sheets': {}}
        for title in sheets:
            rows = list(wb[title].values)
            record['sheets'][title] = {
                'rows_including_header': len(rows), 'columns': len(rows[0]),
                'header': list(rows[0]),
                'nonempty_per_column': [sum(row[j] is not None for row in rows[1:])
                                        for j in range(len(rows[0]))]}
        wb.close()
        inventory['files'][str(name)] = record
    return inventory


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache', type=Path, default=ROOT / 'data_raw/fly_genetic_rescue')
    p.add_argument('--output', type=Path, default=DOCS / 'metadata_inventory.json')
    p.add_argument('--fetch', action='store_true')
    args = p.parse_args()
    inventory = acquire(args.cache, args.fetch)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(inventory, indent=2, ensure_ascii=False) + '\n')
    print('Verified two source workbooks and wrote metadata inventory.')


if __name__ == '__main__':
    main()
