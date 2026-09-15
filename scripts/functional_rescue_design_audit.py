#!/usr/bin/env python3
"""Audit compound and well identities, without calculating phenotypic effects."""
import hashlib
import json
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs/v2/functional_rescue'
CACHE = ROOT / 'data_raw/functional_rescue'
MAIN = 'Fig2/Notebooks/1_Primary_Screen_Data_Classification/'


def well(value):
    match = re.fullmatch(r'([A-P])(\d{1,2})', value.strip())
    if not match:
        raise ValueError(f'Invalid well identifier: {value}')
    return f'{match[1]}{int(match[2]):02d}'


def load_design():
    manifest = json.loads((DOCS / 'input_manifest.json').read_text())
    for item in manifest['files']:
        if hashlib.sha256((CACHE / item['path']).read_bytes()).hexdigest() != item['sha256']:
            raise ValueError(f"Input changed: {item['path']}")
    a = pd.read_csv(CACHE / (MAIN + 'all_raw_data.csv'), usecols=['Plate','Well','tags'], dtype=str)
    m = pd.read_csv(CACHE / (MAIN + 'data_for_compound_mapping.csv'), dtype=str)
    m = m.rename(columns={'Destination Plate Barcode':'Plate', 'Destination Well':'Well',
                          'Concentration (M)':'dose_M', 'CpdID':'compound_name'})
    d = pd.read_csv(CACHE / 'Fig3/Notebooks/Exp72_DRC.csv', dtype=str,
                    usecols=['Plate','Well','tags','SC_cat','SC_name','final_conc'])
    d = d.rename(columns={'SC_name':'compound_name','final_conc':'dose_M'})
    for df in (a, m, d):
        df['Well'] = df['Well'].map(well)
        if df.duplicated(['Plate','Well']).any():
            raise ValueError('Duplicate plate/well keys')
    a = a.merge(m, on=['Plate','Well'], how='left', validate='one_to_one')
    if a.loc[a.tags=='Compound','SC_cat'].isna().any() or a.loc[a.tags!='Compound','SC_cat'].notna().any():
        raise ValueError('Compound map does not match raw well labels')
    for df in (a,d):
        df['compound_name'] = df['compound_name'].str.strip()
        df['dose_M'] = pd.to_numeric(df['dose_M'])
    return a,d


def main():
    a,d = load_design()
    catalogue = a[a.tags=='Compound'][['SC_cat','compound_name','dose_M']].drop_duplicates()
    if catalogue.SC_cat.duplicated().any():
        raise ValueError('Multiple primary identities/doses for catalogue ID')
    rows = []
    for cat, g in d[d.SC_cat.notna()].groupby('SC_cat', sort=True):
        source = catalogue[catalogue.SC_cat==cat]
        item = {'SC_cat':cat,'followup_name':g.compound_name.iloc[0],
                'primary_name':None,'primary_dose_M':None,'followup_dose_M':None,
                'relative_dose_difference':None,'status':'UNMATCHED_CATALOGUE',
                'primary_wells':0,'followup_wells':0}
        if len(source)==1:
            dose = float(source.dose_M.iloc[0])
            ds = sorted(g.dose_M.dropna().unique(), key=lambda v: (abs(v-dose),v))
            matched = ds[0]
            gap = abs(matched-dose)/dose
            item.update(primary_name=source.compound_name.iloc[0],primary_dose_M=dose,
                        followup_dose_M=matched,relative_dose_difference=gap,
                        primary_wells=int((a.SC_cat==cat).sum()),
                        followup_wells=int((g.dose_M==matched).sum()),
                        status='MATCHED' if gap<=0.05 else 'DOSE_MISMATCH')
        rows.append(item)
    table = pd.DataFrame(rows)
    table.to_csv(DOCS / 'matched_compound_design.tsv', sep='\t', index=False, na_rep='NA')
    audit = {
        'scope':'Design labels and input hashes only; no drug effects calculated',
        'primary_wells':len(a),'primary_compounds':int(a.SC_cat.nunique()),
        'primary_plates':sorted(a.Plate.unique()),'followup_wells':len(d),
        'followup_compounds':int(d.SC_cat.nunique()),'followup_plates':sorted(d.Plate.unique()),
        'matched_compounds':int((table.status=='MATCHED').sum()),
        'unmatched':table.loc[table.status!='MATCHED',['SC_cat','followup_name','status']].to_dict('records'),
        'primary_control_wells':a[a.tags!='Compound'].groupby(['Plate','tags']).size().reset_index(name='n').to_dict('records'),
        'followup_control_wells':d[d.SC_cat.isna()].groupby(['Plate','tags']).size().reset_index(name='n').to_dict('records'),
        'biological_replicate_policy':'Plate groups are recorded; biological/donor independence is not inferred from plate names or biochemical rows.',
        'source_exclusions':'The author DRC notebook removes Alectinib explicitly. This analysis retains deposited Alectinib and matches its actual 2.5 uM primary dose; no supplied reason establishes its invalidity.',
    }
    (DOCS / 'design_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(table[['SC_cat','followup_name','status']].to_string(index=False))


if __name__=='__main__':
    main()
