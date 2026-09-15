#!/usr/bin/env python3
"""Verify source/output integrity and scientifically consequential edge cases."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from functional_rescue_analysis import CONFIG, ROOT, VALUES, average_precision, auc, normalize
from functional_rescue_design_audit import CACHE, DOCS, load_design


def main():
    out=ROOT/'results/v2/functional_rescue'
    manifest=json.loads((DOCS/'input_manifest.json').read_text())
    load_design()  # Source SHA256, well uniqueness, and exact catalogue mapping.
    verified=[]
    for line in (out/'SHA256SUMS.txt').read_text().splitlines():
        expected,relative=line.split('  ',1)
        path=ROOT/relative
        assert path.is_file(), relative
        assert hashlib.sha256(path.read_bytes()).hexdigest()==expected, relative
        verified.append(relative)
    parsed=[]
    for path in sorted(out.iterdir()):
        if path.name.endswith('.tsv') or path.name.endswith('.tsv.gz'):
            if path.suffix=='.gz':gzip.decompress(path.read_bytes())
            table=pd.read_csv(path,sep='\t')
            assert not table.columns.duplicated().any(),path.name
            parsed.append({'file':path.name,'rows':len(table)})
    # A hand-derived tied-score example: AP = 1/2*1/2 + 1/2*2/3.
    y=[1,0,1,0]; s=[.9,.9,.2,.1]
    assert np.isclose(average_precision(y,s),7/12,atol=1e-15)
    assert np.isclose(auc(y,s),.625,atol=1e-15)
    assert np.isnan(average_precision([0,0],[.9,.1]))
    assert np.isnan(auc([0,0],[.9,.1]))
    assert np.isnan(auc([1,1],[.9,.1]))
    raw=pd.read_csv(CACHE/'Fig1/Notebooks/Exp67_4plates_controls.csv')
    base,_=normalize(raw,'Mut - DMSO','WT - DMSO')
    altered=raw.copy()
    for i,field in enumerate(VALUES):
        altered[field]=altered[field]*(i+2)+(i*3 if i!=4 else 0)
    affine,_=normalize(altered,'Mut - DMSO','WT - DMSO')
    np.testing.assert_allclose(base[[f'z{i}' for i in range(6)]],affine[[f'z{i}' for i in range(6)]],rtol=1e-11,atol=1e-11)
    np.testing.assert_allclose(base.ratio4,affine.ratio4,rtol=1e-12,atol=1e-12)
    # Check the critical count result directly against original deposited counts.
    raw_d=pd.read_csv(CACHE/'Fig3/Notebooks/Exp72_DRC.csv')
    m=pd.read_csv(out/'matched_benchmark.tsv',sep='\t')
    for _,row in m[m.followup_alpha_response].iterrows():
        g=raw_d[(raw_d.SC_cat==row.SC_cat)&np.isclose(raw_d.final_conc,row.followup_dose_M,rtol=1e-8,atol=0)]
        ratios=[]
        for _,well in g.iterrows():
            vehicle=raw_d[(raw_d.Plate==well.Plate)&(raw_d.tags=='Mut;DMSO')].Nuclei_Number_Living
            ratios.append(well.Nuclei_Number_Living/vehicle.median())
        assert len(ratios)==4
        assert np.isclose(np.median(ratios),row.followup_Nuclei_Number_Living_ratio,rtol=1e-10)
    summary=json.loads((out/'analysis_summary.json').read_text())
    assert summary['primary']['decision']=='BENCHMARK_NOT_EVALUABLE'
    assert summary['primary']['positive']==0 and summary['primary']['n_compounds']==31
    for key in ['alpha_AP','multivariate_AP','AP_gain','alpha_AUC','multivariate_AUC',
                'leave_one_compound_out_positive_gain_fraction']:
        assert summary['primary'][key] is None,key
    sensitivity=pd.read_csv(out/'sensitivity_summary.tsv',sep='\t')
    assert len(sensitivity)==7 and (sensitivity.positive==0).all()
    ctrl=pd.read_csv(out/'control_transfer_diagnostics.tsv',sep='\t')
    wt=ctrl[(ctrl.dataset=='followup')&(ctrl.control_group=='Corrected vehicle')].iloc[0]
    assert wt.n_wells==48 and wt['passing_joint_0.8']==0
    info={'status':'PASS','source_files_verified':len(manifest['files']),
          'release_files_verified':len(verified),'parsed_tables':parsed,
          'checks':['source checksums and one-to-one well/catalogue design',
                    'release checksums and gzip integrity',
                    'tied-score AP/AUC and undefined-class metrics',
                    'positive-affine normalization invariance',
                    'marker-responder cell-count ratios from original source',
                    'fixed primary and seven sensitivity outcomes',
                    'failed corrected-control transfer reported'],
          'scope':'This functional-rescue extension only; does not verify the repository-wide V2/R pipeline.'}
    (out/'verification.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps({k:v for k,v in info.items() if k not in ['parsed_tables','checks']},indent=2))


if __name__=='__main__':
    main()
