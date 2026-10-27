from pathlib import Path
import json,anndata as ad,numpy as np,pandas as pd
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
records=[]
for p in sorted((ROOT/'03_processed').glob('GSE205013_P*_identity_final.json')):
    sample=p.name.split('_')[1];a=ad.read_h5ad(ROOT/'03_processed'/f'GSE205013_{sample}_annotated.h5ad')
    counts=a.layers['counts'];assert np.isfinite(counts.data).all() and (counts.data>=0).all() and (counts.data==np.rint(counts.data)).all()
    assert counts.data.max()<2**24
    mal=(a.obs.malignancy=='Malignant_CNV').to_numpy();epi=a.obs.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy();strict=a.obs.qc_strict.to_numpy()
    masks={'Malignant':mal,'EpiAll':epi,'MalStrictQC':mal&strict,'EpiStrictQC':epi&strict};masks.update({'Type_'+t:(a.obs.major_type==t).to_numpy() for t in a.obs.major_type.unique()})
    pbpath=ROOT/'03_processed'/f'GSE205013_{sample}_pseudobulk.npz';pb=np.load(pbpath);rebuilt={};changed=False
    assert np.array_equal(pb['genes'].astype(str),a.var_names.to_numpy(dtype=str))
    nums=json.loads((ROOT/'03_processed'/f'GSE205013_{sample}_cellnumbers.json').read_text())
    for group,mask in masks.items():
        expected=np.asarray(counts[mask].sum(0,dtype=np.int64)).ravel();rebuilt[group]=expected
        assert int(mask.sum())==nums[group]
        equal=np.array_equal(expected,pb[group]);changed|=not equal
        records.append(dict(dataset='GSE205013',patient_id=sample,cell_group=group,n_cells=int(mask.sum()),sum_counts=int(expected.sum()),max_gene_count=int(expected.max()),exact_sum_match_before=equal,max_rounding_difference=int(np.max(np.abs(expected-pb[group]))),status='PASS' if equal else 'REPAIRED_INT64_ACCUMULATION'))
    if changed:np.savez_compressed(pbpath,genes=a.var_names.to_numpy(dtype=str),**rebuilt)
pd.DataFrame(records).to_csv(ROOT/'04_qc/pseudobulk_count_integrity.tsv',sep='\t',index=False)
assert len({r['patient_id'] for r in records})==27
print('COUNT_INTEGRITY_COMPLETE',len(records),'rows','repairs',sum(r['status']!='PASS' for r in records))
