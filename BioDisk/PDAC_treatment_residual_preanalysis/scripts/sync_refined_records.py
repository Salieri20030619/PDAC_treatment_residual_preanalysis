from pathlib import Path
import h5py,json
from anndata.io import read_elem,write_elem
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
for marker in (ROOT/'03_processed').glob('GSE205013_P*_cnv_refined.json'):
    sample=marker.name.split('_')[1];base=ROOT/'03_processed'
    with h5py.File(base/f'GSE205013_{sample}_annotated.h5ad','r') as f:obs=read_elem(f['obs'])
    cnpath=base/f'GSE205013_{sample}_cnv.h5ad'
    if cnpath.exists():
        with h5py.File(cnpath,'r+') as f:
            cnobs=read_elem(f['obs'])
            for key in ['cnv_nchr','cnv_positive']:
                if key+'_initial' not in cnobs:cnobs[key+'_initial']=cnobs[key]
            for key in ['cnv_nchr','cnv_positive','cnv_coherent','malignancy']:
                if key in obs:cnobs[key]=obs.loc[cnobs.index,key]
            write_elem(f,'obs',cnobs)
    done=base/f'GSE205013_{sample}.complete.json';rec=json.loads(done.read_text())
    rec['n_malignant_initial']=rec.get('n_malignant_initial',rec['n_malignant']);rec['n_malignant']=int((obs.malignancy=='Malignant_CNV').sum());rec['cnv_refined']=True;done.write_text(json.dumps(rec,indent=2))
    print('SYNCED',sample,rec['n_malignant'])
print('REFINED_RECORD_SYNC_COMPLETE')
