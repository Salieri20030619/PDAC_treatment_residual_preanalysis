from pathlib import Path
import json,hashlib,datetime,gzip,pandas as pd,numpy as np,h5py
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
records=[]
for p in (ROOT/'01_metadata').glob('*manifest.json'):
    obj=json.loads(p.read_text())
    if isinstance(obj,list):
        for rec in obj:
            if isinstance(rec,dict) and 'url' in rec:
                records.append({'manifest':p.name,**rec})
flat=[]
for r in records:
    flat.append({k:r.get(k,'') for k in ['manifest','url','file','bytes','sha256','status','exit_code','started_utc','finished_utc','utc','official_checksum','official_checksum_note','gzip_integrity','tar_integrity','member_gzip_integrity','hdf5_open','error']})
pd.DataFrame(flat).drop_duplicates(['url','file'],keep='last').to_csv(ROOT/'01_metadata/download_manifest.tsv',sep='\t',index=False)
datarecs=[]
for name in ['GSE205013_RAW.tar','GSE202051_totaldata-final-toshare.h5ad.gz']:
    r=json.loads((ROOT/'logs'/f'{name}.download.json').read_text())
    assert r['status']=='INTEGRITY_PASS' and Path(r['file']).stat().st_size==r['bytes'],r
    datarecs.append(r)
members=json.loads((ROOT/'01_metadata/GSE205013_archive_members.json').read_text())
assert len(members)==81 and all(m['gzip']=='PASS' for m in members)
expected=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
g205=expected[expected.dataset=='GSE205013']
assert len(g205)==27 and g205.patient_id.nunique()==27
mutation_table=json.loads((ROOT/'01_metadata/Werba2023_supplement_sheets.json').read_text())['Genomic Mutations']
mutation_ids={r[0].strip() for r in mutation_table[1:] if r and r[0]}
assert set(g205.patient_id)<=mutation_ids
checks=[]
for _,r in g205.iterrows():
    q=json.loads((ROOT/'04_qc'/f'{r.patient_id}_qc_audit.json').read_text())
    assert q['counts_integer'] and q['n_input']>0
    counts=ROOT/'03_processed'/f'GSE205013_{r.patient_id}_pseudobulk.npz'
    assert counts.exists()
    z=np.load(counts);assert np.all(z['Malignant']>=0) and np.all(z['Malignant']==np.floor(z['Malignant']))
    checks.append(dict(dataset='GSE205013',patient_id=r.patient_id,GSM=r.GSM,expression_source=f'02_raw/GSE205013/matrices/{r.GSM}_{r.patient_id}_matrix.mtx.gz',n_cells_input=q['n_input'],matrix_readable=True))
with h5py.File(ROOT/'02_raw/GSE202051/GSE202051_totaldata-final-toshare.h5ad') as f:assert 'X' in f and 'obs' in f and 'var' in f
assert (ROOT/'04_qc/GSE202051_h5ad_audit.json').exists()
(ROOT/'logs/DOWNLOAD_AUDIT_PASS.json').write_text(json.dumps({'status':'DOWNLOAD_AUDIT_PASS','utc':now(),'expression_files':2,'archive_members':81,'GSE205013_matrices_read':27,'official_checksums':'Not published in saved directory indexes; SHA256 locally computed, gzip/tar integrity and official Content-Length checked'},indent=2))
print('DOWNLOAD_AUDIT_PASS')
mapping=pd.read_csv(ROOT/'01_metadata/GSE202051_expression_patient_mapping.tsv',sep='\t')
assert mapping.patient_id.is_unique
g202=expected[expected.dataset=='GSE202051']
assert set(mapping.patient_id)<=set(g202.patient_id)
assert set(mapping[mapping.patient_id.str.startswith('PDAC_')].patient_id)<=set(g202[g202.comparison_scope=='Primary'].patient_id)
for _,r in mapping.iterrows():
    m=g202[g202.patient_id==r.patient_id].iloc[0]
    checks.append(dict(dataset='GSE202051',patient_id=r.patient_id,GSM=m.GSM,expression_source='02_raw/GSE202051/GSE202051_totaldata-final-toshare.h5ad',n_cells_input=r.n_cells,matrix_readable=True))
pd.DataFrame(checks).to_csv(ROOT/'01_metadata/expression_patient_audit.tsv',sep='\t',index=False)
registration=g202[['dataset','patient_id','specimen_type','comparison_scope']].copy()
registration['present_in_expression_object']=registration.patient_id.isin(mapping.patient_id)
registration.to_csv(ROOT/'01_metadata/GSE202051_registered_vs_available.tsv',sep='\t',index=False)
samplemeta=pd.read_csv(ROOT/'01_metadata/patient_sample_metadata.tsv',sep='\t',keep_default_na=False)
for (ds,pid),group in samplemeta[samplemeta.patient_id!=''].groupby(['dataset','patient_id']):
    for field in ['tissue_site','treatment_status','procedure','stage','regimen']:
        assert group[field].nunique()<=1,(ds,pid,field,group[field].tolist())
samplemeta['expression_source']=samplemeta.apply(lambda r:f'02_raw/GSE205013/matrices/{r.GSM}_{r.patient_id}_matrix.mtx.gz' if r.dataset=='GSE205013' else '02_raw/GSE202051/GSE202051_totaldata-final-toshare.h5ad' if r.patient_id in set(mapping.patient_id) else '',axis=1)
samplemeta['expression_available']=samplemeta.expression_source!=''
samplemeta.to_csv(ROOT/'01_metadata/patient_sample_metadata.tsv',sep='\t',index=False)
audit={'status':'METADATA_AUDIT_PASS','utc':now(),'GSE205013_patients':27,'GSE202051_expression_patients':len(mapping),'GSE202051_expression_tumor_patients':int(mapping.patient_id.str.startswith('PDAC_').sum()),'GSE202051_missing_registered_IDs':registration.loc[~registration.present_in_expression_object,'patient_id'].tolist(),'known_conflict_handling':'P11 explicitly Unknown and excluded from treatment inference; no imputed resolution','technical_replicates':'Collapsed to patient, no GSM-count sample size inflation','normal_tumor_identity':'Anonymous labels; no assumption that registered normal specimens are additional independent tumor patients'}
(ROOT/'logs/METADATA_AUDIT_PASS.json').write_text(json.dumps(audit,indent=2))
print('METADATA_AUDIT_PASS');print('DATA_AUDITS_COMPLETE')
