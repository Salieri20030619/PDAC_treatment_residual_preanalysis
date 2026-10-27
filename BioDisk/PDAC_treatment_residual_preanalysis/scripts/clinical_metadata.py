from pathlib import Path
import pandas as pd,json,re,hashlib,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
df=pd.read_csv(ROOT/'01_metadata/sample_metadata_GEO.tsv',sep='\t',keep_default_na=False)
df['treatment_status_GEO']=df.treatment_status;df['metadata_conflict']='';df['radiotherapy']='';df['clinical_source']='';df['treatment_duration_reported']='';df['regimen_group']='';df['chemistry']=''
df['procedure']=df.procedure.replace({'Liver-BIopsy':'Liver-Biopsy'})
text=(ROOT/'00_docs/41467_2023_36296_MOESM1_ESM.txt').read_text()
clinical=[]
for line in text.splitlines():
    match=re.match(r'^\s*(P\d{2})\s{2,}(.*?)\s{2,}(.*?)\s{2,}(PR|PD|SD|Non-evaluable)\s*$',line)
    if match:clinical.append(match.groups())
assert len(clinical)==8 and len(set(x[0] for x in clinical))==8
for patient,reg,duration,response in clinical:
    sel=(df.dataset=='GSE205013')&(df.patient_id==patient)
    df.loc[sel,'regimen']=reg.replace('FFX','FOLFIRINOX').replace('G/A','gemcitabine/nab-paclitaxel')
    df.loc[sel,'regimen_group']='FFX+SBRT' if 'SBRT' in reg else 'FFX+GA' if ',' in reg else 'FFX' if reg=='FFX' else 'GA'
    df.loc[sel,'treatment_duration_reported']=duration;df.loc[sel,'response']=response
    df.loc[sel,'radiotherapy']='Yes' if 'SBRT' in reg else 'Not_listed'
    df.loc[sel,'clinical_source']='Werba2023 Supplementary Table 2; PDF page 24; visually verified'
    if (df.loc[sel,'treatment_status_GEO']=='Untreated').any():
        df.loc[sel,'metadata_conflict']='GEO untreated vs clinical table treated; timing at sequenced sampling unresolved'
        df.loc[sel,'treatment_status']='Unknown'
df.loc[(df.dataset=='GSE205013')&(df.treatment_status=='Untreated'),'regimen_group']='Untreated'
df.loc[df.dataset=='GSE205013','preservation']='Fresh'
df.loc[df.dataset=='GSE205013','patient_id_basis']='GEO Pxx matched published patient mutation table and study n=27 unique patients'
sheets=json.loads((ROOT/'01_metadata/Hwang2022_supplement_sheets.json').read_text())
table=sheets['Supplementary Table 1'];clinical202=[]
expansion={'CRT':'FOLFIRINOX + radiotherapy + concurrent capecitabine or 5-FU','CRTL':'FOLFIRINOX + losartan + radiotherapy + concurrent capecitabine or 5-FU','CRTN':'FOLFIRINOX + SBRT + nivolumab','CRTLN':'FOLFIRINOX + SBRT + losartan + nivolumab'}
for n,r in enumerate(table,1):
    patient=r[0].strip()
    if not patient.startswith('PDAC_'):continue
    rec=dict(zip(table[1],r));rec['source_excel_row']=n;clinical202.append(rec)
    if '\u2020' in patient:continue
    sel=(df.dataset=='GSE202051')&(df.patient_id==patient)
    if not sel.any():continue
    reg=r[6].replace('*','')
    df.loc[sel,'stage']=r[3];df.loc[sel,'chemistry']=r[11]
    df.loc[sel,'regimen_group']='Untreated' if reg=='None' else reg
    if reg in expansion:df.loc[sel,'regimen']=expansion[reg]
    if reg!='None':df.loc[sel,'radiotherapy']='Yes'
    df.loc[sel,'clinical_source']=f'Hwang2022 Supplementary Table 1; row {n}'
    df.loc[sel,'patient_id_basis']='GEO specimen ID matched unique clinical cohort ID; replicate GSM collapsed'
df['unknown_fields']=df.apply(lambda r:';'.join(k for k in ['patient_id','regimen','sampling_time','response','procedure','stage','preservation'] if not r[k]),axis=1)
df.to_csv(ROOT/'01_metadata/patient_sample_metadata.tsv',sep='\t',index=False)
pd.DataFrame(clinical,columns=['patient_id','reported_regimen','duration','response']).to_csv(ROOT/'01_metadata/Werba_clinical_table2_extracted.tsv',sep='\t',index=False)
pd.DataFrame(clinical202).to_csv(ROOT/'01_metadata/Hwang_clinical_table1_extracted.tsv',sep='\t',index=False)
df[df.metadata_conflict!=''].to_csv(ROOT/'01_metadata/metadata_conflicts.tsv',sep='\t',index=False)
patients=df[df.patient_id!=''].drop_duplicates(['dataset','patient_id'])
patients.to_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',index=False)
for key in ['procedure','stage','regimen_group','radiotherapy','chemistry']:
    tab=patients.groupby(['dataset','comparison_scope',key,'treatment_status'],dropna=False).size().reset_index(name='n_patients_registered')
    tab.to_csv(ROOT/f'06_tables/patient_{key}_crosstab.tsv',sep='\t',index=False)
    print(key,tab.to_string(index=False),flush=True)
note={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'GSE205013_registered_patients':int((patients.dataset=='GSE205013').sum()),'GSE202051_registered_tumor_patients':int(((patients.dataset=='GSE202051')&(patients.specimen_type=='Tumor')).sum()),'expression_mapping':'PENDING','P11':'Excluded from treatment inference; clinical conflict retained','GSE202051_Normal':'6 normal IDs registered; source says aggregate 45 specimens, actual object inclusion pending','cohort_overlap':'Distinct primary publications and recruiting institutions (NYU vs MGH); no overlapping patient identifiers or reported cohort reuse found; anonymous IDs prevent genetic identity verification','EGA':'Excluded entirely','audit_status':'Clinical reconciliation complete; expression audit pending'}
(ROOT/'01_metadata/clinical_reconciliation.json').write_text(json.dumps(note,indent=2))
print('CLINICAL_RECONCILIATION_COMPLETE')
