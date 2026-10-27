from pathlib import Path
import re,json,collections
import pandas as pd
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
def parse(path):
    samples=[]; current=None; series=collections.defaultdict(list)
    for n,line in enumerate(path.read_text().splitlines(),1):
        if line.startswith('^SAMPLE = '):
            current=collections.defaultdict(list); current['GSM']=[line.split(' = ')[1]]; samples.append(current)
        elif line.startswith('^'): current=None
        elif ' = ' in line:
            key,value=line.split(' = ',1)
            if current is not None: current[key].append(value); current[key+'_line'].append(n)
            elif key.startswith('!Series'): series[key].append(value)
    return samples,series
rows=[]; evidence=[]
for ds in ['GSE205013','GSE202051']:
    path=ROOT/f'01_metadata/{ds}_family.soft'; samples,series=parse(path)
    (ROOT/f'01_metadata/{ds}_series.json').write_text(json.dumps(series,indent=2))
    for s in samples:
        gsm=s['GSM'][0]; title=s['!Sample_title'][0]
        ch={x.split(':',1)[0].lower():x.split(':',1)[1].strip() for x in s['!Sample_characteristics_ch1'] if ':' in x}
        tissue=ch.get('tissue',''); treat=ch.get('treatment','')
        patient=title if ds=='GSE205013' else re.split(r',?\s*replicate',title)[0].strip()
        organoid=ch.get('cell line')=='Patient-derived organoid'
        row=dict(dataset=ds,patient_id=patient if not organoid else '',sample_id=patient,GSM=gsm,tissue_site=tissue,treatment_status=('Untreated' if treat=='Untreated' else 'Treated' if treat else ''),regimen=treat if ds=='GSE202051' and treat!='Untreated' else '',sampling_time='',response='',procedure=ch.get('procedure',''),stage=ch.get('stage',''),preservation='',assay='scRNAseq' if ds=='GSE205013' else 'snRNAseq',source_file=path.name,source_url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}',title=title,regimen_raw=treat,specimen_type='Organoid' if organoid else ch.get('cell type','Tumor'),replicate=re.search(r'replicate (\d+)',title).group(1) if 'replicate ' in title else '',patient_id_basis='GEO sample title; clinical-table verification pending',technical_replicate_basis='GEO title; not independent patients',comparison_scope='Excluded_organoid' if organoid else 'Normal_reference' if ch.get('cell type')=='Normal' else 'Liver_descriptive' if 'liver' in tissue.lower() else 'Primary')
        if ds=='GSE202051' and not organoid:
            row['procedure']='Resection'; row['preservation']='Snap frozen'
        row['unknown_fields']=';'.join(k for k in ['patient_id','regimen','sampling_time','response','procedure','stage','preservation'] if not row[k])
        rows.append(row)
        for key in ['!Sample_title','!Sample_characteristics_ch1','!Sample_extract_protocol_ch1','!Sample_description','!Sample_data_processing']:
            for val,line in zip(s[key],s[key+'_line']): evidence.append(dict(dataset=ds,GSM=gsm,source_file=path.name,line=line,field=key,value=val))
df=pd.DataFrame(rows)
df.to_csv(ROOT/'01_metadata/sample_metadata_GEO.tsv',sep='\t',index=False)
pd.DataFrame(evidence).to_csv(ROOT/'01_metadata/metadata_evidence_GEO.tsv',sep='\t',index=False)
patients=df[df.patient_id!=''].drop_duplicates(['dataset','patient_id'])
for cols,name in [(['dataset','comparison_scope','procedure','treatment_status'],'group'),(['dataset','comparison_scope','stage','treatment_status'],'stage'),(['dataset','comparison_scope','regimen_raw'],'regimen')]:
    cross=patients.groupby(cols,dropna=False).size().reset_index(name='n_candidate_patients')
    cross.to_csv(ROOT/f'06_tables/GEO_patient_{name}_crosstab.tsv',sep='\t',index=False)
    print(cross.to_string(index=False))
counts=df.groupby(['dataset','sample_id','specimen_type']).agg(n_GSM=('GSM','size'),GSM=('GSM',';'.join),treatments=('regimen_raw',lambda x:';'.join(sorted(set(x))))).reset_index()
counts.to_csv(ROOT/'06_tables/GEO_replicate_audit.tsv',sep='\t',index=False)
print('GEO_METADATA_PARSE_COMPLETE; clinical and expression mapping still pending')
