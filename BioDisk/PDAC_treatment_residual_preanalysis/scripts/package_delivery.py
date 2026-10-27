from pathlib import Path
import zipfile,json,hashlib,datetime
R=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
assert json.loads((R/'logs/PREANALYSIS_COMPLETE.json').read_text())['status']=='PREANALYSIS_COMPLETE'
out=R/'08_reports/PDAC_client_delivery_20260922.zip'
paths=[]
for folder in ['01_metadata','06_tables','07_figures','08_reports']:
    for p in sorted((R/folder).iterdir()):
        if p.is_file() and (p.suffix in ['.tsv','.md','.html','.json','.png'] or p.name.endswith('.tsv.gz')) and p.name!='delivery_archive.json':paths.append(p)
for p in (R/'scripts').glob('*'):
    if p.is_file() and p.suffix in ['.py','.R','.sh']:paths.append(p)
for p in (R/'00_docs').iterdir():
    if p.is_file() and p.suffix in ['.md','.json']:paths.append(p)
for name in ['DOWNLOAD_AUDIT_PASS.json','METADATA_AUDIT_PASS.json','PREANALYSIS_COMPLETE.json','execution_ledger.jsonl','final_validation_checks.json']:
    paths.append(R/'logs'/name)
paths.append(R/'05_analysis/scientific_interpretation.json')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(set(paths)):z.write(p,str(p.relative_to(R)))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
h=hashlib.sha256(out.read_bytes()).hexdigest()
record={'file':out.name,'bytes':out.stat().st_size,'sha256':h,'n_files':len(set(paths)),'status':'DELIVERY_ARCHIVE_PASS','scope':'Client proposal, evidence reports, tables, figures, scripts and selected audit logs; raw expression and cell objects remain on sgs10'}
(R/'08_reports/delivery_archive.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record));print('DELIVERY_ARCHIVE_PASS')
