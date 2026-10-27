from pathlib import Path
import requests,json,re,hashlib,datetime
from bs4 import BeautifulSoup
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
tasks=[('PXD017051.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD017051'),('PXD020437.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD020437'),('metabolomics_2025.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC11793008/')]
manifest=[]
for name,url in tasks:
    dest=ROOT/'00_docs'/name
    rec={'url':url,'file':str(dest),'scope':'Next-stage candidate metadata only; no omics raw-data download','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        r=requests.get(url,proxies=proxies,timeout=(20,60));r.raise_for_status();dest.write_bytes(r.content)
        rec.update(status='SOURCE_SAVED',bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
        if name.endswith('.html'):
            text=BeautifulSoup(r.text,'html.parser').get_text(' ',strip=True);dest.with_suffix('.txt').write_text(text)
            for pattern in ['data availability','data access','metabolights','metabolomics workbench','MTBLS']:
                for m in list(re.finditer(pattern,text,re.I))[:3]:print(name,pattern,text[max(0,m.start()-100):m.end()+1200],flush=True)
        else:
            data=r.json();print(name,json.dumps({k:v for k,v in data.items() if any(s in k.lower() for s in ['title','projectdescription','submission','publication','accession'])},ensure_ascii=False)[:3500],flush=True)
    except Exception as e:rec.update(status='CANDIDATE_SOURCE_FETCH_FAILED',error=repr(e))
    manifest.append(rec)
(ROOT/'01_metadata/next_omics_source_manifest.json').write_text(json.dumps(manifest,indent=2));print('NEXT_OMICS_SOURCE_REVIEW_COMPLETE')
