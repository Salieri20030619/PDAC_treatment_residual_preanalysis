from pathlib import Path
import subprocess,re,json,hashlib,datetime,concurrent.futures
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
PROXY='http://127.0.0.1:7897'
TASKS=[
 ('Werba2023.html','https://www.nature.com/articles/s41467-023-36296-4'),
 ('Werba2023_PMC.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC9925748/'),
 ('Hwang2022.html','https://www.nature.com/articles/s41588-022-01134-8'),
 ('hallmark_2025.1.gmt','https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt'),
 ('reactome_2025.1.gmt','https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c2.cp.reactome.v2025.1.Hs.symbols.gmt'),
 ('hg38_refGene.txt.gz','https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz')]
def get(task):
    name,url=task; dest=ROOT/'00_docs'/name
    cmd=['curl','--proxy',PROXY,'-fL','--retry','2','--connect-timeout','20','--max-time','240','-o',str(dest),url]
    if dest.exists() and dest.stat().st_size>1000: p=None
    else: p=subprocess.run(cmd,capture_output=True,text=True)
    r=dict(file=str(dest),url=url,exit_code=p.returncode if p else 0,utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    if r['exit_code']==0:
        r.update(bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
        if name.endswith('.html'):
            text=dest.read_text()
            links=sorted(set(re.findall(r'https://[^"<>\s]+(?:\.xlsx|\.pdf)',text)))
            (ROOT/'00_docs'/f'{name}.supplement_links.json').write_text(json.dumps(links,indent=2))
            r['supplement_links']=links
            plain=re.sub('<[^>]+>',' ',text)
            (ROOT/'00_docs'/name.replace('.html','.txt')).write_text(plain)
    else: r['error']=p.stderr[-1500:]
    print(json.dumps(r),flush=True); return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as p: out=list(p.map(get,TASKS))
(ROOT/'01_metadata/reference_download_manifest.json').write_text(json.dumps(out,indent=2))
print('REFERENCE_FETCH_COMPLETE')
