from pathlib import Path
import os, sys, subprocess, json, datetime, importlib.metadata, gzip, hashlib
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
DIRS='00_docs 01_metadata 02_raw/GSE205013 02_raw/GSE202051 03_processed 04_qc 05_analysis 06_tables 07_figures 08_reports scripts logs'.split()
for d in DIRS: (ROOT/d).mkdir(parents=True,exist_ok=True)
os.environ.update(http_proxy='http://127.0.0.1:7897',https_proxy='http://127.0.0.1:7897',HTTP_PROXY='http://127.0.0.1:7897',HTTPS_PROXY='http://127.0.0.1:7897')
def run(args):
    p=subprocess.run(args,capture_output=True,text=True)
    return {'command':args,'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
env={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.executable,'proxy':'http://127.0.0.1:7897','packages':{}}
for name in ['scanpy','anndata','numpy','pandas','scipy','h5py','matplotlib','seaborn','statsmodels','scrublet','leidenalg','igraph','infercnvpy','gseapy','requests']:
    try: env['packages'][name]=importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError: env['packages'][name]='MISSING'
env['system']=[run(c) for c in [['df','-h',str(ROOT)],['free','-h'],['nproc'],['uptime']]]
(ROOT/'logs/environment.json').write_text(json.dumps(env,indent=2))
print(json.dumps(env,indent=2),flush=True)
manifest=[]
for ds in ['GSE205013','GSE202051']:
    base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{ds[:-3]}nnn/{ds}'
    for suffix,target in [(f'/soft/{ds}_family.soft.gz',ROOT/f'01_metadata/{ds}_family.soft.gz'),('/suppl/',ROOT/f'01_metadata/{ds}_suppl_index.html')]:
        url=base+suffix
        p=run(['curl','--proxy',os.environ['https_proxy'],'-fL','--retry','3','--retry-delay','5','--connect-timeout','20','--max-time','180','-o',str(target),url])
        rec={'url':url,'file':str(target),'exit_code':p['exit_code'],'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stderr':p['stderr']}
        if p['exit_code']==0:
            rec['bytes']=target.stat().st_size
            rec['sha256']=hashlib.sha256(target.read_bytes()).hexdigest()
            if target.suffix=='.gz':
                contents=gzip.decompress(target.read_bytes())
                target.with_suffix('').write_bytes(contents)
                rec['gzip_integrity']='PASS'
        manifest.append(rec)
        print(json.dumps(rec),flush=True)
(ROOT/'01_metadata/initial_download_manifest.json').write_text(json.dumps(manifest,indent=2))
print('BOOTSTRAP_COMPLETE',flush=True)
