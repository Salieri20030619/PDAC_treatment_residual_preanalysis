from pathlib import Path
import subprocess, os, sys, json, datetime, hashlib, gzip, tarfile, concurrent.futures, traceback
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
PROXY='http://127.0.0.1:7897'
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def download(ds, name, perform_transfer=True):
    url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{ds[:-3]}nnn/{ds}/suppl/{name}'
    dest=ROOT/'02_raw'/ds/name
    record=ROOT/'logs'/f'{name}.download.json'
    rec={'url':url,'file':str(dest),'proxy':PROXY,'started_utc':now(),'official_checksum':None,'official_checksum_note':'No checksum listed in saved GEO supplementary directory index','status':'STARTED'}
    if record.exists() and dest.exists():
        old=json.loads(record.read_text())
        if old.get('status')=='INTEGRITY_PASS' and old.get('bytes')==dest.stat().st_size:
            print('SKIP_COMPLETE',name,flush=True); return old
    record.write_text(json.dumps(rec,indent=2))
    try:
        with open(ROOT/'logs'/f'{name}.curl.log','a') as log:
            head=subprocess.run(['curl','--proxy',PROXY,'-fsSI','--connect-timeout','20','--max-time','90',url],capture_output=True,text=True)
            rec['head_response']=head.stdout
            (ROOT/'logs'/f'{name}.headers.txt').write_text(head.stdout+head.stderr)
            for attempt in range(1,4) if perform_transfer else []:
                rec['attempt']=attempt
                command=['curl','--proxy',PROXY,'-fL','-C','-','--connect-timeout','30','--max-time','14400','--speed-limit','1024','--speed-time','180','--retry','2','--retry-delay','15','--retry-max-time','1800','-o',str(dest)+'.part',url]
                log.write('\nRUN_COMMAND '+str(command)+'\n'); log.flush()
                proc=subprocess.run(command,stdout=log,stderr=log)
                rec['curl_exit']=proc.returncode
                if proc.returncode==0: break
            if perform_transfer and proc.returncode: raise RuntimeError(f'curl failed after bounded attempts: {proc.returncode}')
        part=Path(str(dest)+'.part')
        rec['bytes']=part.stat().st_size
        sizes=[int(l.split(':',1)[1].strip()) for l in rec['head_response'].splitlines() if l.lower().startswith('content-length:')]
        if sizes and sizes[-1]!=rec['bytes']: raise ValueError('Content-Length mismatch')
        h=hashlib.sha256()
        with part.open('rb') as f:
            for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
        rec['sha256']=h.hexdigest()
        if name.endswith('.tar'):
            out=dest.parent/'matrices'; out.mkdir(exist_ok=True)
            members=[]
            with tarfile.open(part) as tar:
                for member in tar:
                    if not member.isfile(): continue
                    if Path(member.name).name!=member.name: raise ValueError('Unexpected archive path')
                    with tar.extractfile(member) as src, (out/member.name).open('wb') as dst:
                        while block:=src.read(8*1024*1024): dst.write(block)
                    if member.name.endswith('.gz'):
                        with gzip.open(out/member.name,'rb') as g:
                            while g.read(8*1024*1024): pass
                    members.append({'name':member.name,'bytes':member.size,'gzip':'PASS'})
            (ROOT/'01_metadata'/f'{ds}_archive_members.json').write_text(json.dumps(members,indent=2))
            rec['archive_members']=len(members); rec['tar_integrity']='PASS'; rec['member_gzip_integrity']='PASS'
        else:
            target=dest.with_suffix('')
            with gzip.open(part,'rb') as src, Path(str(target)+'.part').open('wb') as dst:
                while block:=src.read(8*1024*1024): dst.write(block)
            Path(str(target)+'.part').rename(target)
            rec['gzip_integrity']='PASS'
            import h5py
            with h5py.File(target,'r') as h:
                rec['h5ad_root_keys']=list(h.keys())
                assert 'obs' in h and 'var' in h and 'X' in h
            rec['hdf5_open']='PASS'
        part.rename(dest)
        rec['status']='INTEGRITY_PASS'
    except Exception as e:
        rec['status']='FAILED'; rec['error']=str(e); rec['traceback']=traceback.format_exc()
    rec['finished_utc']=now(); record.write_text(json.dumps(rec,indent=2)); print(json.dumps(rec),flush=True)
    return rec
if __name__=='__main__':
    tasks=[('GSE205013','GSE205013_RAW.tar'),('GSE202051','GSE202051_totaldata-final-toshare.h5ad.gz')]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool: results=list(pool.map(lambda x:download(*x),tasks))
    (ROOT/'01_metadata/expression_download_manifest.json').write_text(json.dumps(results,indent=2))
    ok=all(r['status']=='INTEGRITY_PASS' for r in results)
    (ROOT/'logs/download_data.exit_status').write_text('0' if ok else '1')
    print('EXPRESSION_DOWNLOAD_INTEGRITY_COMPLETE' if ok else 'EXPRESSION_DOWNLOAD_PARTIAL',flush=True)
    sys.exit(0 if ok else 1)
