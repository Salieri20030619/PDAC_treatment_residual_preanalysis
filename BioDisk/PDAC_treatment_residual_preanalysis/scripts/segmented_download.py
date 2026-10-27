from pathlib import Path
import os,signal,time,json,concurrent.futures,requests,datetime,sys,hashlib
from download_data import ROOT, PROXY, download
TASKS=[('GSE205013','GSE205013_RAW.tar'),('GSE202051','GSE202051_totaldata-final-toshare.h5ad.gz')]
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def stop_initial():
    own=os.getpid(); stopped=[]
    for d in Path('/proc').iterdir():
        if not d.name.isdigit() or int(d.name)==own: continue
        try:
            cmd=(d/'cmdline').read_bytes().split(b'\0')
            args=[a.decode(errors='replace') for a in cmd]
            is_worker='scripts/download_data.py' in args and any('python' in a for a in args[:1])
            is_curl=args and Path(args[0]).name=='curl' and any(str(ROOT/'02_raw') in a for a in args)
            if is_worker or is_curl: os.kill(int(d.name),signal.SIGTERM); stopped.append({'pid':int(d.name),'command':args})
        except (PermissionError,FileNotFoundError,ProcessLookupError): pass
    (ROOT/'logs/download_transfer_handoff.json').write_text(json.dumps({'utc':now(),'reason':'Resume existing bytes with HTTP Range parallel segments via the same proxy','stopped':stopped},indent=2))
    time.sleep(2)
def get_range(job):
    url,start,end,dest=job
    size=end-start+1
    for attempt in range(1,5):
        offset=dest.stat().st_size if dest.exists() else 0
        if offset==size: return str(dest)
        if offset>size: raise ValueError('Oversize segment')
        try:
            with requests.get(url,headers={'Range':f'bytes={start+offset}-{end}'},proxies={'http':PROXY,'https':PROXY},stream=True,timeout=(30,120)) as r:
                r.raise_for_status()
                if r.status_code!=206: raise ValueError(f'Range request returned {r.status_code}')
                expected=f'bytes {start+offset}-{end}/'
                if not r.headers.get('Content-Range','').startswith(expected): raise ValueError('Content-Range mismatch')
                with dest.open('ab') as f:
                    for b in r.iter_content(1024*1024): f.write(b)
            if dest.stat().st_size!=size: raise ValueError('Incomplete segment')
            print('SEGMENT_COMPLETE',dest.name,size,now(),flush=True); return str(dest)
        except Exception as e:
            print('SEGMENT_RETRY',dest.name,attempt,repr(e),flush=True)
            if attempt==4: raise
            time.sleep(20*attempt)
def main():
    stop_initial(); infos=[]; jobs=[]
    for ds,name in TASKS:
        dest=ROOT/'02_raw'/ds/name
        if dest.exists():
            rec=json.loads((ROOT/'logs'/f'{name}.download.json').read_text())
            if rec.get('status')=='INTEGRITY_PASS': continue
        folder=dest.parent/'segments'; folder.mkdir(exist_ok=True)
        url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{ds[:-3]}nnn/{ds}/suppl/{name}'
        if (folder/'plan.json').exists():
            oldplan=json.loads((folder/'plan.json').read_text());assert oldplan['url']==url;size=int(oldplan['total_bytes'])
        else:
            r=requests.head(url,proxies={'http':PROXY,'https':PROXY},timeout=(20,60)); r.raise_for_status(); size=int(r.headers['Content-Length'])
        prefix=folder/'prefix.bin'; part=Path(str(dest)+'.part')
        if part.exists() and not prefix.exists(): part.rename(prefix)
        prefix_size=prefix.stat().st_size if prefix.exists() else 0
        pieces=[]
        for start in range(prefix_size,size,16*1024*1024):
            end=min(size-1,start+16*1024*1024-1); piece=folder/f'{start:012d}_{end:012d}.bin'
            pieces.append(piece); jobs.append((url,start,end,piece))
        infos.append((ds,name,size,prefix,pieces))
        (folder/'plan.json').write_text(json.dumps({'url':url,'total_bytes':size,'prefix_bytes':prefix_size,'segment_count':len(pieces),'proxy':PROXY,'utc':now()},indent=2))
    # Interleave datasets to avoid starving validation cohort.
    jobs.sort(key=lambda x:x[1])
    failed=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        pending={pool.submit(get_range,j):j for j in jobs}
        for fut in concurrent.futures.as_completed(pending):
            try:fut.result()
            except Exception as e:
                failed.append(pending[fut]);print('SEGMENT_DEFERRED',pending[fut][3].name,repr(e),flush=True)
    # One bounded late retry pass permits recovery from an intermittent proxy outage.
    if failed:
        (ROOT/'logs/segments_deferred.json').write_text(json.dumps([str(j[3]) for j in failed],indent=2))
        print('FINAL_BOUNDED_RETRY_PASS',len(failed),flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            for _ in pool.map(get_range,failed):pass
    results=[]
    for ds,name,size,prefix,pieces in infos:
        part=ROOT/'02_raw'/ds/(name+'.part')
        with part.open('wb') as out:
            for piece in ([prefix] if prefix.exists() else [])+pieces:
                with piece.open('rb') as f:
                    while b:=f.read(8*1024*1024): out.write(b)
        assert part.stat().st_size==size
        results.append(download(ds,name,False))
    (ROOT/'01_metadata/expression_download_manifest.json').write_text(json.dumps(results,indent=2))
    assert all(r['status']=='INTEGRITY_PASS' for r in results)
    print('EXPRESSION_DOWNLOAD_INTEGRITY_COMPLETE',flush=True)
if __name__=='__main__':
    code=1
    try: main(); code=0
    finally: (ROOT/'logs/segmented_download.exit_status').write_text(str(code))
