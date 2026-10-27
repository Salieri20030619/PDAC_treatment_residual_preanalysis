from pathlib import Path
import io,tarfile,gzip,json,os,subprocess,time,concurrent.futures,sys,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
class SegmentedReader(io.RawIOBase):
    def __init__(self,paths):
        self.parts=[];self.pos=0;self.length=0
        for p in paths:
            size=p.stat().st_size;self.parts.append((self.length,self.length+size,p));self.length+=size
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.pos
    def seek(self,offset,whence=0):
        self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.length+offset
        return self.pos
    def read(self,n=-1):
        if n<0:n=self.length-self.pos
        n=min(n,self.length-self.pos)
        out=[]
        for start,end,p in self.parts:
            if n<=0:break
            if start<=self.pos<end:
                take=min(n,end-self.pos)
                with p.open('rb') as f:f.seek(self.pos-start);b=f.read(take)
                out.append(b);self.pos+=len(b);n-=len(b)
        return b''.join(out)
def extract_ready():
    folder=ROOT/'02_raw/GSE205013/segments';prefix=folder/'prefix.bin';paths=[prefix];end=prefix.stat().st_size
    for p in sorted(folder.glob('*_*.bin')):
        start,stop=map(int,p.stem.split('_'))
        if start!=end:break
        size=p.stat().st_size;paths.append(p);end+=size
        if size!=stop-start+1:break
    out=ROOT/'02_raw/GSE205013/early_matrices';out.mkdir(exist_ok=True)
    reader=SegmentedReader(paths)
    try:
        with tarfile.open(fileobj=reader,mode='r:') as tar:
            for m in tar:
                if not m.isfile():continue
                if m.offset_data+m.size>end:break
                assert Path(m.name).name==m.name
                dest=out/m.name
                if dest.exists():continue
                part=Path(str(dest)+'.part')
                with tar.extractfile(m) as src,part.open('wb') as dst:
                    while b:=src.read(8*1024*1024):dst.write(b)
                with gzip.open(part) as f:
                    while f.read(8*1024*1024):pass
                part.rename(dest)
    except tarfile.ReadError:pass
    return sorted({p.name.split('_')[1] for p in out.glob('*matrix.mtx.gz')})
def run_sample(sample):
    p=subprocess.run([PY,str(ROOT/'scripts/run_stage.py'),'sample_'+sample,'process_205013.py','SAMPLE_BATCH_COMPLETE',sample],cwd=ROOT,capture_output=True,text=True)
    print(p.stdout,p.stderr,flush=True);return sample,p.returncode
def main():
    own={'P01','P02'} # already launched, don't launch duplicates
    attempted=set();jobs={};last_scan=0
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        while True:
            completed={p.name.split('_')[1].split('.')[0] for p in (ROOT/'03_processed').glob('GSE205013_P*.complete.json')}
            if len(completed)==27:break
            if time.time()-last_scan>=900:
                ready=extract_ready();last_scan=time.time()
                print('CHECKPOINT',datetime.datetime.now().isoformat(),'ready',len(ready),'complete',len(completed),flush=True)
                for sample in ready:
                    if sample not in completed|attempted|own:
                        jobs[sample]=pool.submit(run_sample,sample);attempted.add(sample)
                (ROOT/'logs/stream_samples_progress.json').write_text(json.dumps({'ready':ready,'complete':sorted(completed),'submitted':sorted(attempted),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
            failed=[k for k,v in jobs.items() if v.done() and v.result()[1]!=0]
            if failed:raise RuntimeError('Sample failures requiring review: '+','.join(failed))
            if '--once' in sys.argv:break
            time.sleep(15)
    print('STREAM_SAMPLE_BATCH_COMPLETE',flush=True)
if __name__=='__main__':main()
