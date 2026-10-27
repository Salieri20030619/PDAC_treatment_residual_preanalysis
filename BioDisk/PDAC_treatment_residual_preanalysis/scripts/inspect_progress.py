from pathlib import Path
import json,tarfile,gzip,shutil
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
for ds in ['GSE205013','GSE202051']:
    folder=ROOT/'02_raw'/ds/'segments'
    files=sorted(folder.glob('*.bin'))
    print(ds,'downloaded_GiB',round(sum(p.stat().st_size for p in files)/1024**3,3))
    if ds!='GSE205013': continue
    plan=json.loads((folder/'plan.json').read_text()); prefix=folder/'prefix.bin'
    paths=[prefix]; end=prefix.stat().st_size
    for p in sorted(folder.glob('*_*.bin')):
        start,stop=map(int,p.stem.split('_'))
        if start!=end: break
        paths.append(p);end+=p.stat().st_size
        if p.stat().st_size!=stop-start+1: break
    print('CONTIGUOUS_TAR_BYTES',end)
    target=ROOT/'02_raw'/ds/'early_prefix.tar'
    with target.open('wb') as dst:
        for p in paths:
            with p.open('rb') as src: shutil.copyfileobj(src,dst,8*1024*1024)
    out=ROOT/'02_raw'/ds/'early_matrices';out.mkdir(exist_ok=True)
    try:
        with tarfile.open(target) as tar:
            for m in tar:
                if not m.isfile(): continue
                if m.offset_data+m.size>end: break
                assert Path(m.name).name==m.name
                dest=out/m.name
                if not dest.exists():
                    with tar.extractfile(m) as src,dest.open('wb') as dst: shutil.copyfileobj(src,dst)
                    with gzip.open(dest) as src:
                        while src.read(8*1024*1024):pass
                print('READY_MEMBER',m.name,m.size)
    except tarfile.ReadError as e: print('PARTIAL_ARCHIVE_EXPECTED',str(e))
print('INITIAL_PROGRESS_INSPECTION_COMPLETE')
