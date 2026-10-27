from pathlib import Path
import os,select,json,subprocess,sys,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
stage=sys.argv[1] if len(sys.argv)>1 else 'segmented_resume'
pid=int((ROOT/'logs'/f'{stage}.pid').read_text())
if hasattr(os,'pidfd_open'):
    try:
        fd=os.pidfd_open(pid)
        try:select.select([fd],[],[])
        finally:os.close(fd)
    except ProcessLookupError:pass
else:
    subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
status=json.loads((ROOT/'logs'/f'{stage}.status.json').read_text())
assert status['exit_code']==0 and status['marker_seen'],status
subprocess.run([sys.executable,str(ROOT/'scripts/run_stage.py'),'h5ad_audit','audit_h5ad.py','H5AD_STRUCTURE_AUDIT_COMPLETE'],cwd=ROOT,check=True)
print('DOWNLOAD_WAIT_AND_H5AD_AUDIT_COMPLETE',datetime.datetime.now().isoformat())
