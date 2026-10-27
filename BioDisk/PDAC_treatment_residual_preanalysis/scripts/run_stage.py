from pathlib import Path
import sys,subprocess,os,json,datetime,hashlib,py_compile
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
stage,script,marker,*args=sys.argv[1:]
if '--background' in args:
    args.remove('--background')
    with (ROOT/'logs'/f'{stage}.launcher.log').open('a') as out:
        p=subprocess.Popen([PY,__file__,stage,script,marker]+args,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=out,stderr=out,start_new_session=True)
    (ROOT/'logs'/f'{stage}.pid').write_text(str(p.pid));print('STARTED',stage,p.pid);sys.exit(0)
path=ROOT/'scripts'/script
env=os.environ.copy()
env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4',http_proxy='http://127.0.0.1:7897',https_proxy='http://127.0.0.1:7897')
record={'stage':stage,'script':str(path),'script_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'EXPECTED_SUCCESS_MARKER':marker}
if path.suffix=='.py':
    py_compile.compile(str(path),doraise=True);cmd=[PY,'-u',str(path)]+args
else:
    p=subprocess.run([R,'--vanilla',str(ROOT/'scripts/check_R_syntax.R'),str(path)],capture_output=True,text=True);p.check_returncode();cmd=[R,'--vanilla',str(path)]+args
record['SYNTAX_CHECK']='PASS';record['RUN_COMMAND']=cmd
with (ROOT/'logs'/f'{stage}.log').open('w') as out:
    out.write(json.dumps(record)+'\n');out.flush()
    p=subprocess.run(cmd,cwd=ROOT,env=env,stdout=out,stderr=out)
record['exit_code']=p.returncode
record['marker_seen']=marker in (ROOT/'logs'/f'{stage}.log').read_text().split('\n',1)[-1]
record['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(ROOT/'logs'/f'{stage}.status.json').write_text(json.dumps(record,indent=2))
with (ROOT/'logs/execution_ledger.jsonl').open('a') as f: f.write(json.dumps(record)+'\n')
print(json.dumps(record));sys.exit(0 if p.returncode==0 and record['marker_seen'] else 1)
