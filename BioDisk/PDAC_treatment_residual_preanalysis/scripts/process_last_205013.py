from pathlib import Path
import sys,subprocess
from stream_samples import extract_ready,ROOT
ready=extract_ready()
assert 'P27' in ready,'Final sample is not yet complete; no partial matrix will be used'
subprocess.run([sys.executable,str(ROOT/'scripts/run_stage.py'),'sample_P27_early','process_205013.py','SAMPLE_BATCH_COMPLETE','P27'],cwd=ROOT,check=True)
print('FINAL_RAW_SAMPLE_PROCESSED')
