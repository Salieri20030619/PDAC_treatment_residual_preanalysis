from pathlib import Path
import subprocess,sys,json
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
stage='identity_batch6';pid=(ROOT/'logs'/f'{stage}.pid').read_text().strip()
subprocess.run(['tail','--pid='+pid,'-f','/dev/null'],check=True)
status=json.loads((ROOT/'logs'/f'{stage}.status.json').read_text());assert status['exit_code']==0 and status['marker_seen']
assert len(list((ROOT/'03_processed').glob('GSE205013_P*_identity_final.json')))==27
for stage,script,marker in [('discovery_count_integrity','audit_pseudobulk_integrity.py','COUNT_INTEGRITY_COMPLETE'),('discovery_QC_final','review_qc.py','QC_SUMMARY_COMPLETE'),('discovery_CNV_figure','figure_cnv_QA.py','CNV_QA_FIGURE_COMPLETE')]:
    subprocess.run([sys.executable,str(ROOT/'scripts/run_stage.py'),stage,script,marker],cwd=ROOT,check=True)
print('DISCOVERY_AGGREGATION_COMPLETE')
