from pathlib import Path
import sys,subprocess,concurrent.futures,json
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
def stage(name,script,marker):
    command=[sys.executable,str(ROOT/'scripts/run_stage.py'),name,script,marker]
    p=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
    print(p.stdout,p.stderr,flush=True);p.check_returncode()
selection=json.loads((ROOT/'logs/select_validation_layers.status.json').read_text())
assert selection['exit_code']==0 and selection['marker_seen']
stage('validation_processing','process_202051.py','VALIDATION_PROCESSING_COMPLETE')
stage('final_data_audits','finalize_data_audits.py','DATA_AUDITS_COMPLETE')
stage('patient_scores','patient_scores.py','PATIENT_SCORE_ANALYSIS_COMPLETE')
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    jobs=[pool.submit(stage,'pseudobulk_DE','pseudobulk_DE.R','PSEUDOBULK_DE_COMPLETE'),pool.submit(stage,'score_sensitivities','score_sensitivities.py','SCORE_SENSITIVITY_COMPLETE')]
    for job in jobs:job.result()
stage('figures','figure_results.py','FIGURES_COMPLETE')
stage('result_summary','summarize_results.py','RESULT_SUMMARY_COMPLETE')
print('PATIENT_ANALYSIS_PIPELINE_COMPLETE')
