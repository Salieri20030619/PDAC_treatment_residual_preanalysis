from pathlib import Path
import json,hashlib,datetime,numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
checks=[]
def require(name,condition,detail=''):
    checks.append({'check':name,'pass':bool(condition),'detail':detail})
    if not condition:raise AssertionError(name+': '+str(detail))
identity=list((ROOT/'03_processed').glob('GSE205013_P*_identity_final.json'))
require('27_final_sample_identity_checkpoints',len(identity)==27,str(len(identity)))
for marker in ['DOWNLOAD_AUDIT_PASS','METADATA_AUDIT_PASS']:
    path=ROOT/'logs'/f'{marker}.json';require(marker,path.exists() and json.loads(path.read_text())['status']==marker)
integrity=pd.read_csv(ROOT/'04_qc/pseudobulk_count_integrity.tsv',sep='\t')
require('raw_count_aggregation_audited_for_27_patients',integrity.patient_id.nunique()==27)
scores=pd.read_csv(ROOT/'06_tables/patient_mechanism_scores.tsv',sep='\t')
require('one_score_row_per_patient_and_compartment',not scores.duplicated(['dataset','patient_id','cell_group']).any())
require('two_cohorts_in_patient_scores',set(scores.dataset)=={'GSE205013','GSE202051'})
stats=pd.read_csv(ROOT/'06_tables/patient_mechanism_statistics.tsv',sep='\t')
for key,group in stats.groupby(['dataset','cell_group','comparison']):
    valid=group.pvalue.notna();p=group.loc[valid,'pvalue'].to_numpy()
    if len(p):
        require('FDR_family_'+str(key),np.allclose(group.loc[valid,'FDR'],multipletests(p,method='fdr_bh')[1],rtol=1e-8,atol=1e-12))
        require('valid_patient_Pvalues_'+str(key),np.all((p>=0)&(p<=1)))
    small=(group.n_untreated<10)|(group.n_treated<10)
    require('small_N_exploratory_'+str(key),(group.loc[small,'analysis_class']=='Exploratory').all())
loo=pd.read_csv(ROOT/'06_tables/leave_one_patient_out.tsv',sep='\t')
loo_groups={key:group for key,group in loo.groupby(['dataset','cell_group','comparison','geneset'])}
for _,row in stats[stats.status=='OK'].iterrows():
    d=loo_groups[(row.dataset,row.cell_group,row.comparison,row.geneset)]
    require('LOO_patient_count',len(d)==row.n_untreated+row.n_treated,f'{row.dataset}/{row.cell_group}/{row.comparison}/{row.geneset}')
    require('LOO_unique_omitted_patients',d.omitted_patient.is_unique)
    finite=d.effect.notna()
    require('LOO_estimable_count',int(finite.sum())==row.loo_n_estimable)
    if finite.any():require('LOO_direction_fraction',np.isclose((np.sign(d.loc[finite,'effect'])==np.sign(row.effect)).mean(),row.loo_same_direction))
de_files=list((ROOT/'06_tables').glob('*_edgeR.tsv.gz'));require('raw_count_gene_DE_present',len(de_files)>0)
for path in de_files:
    d=pd.read_csv(path,sep='\t');require('DE_gene_uniqueness_'+path.name,d.gene.is_unique)
    require('DE_FDR_'+path.name,np.allclose(d.FDR,multipletests(d.PValue,method='fdr_bh')[1],rtol=1e-8,atol=1e-12))
q1=pd.read_csv(ROOT/'06_tables/GSE205013_sample_QC.tsv',sep='\t')
require('QC_27_unique_patients',len(q1)==27 and q1.patient_id.is_unique)
require('QC_cell_accounting',(q1.input_cells>=q1.QC_cells).all() and (q1.QC_cells>=q1.singlet_cells).all() and (q1.singlet_cells>=q1.malignant_cells).all())
visual=json.loads((ROOT/'07_figures/visual_QA.json').read_text());require('figure_visual_review',visual['reviewed'] is True and visual['issues_remaining']==[])
review=json.loads((ROOT/'05_analysis/scientific_interpretation.json').read_text());require('scientific_result_review',review['reviewed'] is True)
for file in ['08_reports/PDAC_preanalysis_zh.md','08_reports/PDAC_preanalysis_zh.html','08_reports/download_integrity_report.md','06_tables/cross_cohort_direction.tsv','06_tables/analysis_patient_inclusion.tsv','01_metadata/patient_sample_metadata.tsv']:
    p=ROOT/file;require('deliverable_'+file,p.exists() and p.stat().st_size>100)
for suffix in ['md','html']:
    p=ROOT/'08_reports'/('PDAC_client_research_proposal_zh.'+suffix)
    require('detailed_client_proposal_'+suffix,p.exists() and p.stat().st_size>10000)
    require('client_proposal_UTF8_'+suffix,'\ufffd' not in p.read_text())
htmltext=(ROOT/'08_reports/PDAC_client_research_proposal_zh.html').read_text()
require('client_proposal_embedded_figures',htmltext.count('src="data:image/png;base64,')==2)
require('client_proposal_real_tables','<table>' in htmltext and '<pre>' not in htmltext)
clientqa=json.loads((ROOT/'08_reports/client_report_QA.json').read_text())
require('client_report_content_and_visual_QA',clientqa['reviewed'] is True and clientqa['issues_remaining']==[])
bg=pd.read_csv(ROOT/'06_tables/background_adjusted_statistics.tsv',sep='\t')
require('both_cohorts_QC_covariate_sensitivity',set(bg[bg.covariate=='median_mt'].dataset)=={'GSE205013','GSE202051'})
for key,group in bg.groupby(['dataset','cell_group','covariate']):
    valid=group.pvalue.notna()
    require('background_FDR_'+str(key),np.allclose(group.loc[valid,'FDR'],multipletests(group.loc[valid,'pvalue'],method='fdr_bh')[1],rtol=1e-8,atol=1e-12))
tables=[]
for folder in ['01_metadata','06_tables','07_figures','08_reports']:
    for p in sorted((ROOT/folder).iterdir()):
        if not p.is_file():continue
        if p.suffix in ['.tsv','.json','.md','.html','.png'] or p.name.endswith('.tsv.gz'):
            h=hashlib.sha256()
            with p.open('rb') as f:
                for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
            tables.append({'relative_path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
pd.DataFrame(tables).to_csv(ROOT/'08_reports/deliverable_manifest.tsv',sep='\t',index=False)
(ROOT/'logs/final_validation_checks.json').write_text(json.dumps(checks,indent=2))
record={'status':'PREANALYSIS_COMPLETE','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_checks':len(checks),'scope':'Two-cohort public transcriptomic preanalysis, patient-level statistics, sensitivity, reports and figures','biological_claim':'See reviewed Chinese report; completion does not imply a positive proteostasis result','limitations_retained':review['pending_items']}
(ROOT/'logs/PREANALYSIS_COMPLETE.json').write_text(json.dumps(record,indent=2))
print('PREANALYSIS_COMPLETE')
