from pathlib import Path
import pandas as pd,numpy as np,json
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
meta=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
qc=pd.read_csv(ROOT/'06_tables/GSE205013_sample_QC.tsv',sep='\t')
df=meta[meta.dataset=='GSE205013'].merge(qc,on=['dataset','patient_id'],validate='one_to_one')
df['eligible_malignant']=df.malignant_cells>=30
df.to_csv(ROOT/'01_metadata/GSE205013_final_analysis_eligibility.tsv',sep='\t',index=False)
for field in ['procedure','stage','regimen_group']:
    cross=df.groupby(['comparison_scope','eligible_malignant',field,'treatment_status']).size().reset_index(name='n_patients')
    cross.to_csv(ROOT/'06_tables'/f'GSE205013_postQC_{field}_crosstab.tsv',sep='\t',index=False)
primary=df[(df.comparison_scope=='Primary')&df.eligible_malignant]
print('PRIMARY_ELIGIBLE_PATIENTS',primary[['patient_id','treatment_status','procedure','stage','malignant_cells']].to_string(index=False))
print('PRIMARY_GROUP_COUNTS',primary.groupby('treatment_status').size().to_dict())
print('PROCEDURE_COUNTS',primary.groupby(['procedure','treatment_status']).size().to_dict())
print('EXCLUDED_PRIMARY',df[(df.comparison_scope=='Primary')&~df.eligible_malignant][['patient_id','treatment_status','malignant_cells']].to_string(index=False))
print('TOTAL_COUNTS',qc[['input_cells','QC_cells','singlet_cells','malignant_cells']].sum().to_dict())
print('DISCOVERY_COMPARABILITY_COMPLETE')
