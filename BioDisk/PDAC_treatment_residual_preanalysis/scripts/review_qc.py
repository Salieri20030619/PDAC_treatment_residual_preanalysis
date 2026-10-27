from pathlib import Path
import json,pandas as pd,numpy as np
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
rows=[];ann=[];malrows=[]
for p in sorted((ROOT/'04_qc').glob('P*_qc_audit.json')):
    q=json.loads(p.read_text());sample=q['sample']
    if 'n_singlets' not in q:continue
    cells=json.loads((ROOT/'03_processed'/f'GSE205013_{sample}_cellnumbers.json').read_text())
    rows.append({'dataset':'GSE205013','patient_id':sample,'input_cells':q['n_input'],'QC_cells':q['n_qc'],'singlet_cells':q['n_singlets'],'doublet_fraction':q['doublet_rate'],'median_genes':q['ngenes']['median'],'median_UMI':q['UMI']['median'],'median_mt':q['mt']['median'],'malignant_cells':cells['Malignant'],'epithelial_cells':cells['EpiAll'],'strict_malignant_cells':cells['MalStrictQC'],'cnv_heldout_falsepositive':q.get('cnv',{}).get('heldout_positive_fraction',np.nan),'identity_final':(ROOT/'03_processed'/f'GSE205013_{sample}_identity_final.json').exists()})
    for key,n in q['cell_types'].items():ann.append({'dataset':'GSE205013','patient_id':sample,'cell_type':key,'n_cells':n})
    for key,n in q['malignancy_counts'].items():malrows.append({'dataset':'GSE205013','patient_id':sample,'malignancy':key,'n_cells':n})
df=pd.DataFrame(rows);df.to_csv(ROOT/'06_tables/GSE205013_sample_QC.tsv',sep='\t',index=False)
pd.DataFrame(ann).to_csv(ROOT/'06_tables/GSE205013_major_cell_types.tsv',sep='\t',index=False)
pd.DataFrame(malrows).to_csv(ROOT/'06_tables/GSE205013_malignancy_counts.tsv',sep='\t',index=False)
print(df.to_string(index=False));print('QC_SUMMARY_COMPLETE')
