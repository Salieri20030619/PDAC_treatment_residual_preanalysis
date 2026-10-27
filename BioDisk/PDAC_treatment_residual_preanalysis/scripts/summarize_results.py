from pathlib import Path
import pandas as pd,json
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
inc=pd.read_csv(ROOT/'06_tables/analysis_patient_inclusion.tsv',sep='\t')
primary=inc[(inc.cell_group=='Malignant')&(inc.comparison_scope=='Primary')]
print('PRIMARY_INCLUSION');print(primary[['dataset','patient_id','treatment_status','procedure','regimen_group','n_cells','analysis_eligible']].to_string(index=False))
stats=pd.read_csv(ROOT/'06_tables/patient_mechanism_statistics.tsv',sep='\t');cross=pd.read_csv(ROOT/'06_tables/cross_cohort_direction.tsv',sep='\t')
columns=['geneset','effect_discovery','FDR_discovery','loo_same_direction_discovery','effect_validation','FDR_validation','loo_same_direction_validation','same_direction','Resect_same_direction','NoRT_same_direction','StrictQC_same_direction','candidate_gate']
print('CROSS_COHORT_TOP');print(cross.sort_values(['candidate_gate','same_direction','FDR_discovery'],ascending=[False,False,True])[columns].head(25).to_string(index=False))
genesets=json.loads((ROOT/'00_docs/genesets_frozen.json').read_text())
proteo=[k for k,v in genesets.items() if v['category']=='Proteostasis_candidate']+['HALLMARK_UNFOLDED_PROTEIN_RESPONSE']
print('ALL_PROTEOSTASIS');print(cross[cross.geneset.isin(proteo)][columns].to_string(index=False))
print('DE_AUDIT');print(pd.read_csv(ROOT/'06_tables/pseudobulk_DE_audit.tsv',sep='\t').to_string(index=False))
anchors=['HSPA5','HSP90B1','HSP90AA1','HSP90AB1','XBP1','ATF4','ATF6','DDIT3','EIF2AK3','EIF2AK4','EIF2S1','DNAJB1','HSPA1A','HSPA1B','PSMB5','PSMA1','SQSTM1','MAP1LC3B','ATG5','LAMP1','LAMP2','GATA6','KRT17','S100A2','MKI67']
for ds,cmp in [('GSE205013','Primary'),('GSE202051','CRT')]:
    p=ROOT/'06_tables'/f'{ds}_Malignant_{cmp}_edgeR.tsv.gz'
    if p.exists():
        d=pd.read_csv(p,sep='\t');print('ANCHOR_GENES',ds);print(d[d.gene.isin(anchors)][['gene','logFC','PValue','FDR']].sort_values('gene').to_string(index=False))
print('RESULT_SUMMARY_COMPLETE')
