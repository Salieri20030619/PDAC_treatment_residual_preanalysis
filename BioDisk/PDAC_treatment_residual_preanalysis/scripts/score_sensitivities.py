from pathlib import Path
import json,numpy as np,pandas as pd,statsmodels.api as sm
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
scores=pd.read_csv(ROOT/'06_tables/patient_mechanism_scores.tsv',sep='\t')
sets=json.loads((ROOT/'00_docs/genesets_frozen.json').read_text())
meta=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
qc=[]
for p in (ROOT/'04_qc').glob('P*_all_cell_annotations.tsv.gz'):
    sample=p.name.split('_')[0];df=pd.read_csv(p,sep='\t')
    for label,mask in [('Malignant',df.malignancy=='Malignant_CNV'),('EpiAll',df.major_type.isin(['Epi','Acinar','Endocrine']))]:
        x=df.loc[mask]
        if len(x):qc.append(dict(dataset='GSE205013',patient_id=sample,cell_group=label,median_mt=x.pct_mt.median(),median_logUMI=np.log1p(x.UMI).median(),median_genes=x.n_genes.median(),median_cnv=x.cnv_score.median()))
v=pd.read_csv(ROOT/'04_qc/GSE202051_all_cell_annotations.tsv.gz',sep='\t',usecols=['patient_id','malignancy','major_type','pct_mt','UMI','n_genes'])
for patient,df in v.groupby('patient_id'):
    for label,mask in [('Malignant',df.malignancy=='Malignant_published_CNV'),('EpiAll',df.major_type.isin(['Epi','Acinar','Endocrine']))]:
        x=df.loc[mask]
        if len(x):qc.append(dict(dataset='GSE202051',patient_id=patient,cell_group=label,median_mt=x.pct_mt.median(),median_logUMI=np.log1p(x.UMI).median(),median_genes=x.n_genes.median(),median_cnv=np.nan))
pd.DataFrame(qc).to_csv(ROOT/'06_tables/patient_cellgroup_QC.tsv',sep='\t',index=False)
s=scores.merge(pd.DataFrame(qc),on=['dataset','patient_id','cell_group'],how='left')
out=[];corr=[]
background=['HALLMARK_HYPOXIA','HALLMARK_E2F_TARGETS','HALLMARK_APOPTOSIS','MOFFITT_BASAL','REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS','median_mt','median_logUMI']
for ds in s.dataset.unique():
    for group in ['Malignant','EpiAll']:
        d=s[(s.dataset==ds)&(s.cell_group==group)&(s.comparison_scope=='Primary')&s.analysis_eligible].copy()
        if ds=='GSE202051':d=d[d.regimen_group.isin(['Untreated','CRT'])]
        d=d.reset_index(drop=True)
        for cov in background:
            if cov not in d or not d[cov].notna().all() or d[cov].nunique()<2:continue
            x=pd.DataFrame({'Intercept':np.ones(len(d)),'Treated':(d.treatment_status=='Treated').astype(float)})
            if ds=='GSE205013' and d.procedure.nunique()>1:x=pd.concat([x,pd.get_dummies(d.procedure,prefix='Proc',drop_first=True,dtype=float)],axis=1)
            x['covariate']=(d[cov]-d[cov].mean())/d[cov].std()
            if np.linalg.matrix_rank(x)<x.shape[1] or len(d)-x.shape[1]<3:continue
            for gs in sets:
                if gs not in d or gs==cov:continue
                f=sm.OLS(d[gs],x).fit(cov_type='HC3',use_t=True)
                out.append(dict(dataset=ds,cell_group=group,covariate=cov,geneset=gs,n_patients=len(d),effect=f.params['Treated'],pvalue=f.pvalues['Treated']))
                rho,p=spearmanr(d[gs],d[cov]);corr.append(dict(dataset=ds,cell_group=group,covariate=cov,geneset=gs,n_patients=len(d),rho=rho,pvalue=p))
for items,file in [(out,'background_adjusted_statistics.tsv'),(corr,'patient_score_QC_correlations.tsv')]:
    df=pd.DataFrame(items)
    if not df.empty:
        df['FDR']=np.nan
        for _,idx in df.groupby(['dataset','cell_group','covariate']).groups.items():
            ok=df.loc[idx,'pvalue'].notna();use=df.loc[idx].index[ok]
            df.loc[use,'FDR']=multipletests(df.loc[use,'pvalue'],method='fdr_bh')[1]
    df.to_csv(ROOT/'06_tables'/file,sep='\t',index=False)
tests=pd.read_csv(ROOT/'06_tables/patient_mechanism_statistics.tsv',sep='\t')
a=tests[(tests.dataset=='GSE205013')&(tests.cell_group=='Malignant')&(tests.comparison=='Primary')].copy()
b=tests[(tests.dataset=='GSE202051')&(tests.cell_group=='Malignant')&(tests.comparison=='CRT')].copy()
cols=['geneset','effect','pvalue','FDR','loo_same_direction','loo_estimable_fraction','n_untreated','n_treated']
cross=a[cols].merge(b[cols],on='geneset',suffixes=('_discovery','_validation'),how='outer')
cross['same_direction']=np.sign(cross.effect_discovery)==np.sign(cross.effect_validation)
for name,group,comparison in [('Resect','Malignant','Resection'),('NoRT','Malignant','NoRT'),('StrictQC','MalStrictQC','Primary'),('EpiAll','EpiAll','Primary')]:
    sub=tests[(tests.dataset=='GSE205013')&(tests.cell_group==group)&(tests.comparison==comparison)].set_index('geneset')
    cross[name+'_effect']=cross.geneset.map(sub.effect)
    cross[name+'_same_direction']=np.sign(cross[name+'_effect'])==np.sign(cross.effect_discovery)
vstrict=tests[(tests.dataset=='GSE202051')&(tests.cell_group=='MalStrictQC')&(tests.comparison=='CRT')].set_index('geneset')
cross['ValidationStrictQC_effect']=cross.geneset.map(vstrict.effect)
cross['ValidationStrictQC_same_direction']=np.sign(cross.ValidationStrictQC_effect)==np.sign(cross.effect_validation)
cross['signature_origin']=[sets[x]['category'] for x in cross.geneset]
cross['candidate_gate']=(cross.same_direction&(cross.loo_same_direction_discovery>=.9)&(cross.loo_same_direction_validation>=.9)&(cross.loo_estimable_fraction_discovery>=.9)&(cross.loo_estimable_fraction_validation>=.9)&cross.Resect_same_direction&cross.NoRT_same_direction&cross.StrictQC_same_direction&((cross.FDR_discovery<.1)|(cross.FDR_validation<.1)))
cross['independent_signature_validation']=cross.signature_origin!='Published_PDAC_program'
cross.to_csv(ROOT/'06_tables/cross_cohort_direction.tsv',sep='\t',index=False)
print('SCORE_SENSITIVITY_COMPLETE',len(cross))
