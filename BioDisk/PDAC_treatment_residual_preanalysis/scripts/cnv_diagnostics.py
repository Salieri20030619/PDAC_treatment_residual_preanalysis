from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,json
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
for sample in ['P01','P02','P03','P04','P05']:
    p=ROOT/'03_processed'/f'GSE205013_{sample}_cnv.h5ad'
    if not p.exists():continue
    a=ad.read_h5ad(p)
    obs=a.obs;X=a.obsm['X_cnv'].toarray()
    print('SAMPLE',sample,'TYPES',obs.major_type.value_counts().to_dict())
    print(obs.groupby('major_type',observed=True).agg(n=('cnv_score','size'),score50=('cnv_score','median'),score90=('cnv_score',lambda x:np.quantile(x,.9)),pos=('cnv_positive','sum'),nchr=('cnv_nchr','max')).to_string())
    print('CNV_SCORE_THRESHOLD',json.loads((ROOT/'04_qc'/f'{sample}_qc_audit.json').read_text())['cnv'])
    coords=sorted(a.uns['cnv']['chr_pos'].items(),key=lambda x:x[1]);cr=[]
    for i,(chrom,start) in enumerate(coords):
        end=coords[i+1][1] if i+1<len(coords) else X.shape[1]
        cr.append(np.sqrt(np.mean(X[:,start:end]**2,axis=1)))
    cr=np.column_stack(cr);held=(obs.cnv_reference=='Observation')&~obs.major_type.isin(['Epi','Acinar','Endocrine'])
    cut=np.quantile(cr[held.to_numpy()],.99,axis=0);n=(cr>np.maximum(cut,1e-6)).sum(1)
    scorecut=json.loads((ROOT/'04_qc'/f'{sample}_qc_audit.json').read_text())['cnv']['score_threshold']
    print('RMS_CHR_CANDIDATES',int(((n>=2)&(obs.cnv_score>scorecut)&obs.major_type.isin(['Epi','Acinar','Endocrine'])).sum()),'WITH_POSITIVE_SCORE',int(((obs.cnv_score>scorecut)&obs.major_type.isin(['Epi','Acinar','Endocrine'])).sum()))
print('CNV_DIAGNOSTIC_COMPLETE')
