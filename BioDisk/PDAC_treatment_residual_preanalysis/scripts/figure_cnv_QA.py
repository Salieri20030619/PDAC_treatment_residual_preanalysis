from pathlib import Path
import numpy as np,pandas as pd,json,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
meta=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
means=[];rms=[];names=[];qrows=[]
for file in sorted((ROOT/'04_qc').glob('P*_chromosome_cnv.tsv.gz')):
    sample=file.name.split('_')[0]
    if not (ROOT/'03_processed'/f'GSE205013_{sample}_identity_final.json').exists():continue
    q=json.loads((ROOT/'04_qc'/f'{sample}_qc_audit.json').read_text())
    if 'score_threshold' not in q['cnv']:continue
    obs=pd.read_csv(ROOT/'04_qc'/f'{sample}_all_cell_annotations.tsv.gz',sep='\t',index_col=0)
    chrom=pd.read_csv(file,sep='\t',index_col=0)
    chrms=pd.read_csv(ROOT/'04_qc'/f'{sample}_chromosome_cnv_rms.tsv.gz',sep='\t',index_col=0)
    ix=chrom.index.intersection(obs.index[obs.malignancy=='Malignant_CNV'])
    if len(ix):
        means.append(chrom.loc[ix].mean());rms.append(chrms.loc[ix].mean());names.append(sample)
    q=json.loads((ROOT/'04_qc'/f'{sample}_qc_audit.json').read_text())
    qrows.append({'sample':sample,'n_malignant':len(ix),'holdout_pos_frac':q['cnv'].get('heldout_positive_fraction',np.nan),'score_threshold':q['cnv']['score_threshold']})
df=pd.DataFrame(means,index=names).reindex(columns=['chr'+str(i) for i in range(1,23)])
df.to_csv(ROOT/'04_qc/patient_mean_CNV.tsv',sep='\t')
pd.DataFrame(rms,index=names).to_csv(ROOT/'04_qc/patient_mean_CNV_RMS.tsv',sep='\t')
pd.DataFrame(qrows).to_csv(ROOT/'04_qc/CNV_holdout_calibration.tsv',sep='\t',index=False)
fig,ax=plt.subplots(figsize=(11,7));lim=np.nanquantile(np.abs(df.to_numpy()),.95)
im=ax.imshow(df.to_numpy(),cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto');ax.set_yticks(np.arange(len(df)),df.index);ax.set_xticks(np.arange(22),[str(i) for i in range(1,23)]);ax.set_xlabel('Chr.');ax.set_ylabel('Patient');fig.colorbar(im,ax=ax,label='Mean CNV');fig.tight_layout();fig.savefig(ROOT/'07_figures/05_patient_CNV.png',dpi=220);plt.close(fig)
print('CNV_QA_FIGURE_COMPLETE')
