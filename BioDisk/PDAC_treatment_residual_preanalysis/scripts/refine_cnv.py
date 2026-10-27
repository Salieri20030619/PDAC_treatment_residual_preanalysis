from pathlib import Path
import sys,json,numpy as np,pandas as pd,anndata as ad,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
for done in sorted((ROOT/'03_processed').glob('GSE205013_P*.complete.json')):
    sample=done.name.split('_')[1].split('.')[0]
    if len(sys.argv)>1 and sample not in sys.argv[1:]:continue
    final=ROOT/'03_processed'/f'GSE205013_{sample}_cnv_refined.json'
    if final.exists():continue
    path=ROOT/'03_processed'/f'GSE205013_{sample}_annotated.h5ad'
    a=ad.read_h5ad(path);before=(a.obs.malignancy=='Malignant_CNV').sum();a.obs['malignancy']=a.obs.malignancy.astype(object)
    cnpath=ROOT/'03_processed'/f'GSE205013_{sample}_cnv.h5ad'
    if cnpath.exists():
        cn=ad.read_h5ad(cnpath);C=cn.obsm['X_cnv'].toarray()
        epi=cn.obs.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy()
        held=(cn.obs.cnv_reference=='Observation').to_numpy()&~epi
        coords=sorted(cn.uns['cnv']['chr_pos'].items(),key=lambda x:x[1]);rms=[];signed=[];names=[]
        for i,(chrom,start) in enumerate(coords):
            end=coords[i+1][1] if i+1<len(coords) else C.shape[1]
            if end-start<3:continue
            rms.append(np.sqrt(np.mean(C[:,start:end]**2,axis=1)));signed.append(np.mean(C[:,start:end],axis=1));names.append(chrom)
        rms=np.column_stack(rms);chrthreshold=np.quantile(rms[held],.99,axis=0)
        nchr=(rms>np.maximum(chrthreshold,1e-6)).sum(1)
        audit=json.loads((ROOT/'04_qc'/f'{sample}_qc_audit.json').read_text())
        score=cn.obs.cnv_score.to_numpy();threshold=audit['cnv']['score_threshold']
        positive=(score>threshold)&(nchr>=2);candidate=positive&epi;coherent=np.zeros(cn.n_obs,bool)
        for cluster in cn.obs.cluster.unique():
            ids=np.flatnonzero((cn.obs.cluster==cluster).to_numpy()&candidate)
            if len(ids)<10:continue
            profile=np.median(C[ids],axis=0);cosine=(C[ids]@profile)/np.maximum(np.linalg.norm(C[ids],axis=1)*np.linalg.norm(profile),1e-12)
            coherent[ids]=cosine>.3
        a.obs.loc[a.obs.major_type=='Epi','malignancy']='Epi_uncertain'
        a.obs.loc[a.obs.major_type=='Acinar','malignancy']='Acinar_like'
        a.obs.loc[a.obs.major_type=='Endocrine','malignancy']='Endocrine_like'
        low=(score<=np.quantile(score[held],.95))&(cn.obs.major_type=='Epi').to_numpy()
        a.obs.loc[cn.obs_names[low],'malignancy']='Ductal_CNVlow'
        a.obs.loc[cn.obs_names[candidate&coherent],'malignancy']='Malignant_CNV'
        a.obs.loc[cn.obs_names,'cnv_nchr']=nchr;a.obs.loc[cn.obs_names,'cnv_positive']=positive
        a.obs['cnv_coherent']=False;a.obs.loc[cn.obs_names,'cnv_coherent']=coherent
        pd.DataFrame(np.column_stack(signed),index=cn.obs_names,columns=names).to_csv(ROOT/'04_qc'/f'{sample}_chromosome_cnv.tsv.gz',sep='\t')
        pd.DataFrame(rms,index=cn.obs_names,columns=names).to_csv(ROOT/'04_qc'/f'{sample}_chromosome_cnv_rms.tsv.gz',sep='\t')
        audit['cnv']['chromosome_statistic']='RMS of 101-gene-window CNV; avoids median zeroing focal/arm events'
        audit['cnv']['chromosome_thresholds']=dict(zip(names,chrthreshold.tolist()))
        audit['cnv']['original_median_based_malignant']=int(before)
        audit['cnv']['n_malignant']=int((candidate&coherent).sum())
        audit['cnv']['heldout_positive_fraction']=float(positive[held].mean())
        audit['malignancy_counts']=a.obs.malignancy.value_counts().to_dict()
        (ROOT/'04_qc'/f'{sample}_qc_audit.json').write_text(json.dumps(audit,indent=2))
        cn.obs['cnv_nchr_initial']=cn.obs.cnv_nchr
        cn.obs['cnv_positive_initial']=cn.obs.cnv_positive
        cn.obs['cnv_nchr']=nchr;cn.obs['cnv_positive']=positive;cn.obs['cnv_coherent']=coherent
        cn.obs['malignancy']=a.obs.loc[cn.obs_names,'malignancy'].astype(str).to_numpy()
        cn.write_h5ad(cnpath,compression='gzip')
        del cn,C
    pb={};norm={};numbers={};strict=a.obs.qc_strict.to_numpy();epi=a.obs.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy();mal=(a.obs.malignancy=='Malignant_CNV').to_numpy()
    masks={'Malignant':mal,'EpiAll':epi,'MalStrictQC':mal&strict,'EpiStrictQC':epi&strict}
    masks.update({'Type_'+t:(a.obs.major_type==t).to_numpy() for t in a.obs.major_type.unique()})
    for group,mask in masks.items():
        numbers[group]=int(mask.sum());pb[group]=np.asarray(a.layers['counts'][mask].sum(0,dtype=np.int64)).ravel()
        norm[group]=np.asarray(a.X[mask].mean(0)).ravel() if mask.any() else np.full(a.n_vars,np.nan)
    np.savez_compressed(ROOT/'03_processed'/f'GSE205013_{sample}_pseudobulk.npz',genes=a.var_names.to_numpy(dtype=str),**pb)
    np.savez_compressed(ROOT/'03_processed'/f'GSE205013_{sample}_meanlog.npz',genes=a.var_names.to_numpy(dtype=str),**norm)
    (ROOT/'03_processed'/f'GSE205013_{sample}_cellnumbers.json').write_text(json.dumps(numbers,indent=2))
    allobs=pd.read_csv(ROOT/'04_qc'/f'{sample}_all_cell_annotations.tsv.gz',sep='\t',index_col=0)
    for col in ['malignancy','cnv_score','cnv_nchr','cnv_positive','cnv_coherent']:
        if col in a.obs:allobs[col]=a.obs[col].reindex(allobs.index)
    allobs.to_csv(ROOT/'04_qc'/f'{sample}_all_cell_annotations.tsv.gz',sep='\t',index_label='cell_id')
    a.write_h5ad(path,compression='gzip')
    rec={'sample':sample,'n_malignant_before':int(before),'n_malignant_after':int(mal.sum()),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'Technical correction before patient treatment testing: per-chromosome RMS preserves focal/arm changes that median of denoised CNV incorrectly zeros','status':'CNV_REFINEMENT_COMPLETE'}
    final.write_text(json.dumps(rec,indent=2))
    checkpoint=json.loads(done.read_text());checkpoint['n_malignant_initial']=checkpoint.get('n_malignant_initial',checkpoint['n_malignant']);checkpoint['n_malignant']=int(mal.sum());checkpoint['cnv_refined']=True;done.write_text(json.dumps(checkpoint,indent=2))
    print(json.dumps(rec),flush=True)
print('CNV_REFINEMENT_BATCH_COMPLETE')
