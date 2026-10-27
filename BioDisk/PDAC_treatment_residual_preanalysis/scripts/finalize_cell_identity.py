from pathlib import Path
import sys,json,datetime,shutil,os
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
sys.path.insert(0,str(ROOT/'scripts/vendor'))
import numpy as np,pandas as pd,anndata as ad,infercnvpy as cnv
from process_205013 import PANELS,REF
corrections=pd.read_csv(ROOT/'01_metadata/annotation_corrections.tsv',sep='\t',dtype={'cluster':str})
archive=ROOT/'03_processed/cnv_previous';archive.mkdir(exist_ok=True)
obsarchive=ROOT/'04_qc/annotations_previous';obsarchive.mkdir(exist_ok=True)
pos=pd.read_csv(ROOT/'00_docs/hg38_refGene.txt.gz',sep='\t',header=None,usecols=[2,4,5,12],names=['chromosome','start','end','gene_symbol'])
pos=pos[pos.chromosome.isin(['chr'+str(i) for i in range(1,23)])].sort_values(['gene_symbol','chromosome','start']).drop_duplicates('gene_symbol').set_index('gene_symbol')
for sample in sys.argv[1:]:
    assert (ROOT/'03_processed'/f'GSE205013_{sample}.complete.json').exists(),f'Initial sample processing not complete: {sample}'
    final=ROOT/'03_processed'/f'GSE205013_{sample}_identity_final.json'
    if final.exists():print('SKIP_FINAL_IDENTITY',sample,flush=True);continue
    path=ROOT/'03_processed'/f'GSE205013_{sample}_annotated.h5ad';a=ad.read_h5ad(path)
    audit=json.loads((ROOT/'04_qc'/f'{sample}_qc_audit.json').read_text())
    previous=int((a.obs.malignancy=='Malignant_CNV').sum())
    a.obs['major_type']=a.obs.major_type.astype(object)
    for _,row in corrections[corrections['sample']==sample].iterrows():a.obs.loc[a.obs.cluster.astype(str)==row.cluster,'major_type']=row.major_type
    a.obs['malignancy']='Not_epithelial'
    for major,label in [('Epi','Epi_uncertain'),('Acinar','Acinar_like'),('Endocrine','Endocrine_like'),('Mixed','Mixed_uncertain')]:a.obs.loc[a.obs.major_type==major,'malignancy']=label
    idx={g:i for i,g in enumerate(a.var_names)}
    epigenes=[idx[g] for g in ['EPCAM','KRT8','KRT18','KRT19','KRT7'] if g in idx]
    epi_expr=np.asarray(a.X[:,epigenes].mean(1)).ravel()
    confident=np.zeros(a.n_obs,bool)
    for label in REF:
        cols=[idx[g] for g in PANELS[label] if g in idx];values=a.X[:,cols].toarray()
        evidence=(values.mean(1)>=.5)&((values>.3).sum(1)>=2)
        confident|=evidence&((a.obs.major_type==label).to_numpy()|(a.obs.major_type=='StromaMix').to_numpy())
    confident&=epi_expr<.5
    refs=np.flatnonzero(confident);epimask=a.obs.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy()
    a.obs['cnv_reference']='Observation';a.obs['cnv_score']=np.nan;a.obs['cnv_nchr']=0;a.obs['cnv_positive']=False;a.obs['cnv_coherent']=False
    if len(refs)>=100 and epimask.sum()>=10:
        rng=np.random.default_rng(20260922);pool=rng.choice(refs,size=min(2000,len(refs)),replace=False);train=pool[:len(pool)//2];hold=pool[len(pool)//2:]
        a.obs.iloc[train,a.obs.columns.get_loc('cnv_reference')]='Reference'
        for key in ['chromosome','start','end']:a.var[key]=a.var.gene_symbol.map(pos[key]).to_numpy()
        selected=np.unique(np.r_[np.flatnonzero(epimask),pool]);c=a[selected].copy()
        cnv.tl.infercnv(c,reference_key='cnv_reference',reference_cat='Reference',window_size=101,step=10,dynamic_threshold=1.5,n_jobs=2,chunksize=1000)
        C=c.obsm['X_cnv'].toarray();score=np.sqrt(np.mean(C*C,axis=1));held=c.obs_names.isin(a.obs_names[hold]);base=score[held]
        threshold=max(float(np.quantile(base,.99)),float(np.median(base)+5*np.median(np.abs(base-np.median(base)))))
        coords=sorted(c.uns['cnv']['chr_pos'].items(),key=lambda x:x[1]);rms=[];signed=[];names=[]
        for j,(chrom,start) in enumerate(coords):
            end=coords[j+1][1] if j+1<len(coords) else C.shape[1]
            if end-start<3:continue
            rms.append(np.sqrt(np.mean(C[:,start:end]**2,axis=1)));signed.append(np.mean(C[:,start:end],axis=1));names.append(chrom)
        rms=np.column_stack(rms);chrthr=np.quantile(rms[held],.99,axis=0);nchr=(rms>np.maximum(chrthr,1e-6)).sum(1)
        positive=(score>threshold)&(nchr>=2);epi=c.obs.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy();candidate=positive&epi;coherent=np.zeros(c.n_obs,bool)
        for cluster in c.obs.cluster.unique():
            ids=np.flatnonzero((c.obs.cluster==cluster).to_numpy()&candidate)
            if len(ids)<10:continue
            profile=np.median(C[ids],axis=0);cosine=C[ids]@profile/np.maximum(np.linalg.norm(C[ids],axis=1)*np.linalg.norm(profile),1e-12);coherent[ids]=cosine>.3
        low=(score<=np.quantile(base,.95))&(c.obs.major_type=='Epi').to_numpy()
        a.obs.loc[c.obs_names[low],'malignancy']='Ductal_CNVlow';a.obs.loc[c.obs_names[candidate&coherent],'malignancy']='Malignant_CNV'
        for key,values in [('cnv_score',score),('cnv_nchr',nchr),('cnv_positive',positive),('cnv_coherent',coherent)]:a.obs.loc[c.obs_names,key]=values
        for key in ['cnv_score','cnv_nchr','cnv_positive','cnv_coherent','malignancy']:c.obs[key]=a.obs.loc[c.obs_names,key]
        cnpath=ROOT/'03_processed'/f'GSE205013_{sample}_cnv.h5ad'
        if cnpath.exists() and not (archive/cnpath.name).exists():cnpath.rename(archive/cnpath.name)
        c.write_h5ad(cnpath,compression='gzip')
        pd.DataFrame(np.column_stack(signed),index=c.obs_names,columns=names).to_csv(ROOT/'04_qc'/f'{sample}_chromosome_cnv.tsv.gz',sep='\t')
        pd.DataFrame(rms,index=c.obs_names,columns=names).to_csv(ROOT/'04_qc'/f'{sample}_chromosome_cnv_rms.tsv.gz',sep='\t')
        audit['cnv']={'method':'infercnvpy 0.6.1;101-gene window;step10;chromosome RMS','reference_eligibility':'Major immune/stromal label AND lineage meanlog>=0.5 AND >=2 markers logexpr>0.3 AND mean epithelial logexpr<0.5; mixed epithelial/immune clusters excluded','n_reference_eligible':len(refs),'n_reference_train':len(train),'n_reference_holdout':len(hold),'score_threshold':threshold,'heldout_q99':float(np.quantile(base,.99)),'heldout_positive_fraction':float(positive[held].mean()),'n_chr_required':2,'coherence_cosine_gt':.3,'coherence_min_cells':10,'n_malignant':int((candidate&coherent).sum()),'chromosome_thresholds':dict(zip(names,chrthr.tolist()))}
        del c,C
    else:
        audit['cnv']={'status':'Insufficient unambiguous normal references or epithelial candidates; no malignancy assignment','n_reference_eligible':len(refs),'n_malignant':0}
        old=ROOT/'03_processed'/f'GSE205013_{sample}_cnv.h5ad'
        if old.exists():
            previous_path=archive/(old.name if not (archive/old.name).exists() else old.stem+'.pre_final.h5ad')
            old.rename(previous_path)
    pb={};norm={};nums={};strict=a.obs.qc_strict.to_numpy();mal=(a.obs.malignancy=='Malignant_CNV').to_numpy()
    masks={'Malignant':mal,'EpiAll':epimask,'MalStrictQC':mal&strict,'EpiStrictQC':epimask&strict}
    masks.update({'Type_'+t:(a.obs.major_type==t).to_numpy() for t in a.obs.major_type.unique()})
    for group,mask in masks.items():
        nums[group]=int(mask.sum());pb[group]=np.asarray(a.layers['counts'][mask].sum(0,dtype=np.int64)).ravel();norm[group]=np.asarray(a.X[mask].mean(0)).ravel() if mask.any() else np.full(a.n_vars,np.nan)
    np.savez_compressed(ROOT/'03_processed'/f'GSE205013_{sample}_pseudobulk.npz',genes=a.var_names.to_numpy(dtype=str),**pb)
    np.savez_compressed(ROOT/'03_processed'/f'GSE205013_{sample}_meanlog.npz',genes=a.var_names.to_numpy(dtype=str),**norm)
    (ROOT/'03_processed'/f'GSE205013_{sample}_cellnumbers.json').write_text(json.dumps(nums,indent=2))
    apath=ROOT/'04_qc'/f'{sample}_all_cell_annotations.tsv.gz'
    if not (obsarchive/apath.name).exists():shutil.copy2(apath,obsarchive/apath.name)
    allobs=pd.read_csv(apath,sep='\t',index_col=0)
    for col in ['major_type','malignancy','cnv_score','cnv_nchr','cnv_positive','cnv_coherent']:allobs[col]=a.obs[col].reindex(allobs.index)
    allobs.to_csv(apath,sep='\t',index_label='cell_id')
    audit['cell_types']=a.obs.major_type.value_counts().to_dict();audit['malignancy_counts']=a.obs.malignancy.value_counts().to_dict()
    (ROOT/'04_qc'/f'{sample}_qc_audit.json').write_text(json.dumps(audit,indent=2))
    temp=path.with_suffix('.tmp.h5ad');a.write_h5ad(temp,compression='gzip');os.replace(temp,path)
    rec={'sample':sample,'n_malignant_previous':previous,'n_malignant_final':int(mal.sum()),'normal_references_eligible':len(refs),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'CELL_IDENTITY_FINAL_COMPLETE'}
    final.write_text(json.dumps(rec,indent=2))
    done=ROOT/'03_processed'/f'GSE205013_{sample}.complete.json';checkpoint=json.loads(done.read_text());checkpoint['n_malignant']=int(mal.sum());checkpoint['identity_final']=True;done.write_text(json.dumps(checkpoint,indent=2))
    annpath=ROOT/'04_qc'/f'{sample}_cluster_annotation.tsv';ann=pd.read_csv(annpath,sep='\t',dtype={'cluster':str});ann['initial_major_type']=ann.major_type
    ann['major_type']=ann.cluster.map(a.obs[['cluster','major_type']].drop_duplicates().assign(cluster=lambda d:d.cluster.astype(str)).set_index('cluster').major_type)
    ann.to_csv(annpath,sep='\t',index=False)
    print(json.dumps(rec),flush=True)
print('FINAL_CELL_IDENTITY_BATCH_COMPLETE')
