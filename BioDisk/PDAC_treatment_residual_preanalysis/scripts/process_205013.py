from pathlib import Path
import os,sys,json,gzip,traceback,datetime
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('NUMBA_NUM_THREADS','4')
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
sys.path.insert(0,str(ROOT/'scripts/vendor'))
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
PANELS={
 'Tcell':['CD3D','CD3E','TRAC','CD247','CD2'],
 'NK':['NKG7','GNLY','KLRD1','PRF1','KLRF1'],
 'Bcell':['CD79A','MS4A1','CD79B','CD37','BANK1'],
 'Plasma':['MZB1','JCHAIN','SDC1','IGHG1','TNFRSF17'],
 'Myeloid':['LST1','TYROBP','FCER1G','CTSS','AIF1','CD68'],
 'Mast':['TPSAB1','TPSB2','CPA3','KIT','MS4A2'],
 'CAF':['COL1A1','COL1A2','DCN','LUM','COL3A1','PDGFRA'],
 'Endo':['PECAM1','VWF','EMCN','KDR','RAMP2'],
 'Pericyte':['RGS5','CSPG4','MCAM','PDGFRB','ACTA2'],
 'Schwann':['SOX10','PLP1','S100B','MPZ','S100A1'],
 'Hepato':['ALB','APOA1','TTR','APOA2','ASGR1'],
 'Acinar':['PRSS1','PRSS2','CPA1','CTRB1','CELA3A'],
 'Endocrine':['CHGA','CHGB','INS','IAPP','GCG','SST'],
 'Epi':['EPCAM','KRT8','KRT18','KRT19','KRT7','MUC1','TACSTD2']}
REF=['Tcell','NK','Bcell','Myeloid','CAF','Endo','Pericyte']
def stats(v):return dict(zip(['min','p01','p10','median','p90','p99','max'],map(float,np.quantile(v,[0,.01,.1,.5,.9,.99,1]))))
def savejson(path,obj):path.write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def run(sample):
    final=ROOT/'03_processed'/f'GSE205013_{sample}.complete.json'
    if final.exists():print('SKIP_SAMPLE_COMPLETE',sample,flush=True);return
    meta=pd.read_csv(ROOT/'01_metadata/sample_metadata_GEO.tsv',sep='\t',keep_default_na=False)
    row=meta[(meta.dataset=='GSE205013')&(meta.sample_id==sample)].iloc[0]
    base=ROOT/'02_raw/GSE205013/matrices'
    prefix=f'{row.GSM}_{sample}_'
    download_record=ROOT/'logs/GSE205013_RAW.tar.download.json'
    archive_ready=download_record.exists() and json.loads(download_record.read_text()).get('status')=='INTEGRITY_PASS'
    if not archive_ready or not (base/(prefix+'matrix.mtx.gz')).exists():base=ROOT/'02_raw/GSE205013/early_matrices'
    for suffix in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']:assert (base/(prefix+suffix)).exists()
    print('START_SAMPLE',sample,datetime.datetime.now().isoformat(),flush=True)
    with threadpool_limits(limits=4),gzip.open(base/(prefix+'matrix.mtx.gz')) as f:x=mmread(f).T.tocsr()
    x.sum_duplicates();x.eliminate_zeros()
    assert x.data.min()>=0 and np.all(np.isfinite(x.data)) and np.all(x.data==np.floor(x.data))
    feat=pd.read_csv(base/(prefix+'features.tsv.gz'),sep='\t',header=None)
    bc=pd.read_csv(base/(prefix+'barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str)
    assert x.shape==(len(bc),len(feat)) and not bc.duplicated().any()
    var=pd.DataFrame({'gene_id':feat[0].astype(str).to_numpy(),'gene_symbol':feat[1].astype(str).to_numpy()},index=feat[1].astype(str))
    var.index.name='gene'
    a=ad.AnnData(x.astype(np.float32),obs=pd.DataFrame(index=(sample+'_'+bc).to_numpy()),var=var);a.var_names_make_unique();a.obs.index.name='cell_id'
    a.obs['sample_id']=sample;a.obs['GSM']=row.GSM
    genes=a.var.gene_symbol.to_numpy();umi=np.asarray(a.X.sum(1)).ravel();ng=a.X.getnnz(1)
    mt=np.asarray(a.X[:,np.char.startswith(genes.astype(str),'MT-')].sum(1)).ravel()/np.maximum(umi,1)*100
    hb=np.asarray(a.X[:,np.isin(genes,['HBA1','HBA2','HBB','HBM','ALAS2'])].sum(1)).ravel()/np.maximum(umi,1)*100
    a.obs['n_genes']=ng;a.obs['UMI']=umi;a.obs['pct_mt']=mt;a.obs['pct_hb']=hb
    # Fixed published floor, applied after recording each sample's distribution.
    a.obs['qc_pass']=(ng>500)&(umi>1500)&(mt<15)&(hb<1)
    a.obs['qc_strict']=(ng>750)&(umi>2500)&(mt<10)&(hb<1)
    audit={'sample':sample,'n_input':a.n_obs,'n_genes':a.n_vars,'counts_integer':True,'ngenes':stats(ng),'UMI':stats(umi),'mt':stats(mt),'hb':stats(hb),'thresholds':{'genes_gt':500,'UMI_gt':1500,'mt_lt':15,'hb_lt':1},'threshold_basis':'Werba2023 Methods; first-sample actual distributions confirmed unfiltered low-quality tail; per-sample quantiles recorded before application','doublet_input_status':'GEO explicitly says filtered CellRanger, not final SoupX/scDblFinder object','n_qc':int(a.obs.qc_pass.sum())}
    savejson(ROOT/'04_qc'/f'{sample}_qc_audit.json',audit)
    all_obs=a.obs.copy();a=a[a.obs.qc_pass].copy()
    if a.n_obs<100:raise ValueError('Insufficient cells for sample annotation')
    # Once per CellRanger sample; no repeat filtering of published final objects.
    sc.pp.scrublet(a,expected_doublet_rate=min(.15,max(.02,.008*a.n_obs/1000)),n_prin_comps=30,random_state=20260922,verbose=True)
    audit['doublet_rate']=float(a.obs.predicted_doublet.mean())
    all_obs['doublet_score']=a.obs.doublet_score.reindex(all_obs.index)
    all_obs['predicted_doublet']=a.obs.predicted_doublet.reindex(all_obs.index).astype('boolean')
    a=a[~a.obs.predicted_doublet].copy();a.layers['counts']=a.X.copy()
    sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a,n_top_genes=2500,flavor='seurat')
    h=a[:,a.var.highly_variable].copy();sc.pp.scale(h,max_value=10)
    sc.tl.pca(h,n_comps=min(40,h.n_obs-2),svd_solver='arpack',random_state=20260922)
    sc.pp.neighbors(h,n_neighbors=20,n_pcs=min(30,h.obsm['X_pca'].shape[1]),random_state=20260922)
    sc.tl.leiden(h,resolution=.8,key_added='cluster',random_state=20260922,flavor='igraph',n_iterations=2,directed=False)
    a.obs['cluster']=h.obs.cluster;a.obsm['X_pca']=h.obsm['X_pca'];del h
    labels=[];means=[];clusters=list(a.obs.cluster.cat.categories)
    for c in clusters:means.append(np.asarray(a.X[(a.obs.cluster==c).to_numpy()].mean(0)).ravel())
    means=np.array(means);gene_index={g:i for i,g in enumerate(a.var_names)}
    z=(means-means.mean(0))/np.maximum(means.std(0),.15)
    panel_scores=np.column_stack([z[:,[gene_index[g] for g in gs if g in gene_index]].mean(1) for gs in PANELS.values()])
    for ci,c in enumerate(clusters):
        scores=panel_scores[ci];best=int(np.argmax(scores));label=list(PANELS)[best]
        inds=[gene_index[g] for g in PANELS[label] if g in gene_index]
        if (means[ci,inds]>.25).sum()<2 or scores[best]<.25:label='Unknown'
        # Distinguish T/NK by lineage-specific receptor genes, not cytotoxicity alone.
        tinds=[gene_index[g] for g in ['CD3D','CD3E','TRAC'] if g in gene_index]
        if label=='NK' and means[ci,tinds].mean()>.5:label='Tcell'
        top=np.argsort(means[ci]-np.mean(np.delete(means,ci,axis=0),axis=0))[-20:][::-1] if len(clusters)>1 else np.argsort(means[ci])[-20:][::-1]
        rec={'sample':sample,'cluster':c,'n_cells':int((a.obs.cluster==c).sum()),'major_type':label,'top_genes':';'.join(a.var_names[top])}
        rec.update({f'score_{k}':float(v) for k,v in zip(PANELS,scores)});labels.append(rec)
    ann=pd.DataFrame(labels);ann.to_csv(ROOT/'04_qc'/f'{sample}_cluster_annotation.tsv',sep='\t',index=False)
    mapping=ann.set_index('cluster').major_type.to_dict();a.obs['major_type']=a.obs.cluster.astype(str).map(mapping)
    a.obs['malignancy']='Not_epithelial';epimask=a.obs.major_type.isin(['Epi','Acinar','Endocrine'])
    a.obs.loc[a.obs.major_type=='Epi','malignancy']='Epi_uncertain'
    a.obs.loc[a.obs.major_type=='Acinar','malignancy']='Acinar_like'
    a.obs.loc[a.obs.major_type=='Endocrine','malignancy']='Endocrine_like'
    refs=np.flatnonzero(a.obs.major_type.isin(REF));rng=np.random.default_rng(20260922)
    a.obs['cnv_score']=np.nan;a.obs['cnv_nchr']=0;a.obs['cnv_positive']=False
    if epimask.sum()>=10 and len(refs)>=100:
        import infercnvpy as cnv
        pos=pd.read_csv(ROOT/'00_docs/hg38_refGene.txt.gz',sep='\t',header=None,usecols=[2,4,5,12],names=['chromosome','start','end','gene_symbol'])
        pos=pos[pos.chromosome.isin(['chr'+str(i) for i in range(1,23)])].sort_values(['gene_symbol','chromosome','start']).drop_duplicates('gene_symbol').set_index('gene_symbol')
        for key in ['chromosome','start','end']:a.var[key]=a.var.gene_symbol.map(pos[key]).to_numpy()
        pool=rng.choice(refs,size=min(2000,len(refs)),replace=False);train=pool[:len(pool)//2];hold=pool[len(pool)//2:]
        a.obs['cnv_reference']='Observation';a.obs.iloc[train,a.obs.columns.get_loc('cnv_reference')]='Reference'
        selected=np.unique(np.r_[np.flatnonzero(epimask),pool]);cn=a[selected].copy()
        cnv.tl.infercnv(cn,reference_key='cnv_reference',reference_cat='Reference',window_size=101,step=10,dynamic_threshold=1.5,n_jobs=2,chunksize=1000)
        C=cn.obsm['X_cnv'].toarray();score=np.sqrt(np.mean(C*C,axis=1))
        held=cn.obs_names.isin(a.obs_names[hold]);base_score=score[held]
        threshold=max(float(np.quantile(base_score,.99)),float(np.median(base_score)+5*np.median(np.abs(base_score-np.median(base_score)))))
        chroms=sorted(cn.uns['cnv']['chr_pos'].items(),key=lambda x:x[1]);chrvals=[];chrnames=[]
        for i,(chrom,start) in enumerate(chroms):
            end=chroms[i+1][1] if i+1<len(chroms) else C.shape[1]
            if end-start<3:continue
            chrvals.append(np.median(C[:,start:end],axis=1));chrnames.append(chrom)
        chrvals=np.column_stack(chrvals);chrthr=np.quantile(np.abs(chrvals[held]),.99,axis=0)
        nchr=(np.abs(chrvals)>np.maximum(chrthr,1e-6)).sum(1)
        positive=(score>threshold)&(nchr>=2)
        cn.obs['cnv_score']=score;cn.obs['cnv_nchr']=nchr;cn.obs['cnv_positive']=positive
        # Require a recurrent within-sample CNV pattern, not an isolated high score.
        candidate=positive & cn.obs.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy()
        cn.obs['cnv_coherent']=False
        for c in clusters:
            ids=np.flatnonzero((cn.obs.cluster.astype(str)==c).to_numpy()&candidate)
            if len(ids)<10:continue
            profile=np.median(C[ids],axis=0);denom=np.linalg.norm(C[ids],axis=1)*np.linalg.norm(profile)
            cosine=(C[ids]@profile)/np.maximum(denom,1e-12)
            cn.obs.iloc[ids,cn.obs.columns.get_loc('cnv_coherent')]=cosine>.3
        for col in ['cnv_score','cnv_nchr','cnv_positive']:a.obs.loc[cn.obs_names,col]=cn.obs[col].to_numpy()
        mal=cn.obs_names[candidate & cn.obs.cnv_coherent.to_numpy()]
        a.obs.loc[mal,'malignancy']='Malignant_CNV'
        low=cn.obs_names[(score<=np.quantile(base_score,.95))&(cn.obs.major_type=='Epi').to_numpy()]
        a.obs.loc[low,'malignancy']='Ductal_CNVlow'
        pd.DataFrame(chrvals,index=cn.obs_names,columns=chrnames).to_csv(ROOT/'04_qc'/f'{sample}_chromosome_cnv.tsv.gz',sep='\t')
        audit['cnv']={'method':'infercnvpy 0.6.1, 101-gene window, step 10','n_reference_train':len(train),'n_reference_holdout':len(hold),'score_threshold':threshold,'heldout_q99':float(np.quantile(base_score,.99)),'n_chr_required':2,'coherence_cosine_gt':.3,'coherence_min_cells':10,'n_malignant':len(mal)}
        cn.write_h5ad(ROOT/'03_processed'/f'GSE205013_{sample}_cnv.h5ad',compression='gzip');del cn,C
    else:audit['cnv']={'status':'insufficient epithelial cells or normal reference; uncertain labels retained'}
    pb={};norm={};numbers={};strict=a.obs.qc_strict.to_numpy();epi=a.obs.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy();mal=(a.obs.malignancy=='Malignant_CNV').to_numpy()
    masks={'Malignant':mal,'EpiAll':epi,'MalStrictQC':mal&strict,'EpiStrictQC':epi&strict}
    masks.update({'Type_'+t:(a.obs.major_type==t).to_numpy() for t in a.obs.major_type.unique()})
    for group,mask in masks.items():
        numbers[group]=int(mask.sum());pb[group]=np.asarray(a.layers['counts'][mask].sum(0,dtype=np.int64)).ravel()
        norm[group]=np.asarray(a.X[mask].mean(0)).ravel() if mask.any() else np.full(a.n_vars,np.nan)
    np.savez_compressed(ROOT/'03_processed'/f'GSE205013_{sample}_pseudobulk.npz',genes=a.var_names.to_numpy(dtype=str),**pb)
    np.savez_compressed(ROOT/'03_processed'/f'GSE205013_{sample}_meanlog.npz',genes=a.var_names.to_numpy(dtype=str),**norm)
    savejson(ROOT/'03_processed'/f'GSE205013_{sample}_cellnumbers.json',numbers)
    audit['n_singlets']=a.n_obs;audit['cell_types']=a.obs.major_type.value_counts().to_dict();audit['malignancy_counts']=a.obs.malignancy.value_counts().to_dict()
    for col in ['cluster','major_type','malignancy','cnv_score','cnv_nchr','cnv_positive']:all_obs[col]=a.obs[col].reindex(all_obs.index)
    all_obs.to_csv(ROOT/'04_qc'/f'{sample}_all_cell_annotations.tsv.gz',sep='\t',index_label='cell_id')
    a.write_h5ad(ROOT/'03_processed'/f'GSE205013_{sample}_annotated.h5ad',compression='gzip')
    savejson(ROOT/'04_qc'/f'{sample}_qc_audit.json',audit);savejson(final,{'sample':sample,'status':'SAMPLE_ANALYSIS_COMPLETE','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_singlets':a.n_obs,'n_malignant':int(mal.sum())})
    print('SAMPLE_ANALYSIS_COMPLETE',sample,a.n_obs,int(mal.sum()),flush=True)
if __name__=='__main__':
    for sample in sys.argv[1:]:run(sample)
    print('SAMPLE_BATCH_COMPLETE',flush=True)
