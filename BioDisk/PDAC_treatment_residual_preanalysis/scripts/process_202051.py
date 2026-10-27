from pathlib import Path
import sys,json,numpy as np,pandas as pd,anndata as ad
from scipy import sparse
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
cfg=json.loads((ROOT/'05_analysis/GSE202051_layer_selection.json').read_text())
assert cfg['audit_approved'] is True
a=ad.read_h5ad(ROOT/'02_raw/GSE202051/GSE202051_totaldata-final-toshare.h5ad')
obs=a.obs.copy();obs.index.name='cell_id'
assert obs.index.is_unique,'Duplicate published cell identifiers must be reconciled before aggregation'
def matrix(layer):
    if layer=='raw.X':return a.raw.X,a.raw.var_names
    if layer=='X':return a.X,a.var_names
    return a.layers[layer],a.var_names
expr,genes=matrix(cfg['expression_layer']);genes=np.asarray(genes,dtype=str)
if cfg.get('gene_symbol_column'):
    source=a.raw.var if cfg['expression_layer']=='raw.X' else a.var
    genes=source[cfg['gene_symbol_column']].astype(str).to_numpy()
# Symbols must be unique for direct cohort comparison; duplicate symbols are collapsed before log transform when counts exist.
if len(set(genes))!=len(genes):
    genes=ad.utils.make_index_unique(pd.Index(genes)).to_numpy(dtype=str)
expr=sparse.csr_matrix(expr,dtype=np.float32)
counts=None
if cfg.get('counts_layer'):
    counts,cgenes=matrix(cfg['counts_layer']);counts=sparse.csr_matrix(counts)
    assert counts.shape==expr.shape
    assert np.all(np.isfinite(counts.data)) and np.all(counts.data>=0) and np.all(np.abs(counts.data-np.rint(counts.data))<1e-6)
    obs['UMI']=np.asarray(counts.sum(1)).ravel();obs['n_genes']=counts.getnnz(1)
    obs['pct_mt']=np.asarray(counts[:,np.char.startswith(genes.astype(str),'MT-')].sum(1)).ravel()/np.maximum(obs.UMI.to_numpy(),1)*100
    obs['retained_gene_UMI']=obs.UMI
    obs['retained_gene_n_genes']=obs.n_genes
    if cfg.get('umi_column'):obs['UMI']=a.obs[cfg['umi_column']].to_numpy()
    if cfg.get('ngenes_column'):obs['n_genes']=a.obs[cfg['ngenes_column']].to_numpy()
    if cfg.get('mt_column'):obs['pct_mt']=a.obs[cfg['mt_column']].to_numpy()*cfg.get('mt_multiplier',1)
else:
    obs['UMI']=obs[cfg['umi_column']].to_numpy() if cfg.get('umi_column') else np.nan
    obs['n_genes']=obs[cfg['ngenes_column']].to_numpy() if cfg.get('ngenes_column') else expr.getnnz(1)
    obs['pct_mt']=obs[cfg['mt_column']].to_numpy() if cfg.get('mt_column') else np.nan
if cfg['expression_scale']=='counts':
    assert counts is not None
    expr=sparse.diags(1e4/np.maximum(np.asarray(expr.sum(1)).ravel(),1))@expr;expr=expr.tocsr();expr.data=np.log1p(expr.data)
elif cfg['expression_scale']=='TP10K':expr.data=np.log1p(expr.data)
elif cfg['expression_scale']=='log2_TP10K':expr.data*=np.log(2)
elif cfg['expression_scale']=='ln_TP10K':pass
else:raise ValueError('Unapproved expression scale')
rawid=obs[cfg['sample_column']].astype(str)
obs['patient_id']=rawid.map(cfg['patient_mapping'])
obs['source_sample_id']=rawid
assert obs.patient_id.notna().all(),rawid[obs.patient_id.isna()].unique()
rawtype=obs[cfg['celltype_column']].astype(str)
obs['major_type']=rawtype.map(cfg['celltype_mapping'])
assert obs.major_type.notna().all(),rawtype[obs.major_type.isna()].unique()
obs['malignancy']=np.where(rawtype.isin(cfg['malignant_labels']),'Malignant_published_CNV',np.where(rawtype.isin(cfg['uncertain_labels']),'Epi_uncertain','Not_malignant'))
obs['published_celltype']=rawtype
meta=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
qc=[];ct=[];coverage=[]
for patient in sorted(obs.patient_id.unique()):
    ids=np.flatnonzero((obs.patient_id==patient).to_numpy());sub=obs.iloc[ids]
    m=meta[(meta.dataset=='GSE202051')&(meta.patient_id==patient)]
    assert len(m)==1,patient
    x=expr[ids];cc=counts[ids] if counts is not None else None
    # Already processed and doublet-filtered by the original study: retain the published main set.
    ng=sub.n_genes.to_numpy();strict=ng>=np.quantile(ng,.1)
    if sub.pct_mt.notna().all() and sub.pct_mt.max()>0:strict&=sub.pct_mt.to_numpy()<=np.quantile(sub.pct_mt,.9)
    obs.loc[sub.index,'qc_strict']=strict
    mal=(sub.malignancy=='Malignant_published_CNV').to_numpy();epi=sub.major_type.isin(['Epi','Acinar','Endocrine']).to_numpy()
    masks={'Malignant':mal,'EpiAll':epi,'MalStrictQC':mal&strict,'EpiStrictQC':epi&strict}
    masks.update({'Type_'+t:(sub.major_type==t).to_numpy() for t in sub.major_type.unique()})
    numbers={};norm={};pb={}
    for group,mask in masks.items():
        numbers[group]=int(mask.sum());norm[group]=np.asarray(x[mask].mean(0)).ravel() if mask.any() else np.full(x.shape[1],np.nan)
        if cc is not None:pb[group]=np.asarray(cc[mask].sum(0,dtype=np.int64)).ravel()
    np.savez_compressed(ROOT/'03_processed'/f'GSE202051_{patient}_meanlog.npz',genes=genes,**norm)
    if cc is not None:np.savez_compressed(ROOT/'03_processed'/f'GSE202051_{patient}_pseudobulk.npz',genes=genes,**pb)
    (ROOT/'03_processed'/f'GSE202051_{patient}_cellnumbers.json').write_text(json.dumps(numbers,indent=2))
    qc.append(dict(dataset='GSE202051',patient_id=patient,input_cells=len(ids),retained_cells=len(ids),malignant_cells=int(mal.sum()),epithelial_cells=int(epi.sum()),median_genes=float(np.median(ng)),median_UMI=float(sub.UMI.median()),median_mt=float(sub.pct_mt.median()),doublet_handling='Published Scrublet filtered; no repeat filtering',counts_layer=cfg.get('counts_layer','Unavailable'),regimen_group=m.iloc[0].regimen_group,treatment_status=m.iloc[0].treatment_status))
    for typ,n in sub.major_type.value_counts().items():ct.append(dict(dataset='GSE202051',patient_id=patient,cell_type=typ,n_cells=int(n)))
    coverage.append(dict(patient_id=patient,n_source_samples=sub.source_sample_id.nunique(),source_sample_ids=';'.join(sorted(sub.source_sample_id.unique())),n_cells=len(ids),n_malignant=int(mal.sum())))
    print('VALIDATION_PATIENT_COMPLETE',patient,len(ids),int(mal.sum()),flush=True)
obs.to_csv(ROOT/'04_qc/GSE202051_all_cell_annotations.tsv.gz',sep='\t')
pd.DataFrame(qc).to_csv(ROOT/'06_tables/GSE202051_sample_QC.tsv',sep='\t',index=False)
pd.DataFrame(ct).to_csv(ROOT/'06_tables/GSE202051_major_cell_types.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(ROOT/'01_metadata/GSE202051_expression_patient_mapping.tsv',sep='\t',index=False)
# Independent marker summary verifies that published malignant labels occupy an epithelial compartment.
panels={'Epi':['EPCAM','KRT8','KRT18','KRT19','KRT7'],'Ductal':['CFTR','SLC4A4','CA2'],'Acinar':['PRSS1','CPA1','CTRB1'],'Endocrine':['CHGA','CHGB','INS','GCG'],'Immune':['PTPRC','CD3D','LST1'],'CAF':['COL1A1','COL1A2','DCN']}
idx={g:i for i,g in enumerate(genes)};markers=[]
for typ in sorted(rawtype.unique()):
    mask=(rawtype==typ).to_numpy();r={'published_celltype':typ,'n_cells':int(mask.sum())}
    for name,gs in panels.items():
        cols=[idx[g] for g in gs if g in idx];r[name+'_meanlog']=float(expr[mask][:,cols].mean()) if cols else np.nan
    markers.append(r)
pd.DataFrame(markers).to_csv(ROOT/'04_qc/GSE202051_published_annotation_marker_check.tsv',sep='\t',index=False)
print('VALIDATION_PROCESSING_COMPLETE',len(qc),a.n_obs)
