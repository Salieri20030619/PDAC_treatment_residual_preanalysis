from pathlib import Path
import json,re,numpy as np,pandas as pd,h5py,datetime
from anndata.io import read_elem
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
path=ROOT/'02_raw/GSE202051/GSE202051_totaldata-final-toshare.h5ad'
with h5py.File(path) as f:
    obs=read_elem(f['obs']);var=read_elem(f['var'])
    assert obs.index.is_unique and var.index.is_unique
    assert np.array_equal(f['X/indptr'][:],f['layers/counts/indptr'][:])
    rng=np.random.default_rng(20260922);chosen=[]
    for pid in obs.pid.unique():
        ids=np.flatnonzero((obs.pid==pid).to_numpy());chosen.extend(rng.choice(ids,size=min(5,len(ids)),replace=False).tolist())
    errors=[];countsums=[]
    for row in chosen:
        start,end=f['X/indptr'][row:row+2];ct=f['layers/counts/data'][start:end];x=f['X/data'][start:end]
        assert np.array_equal(f['X/indices'][start:end],f['layers/counts/indices'][start:end])
        assert (ct==np.rint(ct)).all() and (ct>=0).all()
        total=float(obs.n_counts.iloc[row]);target=np.log1p(ct*1e4/total)
        errors.append(float(np.max(np.abs(x-target))))
        countsums.append({'cell_id':str(obs.index[row]),'pid':str(obs.pid.iloc[row]),'retained_gene_count_sum':float(ct.sum()),'reported_original_UMI':total,'X_lnTP10K_max_abs_error':errors[-1]})
    assert max(errors)<1e-5,max(errors)
    obsm={k:list(f['obsm'][k].attrs.get('shape',f['obsm'][k].shape if isinstance(f['obsm'][k],h5py.Dataset) else [])) for k in f['obsm']}
pd.DataFrame(countsums).to_csv(ROOT/'04_qc/GSE202051_layer_numeric_verification.tsv',sep='\t',index=False)
meta=pd.read_csv(ROOT/'01_metadata/patient_sample_metadata.tsv',sep='\t',keep_default_na=False)
patients=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
registered=set(patients[(patients.dataset=='GSE202051')&(patients.comparison_scope=='Primary')].patient_id)
samplemap={};mapping=[]
for sample,group in obs.groupby('sampleid',observed=True):
    pids=group.pid.astype(str).unique();assert len(pids)==1
    pid=pids[0];assert re.fullmatch('[UT][0-9]+',pid)
    canonical=f'PDAC_{pid[0]}_{int(pid[1:])}';assert canonical in registered
    treatments=group.treatment_status.astype(str).unique();responses=group.response.astype(str).unique();assert len(treatments)==len(responses)==1
    md=patients[(patients.dataset=='GSE202051')&(patients.patient_id==canonical)].iloc[0]
    observed=treatments[0].upper();expected=md.regimen_group.upper()
    assert observed==expected or (expected=='OTHER' and observed in ['CRTX','GART','RT']),(canonical,observed,expected)
    samplemap[str(sample)]=canonical
    mapping.append({'source_sample_id':str(sample),'source_pid':pid,'patient_id':canonical,'n_cells':len(group),'treatment_in_object':treatments[0],'regimen_group_clinical':md.regimen_group,'response_in_object':responses[0]})
assert set(samplemap.values())==registered
pd.DataFrame(mapping).to_csv(ROOT/'01_metadata/GSE202051_sampleid_pid_crosswalk.tsv',sep='\t',index=False)
cellmap={'Malignant':'Epi','Ductal':'Epi','Ductal (atypical)':'Epi','ADM':'Epi','Acinar':'Acinar','CAF':'CAF','myCAF':'CAF','Pericyte':'Pericyte','Vascular smooth muscle':'VSMC','Vascular':'Endo','Lymphatic':'Endo','Macrophage':'Myeloid','Dendritic':'Myeloid','Neutrophil':'Myeloid','CD8+ T':'Tcell','CD4+ T':'Tcell','Treg':'Tcell','Natural killer':'NK','B':'Bcell','Plasma':'Plasma','Mast':'Mast','Schwann':'Schwann','Adipocyte':'Adipocyte','Intra-pancreatic neurons':'Neuron','Beta':'Endocrine','Alpha':'Endocrine','Gamma':'Endocrine','Delta':'Endocrine','Epsilon':'Endocrine','Hormone-negative neuroendocrine':'Endocrine'}
assert set(obs['Level 2 Annotation'].astype(str))==set(cellmap)
cfg={'audit_approved':True,'expression_layer':'X','expression_scale':'ln_TP10K','counts_layer':'counts','gene_symbol_column':None,'sample_column':'sampleid','patient_column_in_source':'pid','patient_mapping':samplemap,'celltype_column':'Level 2 Annotation','celltype_mapping':cellmap,'malignant_labels':['Malignant'],'uncertain_labels':['Ductal (atypical)','ADM'],'umi_column':'n_counts','ngenes_column':'n_genes','mt_column':'percent_mito','mt_multiplier':100,'source':'GSE202051 totaldata-final-toshare.h5ad; numeric verification of X=ln(1+counts*10000/reported UMI) across five cells per patient; count integer integrity checked fully during processing','numeric_max_abs_error':max(errors),'n_numeric_check_cells':len(errors),'published_CNV_evidence':{'obs':['celltype_infercnv','cnv_leiden','cnv_score','Level 1 Annotation','Level 2 Annotation'],'obsm':obsm,'note':'cnv_score has 22 cluster-level values, not 224988 independent CNV measurements'},'scope_note':'Object has 44 sampleid labels but 43 pid labels; all 43 are PDAC tumor patients. Normal GEO samples are absent. T13 has two sample aliases, which are merged.'}
(ROOT/'05_analysis/GSE202051_layer_selection.json').write_text(json.dumps(cfg,indent=2,default=lambda x:int(x)))
# Enrich response metadata from the published object, keeping RECIST and histological categories distinct.
for frame in [meta,patients]:
    frame['response_definition']=np.where(frame.dataset=='GSE205013','RECIST reported in clinical supplement','')
    for pid,group in obs.groupby('pid',observed=True):
        response=group.response.astype(str).unique();assert len(response)==1
        canonical=f'PDAC_{str(pid)[0]}_{int(str(pid)[1:])}'
        sel=(frame.dataset=='GSE202051')&(frame.patient_id==canonical)
        if str(pid).startswith('T'):
            frame.loc[sel,'response']=response[0];frame.loc[sel,'response_definition']='Published histological response category; not RECIST'
    frame['unknown_fields']=frame.apply(lambda r:';'.join(k for k in ['patient_id','regimen','sampling_time','response','procedure','stage','preservation'] if not r[k]),axis=1)
meta.to_csv(ROOT/'01_metadata/patient_sample_metadata.tsv',sep='\t',index=False);patients.to_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',index=False)
evidence=pd.read_csv(ROOT/'01_metadata/metadata_field_provenance.tsv',sep='\t',keep_default_na=False)
mask=(evidence.dataset=='GSE202051')&(evidence.field=='response')
lookup=meta.set_index('GSM')['response'];evidence.loc[mask,'value']=evidence.loc[mask,'GSM'].map(lookup)
evidence.loc[mask,'source_url']='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202051'
evidence.loc[mask,'basis']='Published totaldata-final-toshare.h5ad obs.response, unique within author pid; histologic response not RECIST; untreated left blank'
evidence.loc[mask,'unknown']=evidence.loc[mask,'value']==''
evidence.to_csv(ROOT/'01_metadata/metadata_field_provenance.tsv',sep='\t',index=False)
print('OBS_SAMPLE_COUNT',obs.sampleid.nunique(),'PATIENT_COUNT',obs.pid.nunique(),'MAX_LN_NORMALIZATION_ERROR',max(errors))
print(pd.DataFrame(mapping).to_string(index=False));print('VALIDATION_LAYER_SELECTION_COMPLETE')
