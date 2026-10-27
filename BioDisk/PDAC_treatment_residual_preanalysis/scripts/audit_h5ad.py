from pathlib import Path
import h5py,json,numpy as np,pandas as pd
from anndata.io import read_elem
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
path=ROOT/'02_raw/GSE202051/GSE202051_totaldata-final-toshare.h5ad'
def audit(node):
    info={'encoding':str(node.attrs.get('encoding-type',''))}
    if isinstance(node,h5py.Dataset):
        info.update(shape=list(node.shape),dtype=str(node.dtype))
        x=node[:min(100,node.shape[0]),:min(5000,node.shape[1])].ravel() if len(node.shape)==2 else node[:100000]
    else:
        info['keys']=list(node.keys());info['shape']=list(map(int,node.attrs.get('shape',[])))
        if 'data' not in node:return info
        x=node['data'][:min(1000000,node['data'].shape[0])];info['nnz']=int(node['data'].shape[0]);info['dtype']=str(node['data'].dtype)
    if np.issubdtype(x.dtype,np.number):
        info.update(min=float(np.min(x)),max=float(np.max(x)),nonzero_fraction=float(np.mean(x!=0)),integer_fraction=float(np.mean(np.abs(x-np.rint(x))<1e-6)),examples=x[x!=0][:20].tolist())
    return info
with h5py.File(path) as f:
    out={'root_keys':list(f.keys()),'layers':{},'X':audit(f['X'])}
    if 'raw' in f and 'X' in f['raw']:out['raw_X']=audit(f['raw/X'])
    if 'layers' in f:
        for key in f['layers']:out['layers'][key]=audit(f['layers'][key])
    obs=read_elem(f['obs']);var=read_elem(f['var'])
    obs.to_csv(ROOT/'01_metadata/GSE202051_obs_raw.tsv.gz',sep='\t',index_label='cell_id')
    var.to_csv(ROOT/'01_metadata/GSE202051_var_raw.tsv.gz',sep='\t',index_label='gene_id')
    out['obs_shape']=list(obs.shape);out['var_shape']=list(var.shape);out['var_columns']=var.columns.tolist();out['var_preview']=var.head(5).astype(str).to_dict()
    out['obs_columns']={}
    for col in obs:
        v=obs[col];d={'dtype':str(v.dtype),'nunique':int(v.nunique()),'missing':int(v.isna().sum())}
        if v.nunique()<120:d['value_counts']={str(k):int(n) for k,n in v.value_counts(dropna=False).items()}
        elif pd.api.types.is_numeric_dtype(v):d['quantiles']=v.quantile([0,.01,.5,.99,1]).to_dict()
        else:d['examples']=v.astype(str).head(8).tolist()
        out['obs_columns'][col]=d
    if 'raw' in f and 'var' in f['raw']:
        rv=read_elem(f['raw/var']);out['raw_var_shape']=list(rv.shape);out['raw_var_preview']=rv.head(5).astype(str).to_dict();rv.to_csv(ROOT/'01_metadata/GSE202051_raw_var.tsv.gz',sep='\t')
    if 'uns' in f:
        out['uns_keys']=list(f['uns'].keys())
        for key in ['log1p','cnv','scrublet']:
            if key in f['uns']:out['uns_'+key]=str(read_elem(f['uns'][key]))[:3000]
(ROOT/'04_qc/GSE202051_h5ad_audit.json').write_text(json.dumps(out,indent=2,default=str))
print(json.dumps(out,indent=2,default=str));print('H5AD_STRUCTURE_AUDIT_COMPLETE')
