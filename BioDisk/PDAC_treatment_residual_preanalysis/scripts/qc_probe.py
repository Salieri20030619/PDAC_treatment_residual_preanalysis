from pathlib import Path
import gzip,json,numpy as np,pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
rows=[]
for p in sorted((ROOT/'02_raw/GSE205013/early_matrices').glob('*matrix.mtx.gz')):
    prefix=p.name.replace('matrix.mtx.gz',''); sample=prefix.split('_')[1]
    with threadpool_limits(limits=4):
        with gzip.open(p) as f: x=mmread(f).tocsr().T.tocsr()
    genes=pd.read_csv(p.parent/(prefix+'features.tsv.gz'),sep='\t',header=None)[1].astype(str).to_numpy()
    umi=np.asarray(x.sum(1)).ravel(); ng=x.getnnz(1); mt=np.asarray(x[:,np.char.startswith(genes.astype(str),'MT-')].sum(1)).ravel()/np.maximum(umi,1)*100
    hb=np.asarray(x[:,np.isin(genes,['HBA1','HBA2','HBB','HBM','ALAS2'])].sum(1)).ravel()/np.maximum(umi,1)*100
    r={'sample':sample,'n_cells':x.shape[0],'n_genes':x.shape[1],'integer_counts':bool(np.all(x.data==np.floor(x.data))),'ngenes_quantile':np.quantile(ng,[0,.01,.1,.5,.9,.99,1]).tolist(),'umi_quantile':np.quantile(umi,[0,.01,.1,.5,.9,.99,1]).tolist(),'mt_quantile':np.quantile(mt,[0,.01,.1,.5,.9,.99,1]).tolist(),'paper_qc_kept':int(((ng>500)&(umi>1500)&(mt<15)&(hb<1)).sum())}
    rows.append(r);print(json.dumps(r),flush=True)
(ROOT/'04_qc/initial_matrix_qc.json').write_text(json.dumps(rows,indent=2))
import scanpy as sc,inspect,importlib.util
print('scrublet_signature',inspect.signature(sc.pp.scrublet));print('skimage_available',importlib.util.find_spec('skimage') is not None)
print('QC_PROBE_COMPLETE')
