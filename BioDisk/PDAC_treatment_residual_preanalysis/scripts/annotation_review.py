from pathlib import Path
import pandas as pd,numpy as np,json,sys
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
tables=[]
for file in sorted((ROOT/'04_qc').glob('P*_cluster_annotation.tsv')):
    df=pd.read_csv(file,sep='\t');scores=df.filter(regex='^score_')
    order=np.argsort(scores.to_numpy(),axis=1)
    df['score_top']=np.take_along_axis(scores.to_numpy(),order[:,-1:],axis=1).ravel()
    df['score_second']=np.take_along_axis(scores.to_numpy(),order[:,-2:-1],axis=1).ravel()
    df['second_marker_panel']=np.array(scores.columns)[order[:,-2]]
    df['annotation_margin']=df.score_top-df.score_second
    df['review_flag']=np.where(df.major_type=='Unknown','Uncertain',np.where(df.annotation_margin<.25,'Close_marker_scores','Marker_supported'))
    tables.append(df)
combined=pd.concat(tables,ignore_index=True)
combined.to_csv(ROOT/'04_qc/all_clusters_annotation_review.tsv',sep='\t',index=False)
print(combined.groupby(['major_type','review_flag'],observed=True).agg(n_clusters=('cluster','size'),n_cells=('n_cells','sum')).to_string())
print('FLAGGED_CLUSTERS')
show=(combined.review_flag!='Marker_supported')
if len(sys.argv)>1:show&=combined['sample'].isin(sys.argv[1:])
print(combined.loc[show,['sample','cluster','n_cells','major_type','second_marker_panel','annotation_margin','top_genes']].to_string(index=False))
print('ANNOTATION_REVIEW_COMPLETE')
