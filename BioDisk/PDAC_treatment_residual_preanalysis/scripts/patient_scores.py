from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import rankdata,mannwhitneyu,spearmanr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
SETS=json.loads((ROOT/'00_docs/genesets_frozen.json').read_text())
META=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
def read_cohort(ds):
    rows=[];vectors=[];counts=[];genes=None
    for p in sorted((ROOT/'03_processed').glob(ds+'_*_meanlog.npz')):
        patient=p.name[len(ds)+1:-len('_meanlog.npz')]
        match=META[(META.dataset==ds)&(META.patient_id==patient)]
        if len(match)!=1:raise ValueError('Patient mapping not unique: '+patient)
        md=match.iloc[0].to_dict();data=np.load(p);g=data['genes'].astype(str)
        if genes is None:genes=g
        assert np.array_equal(g,genes)
        nums=json.loads((ROOT/'03_processed'/f'{ds}_{patient}_cellnumbers.json').read_text())
        pbpath=ROOT/'03_processed'/f'{ds}_{patient}_pseudobulk.npz';pb=np.load(pbpath) if pbpath.exists() else None
        for group in ['Malignant','EpiAll','MalStrictQC','EpiStrictQC']:
            if group not in data:continue
            row=dict(md,cell_group=group,n_cells=nums[group],analysis_eligible=nums[group]>=30)
            rows.append(row);vectors.append(data[group]);counts.append(pb[group] if pb is not None and group in pb else None)
    return pd.DataFrame(rows),np.array(vectors),genes,counts
def design(meta,mode='procedure'):
    d=pd.DataFrame({'Intercept':1.,'Treated':(meta.treatment_status=='Treated').astype(float).to_numpy()},index=meta.index)
    if mode=='procedure' and meta.procedure.nunique()>1:
        z=pd.get_dummies(meta.procedure,prefix='Proc',drop_first=True,dtype=float);d=pd.concat([d,z],axis=1)
    if mode=='stage' and meta.stage.nunique()>1:
        z=pd.get_dummies(meta.stage,prefix='Stage',drop_first=True,dtype=float);d=pd.concat([d,z],axis=1)
        if meta.procedure.nunique()>1:
            candidate=pd.concat([d,pd.get_dummies(meta.procedure,prefix='Proc',drop_first=True,dtype=float)],axis=1)
            if np.linalg.matrix_rank(candidate.to_numpy())==candidate.shape[1]:d=candidate
    return d
def fit(y,m,mode,inference=True):
    x=design(m,mode)
    if np.linalg.matrix_rank(x.to_numpy())<x.shape[1] or len(m)-x.shape[1]<(3 if inference else 1):return None,x
    if inference and np.max(np.sum(x.to_numpy()*np.linalg.pinv(x.to_numpy()).T,axis=1))>1-1e-8:return None,x
    mod=sm.OLS(y,x).fit(cov_type='HC3',use_t=True) if inference else sm.OLS(y,x).fit()
    return mod,x
def tests(m,scores,ds,group,name,mode):
    out=[];loo=[]
    m=m.reset_index(drop=True);scores=scores.reset_index(drop=True)
    n0=int((m.treatment_status=='Untreated').sum());n1=int((m.treatment_status=='Treated').sum())
    if min(n0,n1)<2:return out,loo
    for key in scores:
        y=scores[key].to_numpy();f,x=fit(y,m,mode)
        rec=dict(dataset=ds,cell_group=group,comparison=name,geneset=key,n_untreated=n0,n_treated=n1,analysis_class='Exploratory' if min(n0,n1)<10 else 'Patient_level_observational',design=' + '.join(x.columns),design_rank=int(np.linalg.matrix_rank(x)),n_parameters=x.shape[1],mean_untreated=float(y[m.treatment_status=='Untreated'].mean()),mean_treated=float(y[m.treatment_status=='Treated'].mean()))
        rec['raw_delta']=rec['mean_treated']-rec['mean_untreated']
        if f is None:rec.update(status='Not_estimable',effect=np.nan,pvalue=np.nan);out.append(rec);continue
        effect=float(f.params['Treated']);pv=float(f.pvalues['Treated']);ci=f.conf_int().loc['Treated']
        rec.update(status='OK',effect=effect,pvalue=pv,ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]))
        u=mannwhitneyu(y[m.treatment_status=='Treated'],y[m.treatment_status=='Untreated'],alternative='two-sided',method='auto')
        rec['rank_biserial']=float(2*u.statistic/(n0*n1)-1);rec['MW_pvalue_sensitivity']=float(u.pvalue)
        effects=[]
        for i,patient in enumerate(m.patient_id):
            keep=np.arange(len(m))!=i;lf,lx=fit(y[keep],m.loc[keep].reset_index(drop=True),mode,inference=False)
            val=float(lf.params['Treated']) if lf is not None else np.nan
            loo.append(dict(dataset=ds,cell_group=group,comparison=name,geneset=key,omitted_patient=patient,effect=val,estimable=lf is not None))
            if np.isfinite(val):effects.append(val)
        rec.update(loo_n_estimable=len(effects),loo_estimable_fraction=len(effects)/len(m),loo_same_direction=float(np.mean(np.sign(effects)==np.sign(effect))) if effects else np.nan,loo_effect_min=min(effects) if effects else np.nan,loo_effect_max=max(effects) if effects else np.nan)
        out.append(rec)
    valid=[i for i,r in enumerate(out) if np.isfinite(r['pvalue'])]
    if valid:
        fdr=multipletests([out[i]['pvalue'] for i in valid],method='fdr_bh')[1]
        for i,q in zip(valid,fdr):out[i]['FDR']=q
    return out,loo
cohorts={ds:read_cohort(ds) for ds in ['GSE205013','GSE202051']}
cohorts={ds:v for ds,v in cohorts.items() if not v[0].empty}
if not cohorts:raise ValueError('No completed cohorts')
# A fixed shared gene universe makes pathway scores comparable across available cohorts.
universes=[]
for ds,(m,a,g,c) in cohorts.items():
    eligible=m.analysis_eligible.to_numpy()&(m.cell_group=='Malignant').to_numpy()
    detectable=np.nanmean(a[eligible],axis=0)>0 if eligible.any() else np.nanmean(a,axis=0)>0
    universes.append(set(g[detectable]))
common=sorted(set.intersection(*universes));assert len(common)>5000
(ROOT/'05_analysis/score_gene_universe.txt').write_text('\n'.join(common))
univ=set(common);coverage=[];usable={}
for key,gs in SETS.items():
    present=sorted(set(gs['genes'])&univ);ratio=len(present)/len(gs['genes'])
    use=len(present)>=5 and ratio>=.3
    coverage.append(dict(geneset=key,n_source=len(gs['genes']),n_detected=len(present),coverage=ratio,used=use,source=gs['source'],version=gs['version'],category=gs['category'],detected_genes=';'.join(present)))
    if use:usable[key]=present
pd.DataFrame(coverage).to_csv(ROOT/'06_tables/geneset_coverage.tsv',sep='\t',index=False)
all_score=[];all_tests=[];all_loo=[];all_inclusion=[];all_design=[]
for ds,(meta,expr,genes,counts) in cohorts.items():
    ix={g:i for i,g in enumerate(genes)};expr=expr[:,[ix[g] for g in common]]
    ranks=np.array([rankdata(x,method='average')/len(common) if np.isfinite(x).all() else np.full(len(common),np.nan) for x in expr])
    ci={g:i for i,g in enumerate(common)}
    scores=pd.DataFrame({key:np.mean(ranks[:,[ci[g] for g in gs]],axis=1) for key,gs in usable.items()})
    combined=pd.concat([meta[['dataset','patient_id','cell_group','n_cells','treatment_status','procedure','stage','regimen_group','comparison_scope','analysis_eligible']].reset_index(drop=True),scores],axis=1)
    all_score.append(combined);all_inclusion.append(meta)
    # Counts are never reconstructed from normalized expression.
    for group in ['Malignant','EpiAll','MalStrictQC','EpiStrictQC']:
        keep=(meta.cell_group==group)&meta.analysis_eligible&(meta.comparison_scope=='Primary')&meta.treatment_status.isin(['Untreated','Treated'])
        groupmeta=meta.loc[keep].reset_index(drop=True);groupscores=scores.loc[keep].reset_index(drop=True)
        if ds=='GSE205013':
            comparisons=[('Primary','procedure',np.ones(len(groupmeta),bool)),('Resection','none',(groupmeta.procedure=='Resection').to_numpy()),('NoRT','procedure',(groupmeta.patient_id!='P10').to_numpy()),('StageAdj','stage',np.ones(len(groupmeta),bool))]
        else:
            comparisons=[('CRT','none',groupmeta.regimen_group.isin(['Untreated','CRT']).to_numpy()),('CRTL','none',groupmeta.regimen_group.isin(['Untreated','CRTL']).to_numpy()),('AllTx','none',np.ones(len(groupmeta),bool))]
        for name,mode,sel in comparisons:
            m=groupmeta.loc[sel].reset_index(drop=True);s=groupscores.loc[sel].reset_index(drop=True)
            t,l=tests(m,s,ds,group,name,mode);all_tests+=t;all_loo+=l
            d=design(m,mode);requested=design(m,'procedure')
            if mode=='stage' and m.stage.nunique()>1:requested=pd.concat([requested,pd.get_dummies(m.stage,prefix='Stage',drop_first=True,dtype=float)],axis=1)
            elif mode!='stage':requested=d
            all_design.append(dict(dataset=ds,cell_group=group,comparison=name,n=len(m),columns=';'.join(d.columns),rank=int(np.linalg.matrix_rank(d)) if len(m) else 0,parameters=d.shape[1],residual_df=len(m)-d.shape[1],requested_columns=';'.join(requested.columns),requested_rank=int(np.linalg.matrix_rank(requested)) if len(m) else 0,requested_parameters=requested.shape[1],design_note='Stage replaces collinear procedure' if requested.shape[1]>d.shape[1] else '',patient_ids=';'.join(m.patient_id)))
        idx=np.flatnonzero(keep)
        if len(idx) and all(counts[i] is not None for i in idx):
            frame=pd.DataFrame(np.column_stack([counts[i] for i in idx]),index=genes,columns=meta.iloc[idx].patient_id)
            frame.to_csv(ROOT/'05_analysis'/f'{ds}_{group}_raw_pseudobulk.tsv.gz',sep='\t',index_label='gene')
            meta.iloc[idx].to_csv(ROOT/'05_analysis'/f'{ds}_{group}_pseudobulk_metadata.tsv',sep='\t',index=False)
score_table=pd.concat(all_score)
score_table.to_csv(ROOT/'06_tables/patient_mechanism_scores.tsv',sep='\t',index=False)
score_table[score_table.comparison_scope=='Liver_descriptive'].to_csv(ROOT/'06_tables/liver_metastasis_descriptive_scores.tsv',sep='\t',index=False)
pd.concat(all_inclusion).to_csv(ROOT/'06_tables/analysis_patient_inclusion.tsv',sep='\t',index=False)
pd.DataFrame(all_tests).to_csv(ROOT/'06_tables/patient_mechanism_statistics.tsv',sep='\t',index=False)
pd.DataFrame(all_loo).to_csv(ROOT/'06_tables/leave_one_patient_out.tsv',sep='\t',index=False)
pd.DataFrame(all_design).to_csv(ROOT/'06_tables/design_matrix_audit.tsv',sep='\t',index=False)
print('PATIENT_SCORE_ANALYSIS_COMPLETE',len(all_tests),flush=True)
