from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
OUT=ROOT/'07_figures';OUT.mkdir(exist_ok=True)
sns.set_theme(style='ticks',context='paper',font='DejaVu Sans',font_scale=1.05)
plt.rcParams.update({'figure.dpi':140,'savefig.dpi':220,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
colors={'Untreated':'#4477AA','Treated':'#CC6677','Unknown':'#888888'}
meta=pd.read_csv(ROOT/'01_metadata/patient_metadata.tsv',sep='\t',keep_default_na=False)
q=pd.read_csv(ROOT/'06_tables/GSE205013_sample_QC.tsv',sep='\t').merge(meta[['dataset','patient_id','treatment_status','procedure','comparison_scope']],on=['dataset','patient_id'])
fig,axes=plt.subplots(2,2,figsize=(13,7.5),constrained_layout=True)
for ax,key,ylabel in zip(axes.ravel(),['median_genes','median_UMI','median_mt','malignant_cells'],['Genes','UMI','MT (%)','Mal. cells']):
    for tx in colors:
        d=q[q.treatment_status==tx];ax.scatter(d.patient_id,d[key],color=colors[tx],s=30,label=tx if tx!='Untreated' else 'NTx')
    ax.set_ylabel(ylabel);ax.tick_params(axis='x',rotation=90)
axes[0,0].legend(frameon=False,ncol=3,loc='lower left',bbox_to_anchor=(0,1.01),fontsize=8)
fig.savefig(OUT/'01_patient_QC.png',bbox_inches='tight');plt.close(fig)
ct=pd.read_csv(ROOT/'06_tables/GSE205013_major_cell_types.tsv',sep='\t')
wide=ct.pivot(index='patient_id',columns='cell_type',values='n_cells').fillna(0);wide=wide.div(wide.sum(1),axis=0)
fig,ax=plt.subplots(figsize=(13,4.5));wide.plot.bar(stacked=True,ax=ax,width=.85,colormap='tab20');ax.set_ylabel('Cell frac.');ax.set_xlabel('Patient');ax.legend(frameon=False,bbox_to_anchor=(1.01,1),loc='upper left',fontsize=8);fig.tight_layout();fig.savefig(OUT/'02_major_types.png',bbox_inches='tight');plt.close(fig)
scores=pd.read_csv(ROOT/'06_tables/patient_mechanism_scores.tsv',sep='\t')
cross=pd.read_csv(ROOT/'06_tables/cross_cohort_direction.tsv',sep='\t')
stats=pd.read_csv(ROOT/'06_tables/patient_mechanism_statistics.tsv',sep='\t')
sets=json.loads((ROOT/'00_docs/genesets_frozen.json').read_text())
abbr={'HALLMARK_UNFOLDED_PROTEIN_RESPONSE':'H_UPR','REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR':'R_UPR','REACTOME_PROTEIN_FOLDING':'Fold','REACTOME_PROTEASOME_ASSEMBLY':'Proteasome','REACTOME_PROTEIN_UBIQUITINATION':'Ubiquitin','REACTOME_AUTOPHAGY':'Autophagy','REACTOME_LYSOSOME_VESICLE_BIOGENESIS':'Lysosome','REACTOME_TRANSLATION':'Translat','HALLMARK_HYPOXIA':'Hypoxia','HALLMARK_APOPTOSIS':'Apoptosis','HALLMARK_E2F_TARGETS':'E2F','HALLMARK_G2M_CHECKPOINT':'G2M','MOFFITT_CLASSICAL':'Classical','MOFFITT_BASAL':'Basal','REACTOME_ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS':'ATF4','REACTOME_ATF6_ATF6_ALPHA_ACTIVATES_CHAPERONE_GENES':'ATF6','REACTOME_IRE1ALPHA_ACTIVATES_CHAPERONES':'IRE1','REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY':'GCN2_ISR','REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY':'CMA','REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND':'HSP90','REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS':'HeatStress'}
abbr.update({'HWANG_NEURAL_LIKE_PROGENITOR':'NRP','HWANG_NEUROENDOCRINE_LIKE':'NeuroE','HWANG_MESENCHYMAL':'Mesench','HALLMARK_IL6_JAK_STAT3_SIGNALING':'IL6_STAT3'})
for i,key in enumerate(sets,1):abbr.setdefault(key,'GS'+str(i).zfill(2))
pd.DataFrame([{'abbreviation':v,'geneset':k} for k,v in abbr.items()]).to_csv(OUT/'figure_abbreviations.tsv',sep='\t',index=False)
priority=['HALLMARK_UNFOLDED_PROTEIN_RESPONSE','REACTOME_PROTEIN_FOLDING','REACTOME_PROTEASOME_ASSEMBLY','REACTOME_AUTOPHAGY','REACTOME_TRANSLATION','MOFFITT_BASAL']
rank=cross.sort_values(['candidate_gate','same_direction','FDR_validation','FDR_discovery'],ascending=[False,False,True,True])
leaders=['REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY','HWANG_MESENCHYMAL','HWANG_NEURAL_LIKE_PROGENITOR']
shown=list(dict.fromkeys(leaders+priority))[:9]
fig,axes=plt.subplots(3,3,figsize=(12,9),constrained_layout=True)
rng=np.random.default_rng(20260922)
for ax,key in zip(axes.ravel(),shown):
    for di,ds in enumerate(['GSE205013','GSE202051']):
        d=scores[(scores.dataset==ds)&(scores.cell_group=='Malignant')&scores.analysis_eligible&(scores.comparison_scope=='Primary')].copy()
        if ds=='GSE202051':d=d[d.regimen_group.isin(['Untreated','CRT'])]
        for ti,tx in enumerate(['Untreated','Treated']):
            t=d[d.treatment_status==tx];x=di*3+ti
            ax.scatter(x+rng.uniform(-.12,.12,len(t)),t[key],s=24,color=colors[tx],alpha=.85)
            if len(t):ax.plot([x-.18,x+.18],[t[key].median()]*2,color='black',lw=1.4)
    ax.set_title(abbr[key]);ax.set_xticks([0,1,3,4],['D_NTx','D_Tx','V_NTx','V_CRT'],rotation=30);ax.set_ylabel('Rank score')
for ax in axes.ravel()[len(shown):]:ax.set_visible(False)
fig.savefig(OUT/'03_patient_states.png',bbox_inches='tight');plt.close(fig)
proteo=[k for k,v in sets.items() if v['category']=='Proteostasis_candidate']+['HALLMARK_UNFOLDED_PROTEIN_RESPONSE']
shown=list(dict.fromkeys(leaders+proteo+['HALLMARK_HYPOXIA','HALLMARK_APOPTOSIS','MOFFITT_CLASSICAL','MOFFITT_BASAL']))
d=cross.set_index('geneset').reindex(shown).dropna(subset=['effect_discovery'])
fig,axes=plt.subplots(1,2,figsize=(12,max(5,len(d)*.27)),gridspec_kw={'width_ratios':[1.3,1]},constrained_layout=True)
y=np.arange(len(d))
axes[0].axvline(0,color='.6',lw=.8)
for offset,ds,cmp,color,label in [(-.13,'GSE205013','Primary','#4477AA','Discovery'),(.13,'GSE202051','CRT','#CC6677','Validation')]:
    st=stats[(stats.dataset==ds)&(stats.cell_group=='Malignant')&(stats.comparison==cmp)].set_index('geneset').reindex(d.index)
    axes[0].errorbar(st.effect,y+offset,xerr=[st.effect-st.ci_low,st.ci_high-st.effect],fmt='o',ms=3,color=color,lw=.7,label=label)
axes[0].set_yticks(y,[abbr[x] for x in d.index]);axes[0].set_xlabel('Tx effect');axes[0].invert_yaxis();axes[0].legend(frameon=False,fontsize=8)
matrix=d[['loo_same_direction_discovery','loo_same_direction_validation','same_direction','Resect_same_direction','StrictQC_same_direction','ValidationStrictQC_same_direction']].astype(float)
im=axes[1].imshow(matrix.to_numpy(),vmin=0,vmax=1,cmap='Blues',aspect='auto');axes[1].set_yticks(y,[]);axes[1].set_xticks(np.arange(6),['D_LOO','V_LOO','Cross','Resect','D_QC','V_QC'],rotation=35);fig.colorbar(im,ax=axes[1],label='Consist.',shrink=.65)
fig.savefig(OUT/'04_effect_stability.png',bbox_inches='tight');plt.close(fig)
captions='''# Figure captions\n\n01_patient_QC.png: Each point is one GSE205013 patient. Genes/UMI/MT are input-cell medians. Mal. cells are retained CNV-supported malignant candidates; no cell-level significance test is shown. NTx: untreated. MT: mitochondrial.\n\n02_major_types.png: Each bar represents one patient; fractions refer to retained singlets. Epi: epithelial, CAF: cancer-associated fibroblast, Endo: endothelial, Hepato: hepatocyte, NK: natural killer.\n\n03_patient_states.png: Each point is one eligible patient (at least 30 cells of the specified compartment); bars are medians. D: GSE205013 discovery; V: GSE202051 validation; NTx: untreated; Tx: treated; CRT: standard chemotherapy/radiotherapy. No paired observations are implied.\n\n04_effect_stability.png: Patient-level treatment effects and HC3 95% confidence intervals. D_LOO/V_LOO: fraction of estimable leave-one-patient-out fits retaining the full-model direction. Resect: resection subset; NoRT: exclusion of P10; StrictQC: stricter cell quality threshold. Direction agreement is not a significance claim. Missing comparisons remain blank.\n\nAll signature abbreviations are defined in figure_abbreviations.tsv. Genesets and versions are in genesets_frozen.json. Hwang-derived program comparisons within the validation cohort are not independent signature discovery validation.\n'''
captions=captions.replace('NoRT: exclusion of P10; StrictQC: stricter cell quality threshold.','Cross: direction agreement between cohorts; D_QC/V_QC: stricter quality subset in discovery/validation. NoRT is retained in the statistics table; P10 was already excluded from the malignant main set, so it provides no additional independent stability evidence.')
captions+='\n05_patient_CNV.png: Mean inferred CNV per chromosome and patient among CNV-supported malignant candidates. P08 has no retained malignant candidate and is absent. Mean CNV is a transcriptome-derived score, not DNA copy number. Chr.: chromosome.\n\nScores are compared within cohorts; absolute score levels between whole-cell and nuclear assays are not treatment effects. CMA: chaperone-mediated autophagy; Mesench: mesenchymal program; NRP: neural-like progenitor program.\n'
(OUT/'figure_captions.md').write_text(captions)
print('FIGURES_COMPLETE')
