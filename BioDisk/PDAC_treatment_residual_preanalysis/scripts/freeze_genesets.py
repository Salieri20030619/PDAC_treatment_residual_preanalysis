from pathlib import Path
import json,hashlib,datetime,re
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
def readgmt(p):
    out={}
    for line in p.read_text().splitlines():
        name,source,*genes=line.split('\t');out[name]={'source':source,'genes':sorted(set(genes)),'version':'MSigDB 2025.1.Hs','category':'Hallmark_background' if name.startswith('HALLMARK') else 'Proteostasis_candidate'}
    return out
hall=readgmt(ROOT/'00_docs/hallmark_2025.1.gmt');react=readgmt(ROOT/'00_docs/reactome_2025.1.gmt')
selected=['UNFOLDED_PROTEIN_RESPONSE_UPR','PROTEIN_FOLDING','PROTEASOME_ASSEMBLY','PROTEIN_UBIQUITINATION','AUTOPHAGY','LYSOSOME_VESICLE_BIOGENESIS','TRANSLATION','RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY','ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS','ATF6_ATF6_ALPHA_ACTIVATES_CHAPERONE_GENES','IRE1ALPHA_ACTIVATES_CHAPERONES','CHAPERONE_MEDIATED_AUTOPHAGY','HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND','CELLULAR_RESPONSE_TO_HEAT_STRESS']
out=dict(hall)
for key in selected:
    name='REACTOME_'+key
    assert name in react,name
    out[name]=react[name]
    if key=='CELLULAR_RESPONSE_TO_HEAT_STRESS':out[name]['category']='Stress_background'
werba=json.loads((ROOT/'01_metadata/Werba2023_genesets.json').read_text())['Figure 2']
for j,name in enumerate(['MOFFITT_CLASSICAL','MOFFITT_BASAL']):
    genes=[row[j] for row in werba[2:] if row[j] and re.fullmatch('[A-Za-z0-9.-]+',str(row[j]))]
    out[name]={'genes':sorted(set(genes)),'version':'Werba2023 Supplementary Data 3 / Moffitt et al.','source':'https://www.nature.com/articles/s41467-023-36296-4','category':'Subtype_background'}
hwang=json.loads((ROOT/'01_metadata/Hwang2022_supplement_sheets.json').read_text())['Supplementary Table 2']
# Preserve all seven malignant lineage and seven state programs, not selected by current significance.
for j in range(1,15):
    name='HWANG_'+re.sub('[^A-Z0-9]+','_',hwang[3][j].upper()).strip('_')
    genes=[row[j] for row in hwang[4:204] if row[j]]
    out[name]={'genes':sorted(set(genes)),'version':'Hwang2022 Supplementary Table 2 top 200 weighted genes','source':'https://www.nature.com/articles/s41588-022-01134-8','category':'Published_PDAC_program','validation_caveat':'Programs derived in GSE202051; replication there is not independent signature discovery validation'}
dest=ROOT/'00_docs/genesets_frozen.json'
dest.write_text(json.dumps(out,indent=2))
(ROOT/'00_docs/genesets_freeze_record.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'n_sets':len(out),'selection':'All 50 Hallmarks, 14 mechanistic/stress Reactome sets, both Moffitt subtypes, all 14 Hwang malignant programs; fixed before patient treatment tests','min_detected_genes':5,'min_coverage':.3},indent=2))
for name,gs in out.items():print(name,len(gs['genes']))
print('GENESETS_FROZEN')
