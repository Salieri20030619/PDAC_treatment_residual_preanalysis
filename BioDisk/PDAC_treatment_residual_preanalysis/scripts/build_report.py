from pathlib import Path
import json,pandas as pd,numpy as np,datetime,html
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
interpret=json.loads((ROOT/'05_analysis/scientific_interpretation.json').read_text())
assert interpret['reviewed'] is True
q1=pd.read_csv(ROOT/'06_tables/GSE205013_sample_QC.tsv',sep='\t')
q2=pd.read_csv(ROOT/'06_tables/GSE202051_sample_QC.tsv',sep='\t')
inc=pd.read_csv(ROOT/'06_tables/analysis_patient_inclusion.tsv',sep='\t',keep_default_na=False)
stat=pd.read_csv(ROOT/'06_tables/patient_mechanism_statistics.tsv',sep='\t')
cross=pd.read_csv(ROOT/'06_tables/cross_cohort_direction.tsv',sep='\t')
de=pd.read_csv(ROOT/'06_tables/pseudobulk_DE_audit.tsv',sep='\t')
ma=json.loads((ROOT/'logs/METADATA_AUDIT_PASS.json').read_text())
da=json.loads((ROOT/'logs/DOWNLOAD_AUDIT_PASS.json').read_text())
layer=json.loads((ROOT/'05_analysis/GSE202051_layer_selection.json').read_text())
manifest=pd.read_csv(ROOT/'01_metadata/download_manifest.tsv',sep='\t').fillna('')
def table(df):
    d=df.copy()
    for c in d:
        if pd.api.types.is_float_dtype(d[c]):d[c]=d[c].map(lambda x:'' if pd.isna(x) else f'{x:.4g}')
    columns=list(d.columns);lines=['| '+' | '.join(columns)+' |','|'+'|'.join(['---']*len(columns))+'|']
    for row in d.itertuples(index=False,name=None):lines.append('| '+' | '.join(str(x).replace('|','/').replace('\n',' ') for x in row)+' |')
    return '\n'.join(lines)
counts=[]
for ds in ['GSE205013','GSE202051']:
    d=inc[(inc.dataset==ds)&(inc.cell_group=='Malignant')&(inc.comparison_scope=='Primary')]
    for tx,g in d.groupby('treatment_status' if ds=='GSE205013' else 'regimen_group'):
        counts.append({'队列':ds,'组别':tx,'实际对象患者':g.patient_id.nunique(),'恶性细胞≥30患者':int(g.analysis_eligible.astype(str).str.lower().eq('true').sum()),'恶性候选细胞':int(g.n_cells.sum())})
counttable=pd.DataFrame(counts)
selected=interpret['reported_genesets']
result=stat[(stat.cell_group=='Malignant')&(((stat.dataset=='GSE205013')&(stat.comparison=='Primary'))|((stat.dataset=='GSE202051')&(stat.comparison=='CRT')))&stat.geneset.isin(selected)]
result=result[['dataset','comparison','geneset','n_untreated','n_treated','effect','pvalue','FDR','loo_same_direction']]
excluded=inc[(inc.cell_group=='Malignant')&(inc.comparison_scope=='Primary')&~inc.analysis_eligible.astype(str).str.lower().eq('true')][['dataset','patient_id','treatment_status','n_cells']]
downloadrows=[]
for name in ['GSE205013_RAW.tar','GSE202051_totaldata-final-toshare.h5ad.gz']:
    r=json.loads((ROOT/'logs'/f'{name}.download.json').read_text());downloadrows.append({'文件':name,'字节数':r['bytes'],'SHA256':r['sha256'],'状态':r['status']})
downloadtext='# 下载与完整性报告\n\n'+table(pd.DataFrame(downloadrows))+'''\n\n所有公网下载经服务器代理 http://127.0.0.1:7897。GSE205013 采用整包分段续传，未再逐 GSM 重复下载。81 个 tar 成员均通过 gzip 校验，27 个矩阵已实际读入并核对维度、条码与非负整数计数。GSE202051 仅下载组织对象 totaldata-final-toshare；类器官 010nuc/010orgCRT 文件只登记，未默认下载。文件大小与官方 HTTP Content-Length 对照，完整文件计算 SHA256，gzip/tar 和 HDF5 读入分别检查。官方目录未列出校验和，本地 SHA256 用于固定本次文件，不能表述为官方摘要核验。

下载曾发生 HTTP 503、超时与 TLS EOF，均保留在 logs。完整分段跳过，未完整分段断点续传，重试次数有上限。原始文件保留，不删除源数据。所有原始 GEO Series/Sample 元数据、论文数据声明、临床补充资料、有版本基因集和作者仓库提交记录均保留下载清单。

完整记录：01_metadata/download_manifest.tsv；表达文件各自的 headers、curl/segment 日志与 .download.json 位于 logs。DOWNLOAD_AUDIT_PASS 仅在两个表达文件完整、tar/gzip 校验及矩阵可读性核实后产生。
'''
(ROOT/'08_reports/download_integrity_report.md').write_text(downloadtext)
parts=[f'# 胰腺癌治疗相关恶性细胞状态：公开数据预分析\n\n生成日期：{datetime.datetime.now().isoformat(timespec="seconds")}。\n',
f'## 结论\n\n{interpret["headline"]}\n\n{interpret["proteostasis_assessment"]}\n',
'## 可用数据与实际有效患者\n\n'+table(counttable)+f'\n\nGSE205013 实际读入 {int(q1.input_cells.sum()):,} 个 Cell Ranger 条码，QC 与双细胞评估后保留 {int(q1.singlet_cells.sum()):,} 个细胞。GSE202051 实际组织对象包含 {int(q2.input_cells.sum()):,} 个已发表过滤细胞/核，映射到 {ma["GSE202051_expression_tumor_patients"]} 个肿瘤患者 ID。技术重复合并到患者，不把细胞、文库或空间 ROI 当作生物学样本量。\n\n'+interpret['data_availability_note'],
'\n### 未进入恶性细胞主比较的原发灶患者\n\n'+table(excluded)+'\n\n排除仅依据预先固定的细胞质量与每患者至少 30 个 CNV 支持恶性细胞/发表恶性标签细胞的要求；未知或证据不足细胞仍保留在注释和上皮整体敏感性分析中。\n',
'## 治疗相关状态的患者层面证据\n\n'+interpret['state_results']+'\n\n'+table(result)+'\n\n表中 effect 为患者基因秩评分的治疗效应，方向为已治疗减未治疗。GSE205013 Primary 模型处理可辨识的取材方式；GSE202051 CRT 仅比较标准 CRT 与未治疗。各比较内部全部合格基因集做 BH-FDR；完整检验表包括所有阴性结果。任一组少于 10 人的比较定位为探索性，以方向、效应量和留一患者稳定性为主。LOO 列为可估计留一患者模型保留原方向的比例。\n',
'## 跨队列、取材和质量敏感性\n\n'+interpret['stability_results']+'\n\n'+table(cross[cross.geneset.isin(selected)][['geneset','same_direction','loo_same_direction_discovery','loo_same_direction_validation','Resect_same_direction','NoRT_same_direction','StrictQC_same_direction','EpiAll_same_direction']])+'\n\n留一患者完整效应保存在 leave_one_patient_out.tsv；取材、分期、严格 QC、上皮整体及背景状态调整分别报告。某患者原本未进入主模型时，剔除该患者不算额外独立稳定性证据。Hwang 状态程序源于验证队列，其在原队列内的方向重复不视为独立程序发现验证。\n',
'## 下一轮机制与多组学衔接\n\n'+interpret['mechanism_priorities']+'\n\n'+interpret['next_experiments']+'\n',
'## 方法\n\n数据来自 GEO 官方系列/样本记录和公开补充文件，参考原论文及临床附表。先逐条审计患者、样本、文库、取材方式和用药，随后分别处理两队列。GSE205013 按样本记录实际 QC 分布，使用 genes >500、UMI >1500、MT <15%、红细胞相关转录本 <1% 的阈值；每个 Cell Ranger 样本运行一次 Scanpy Scrublet。全细胞进行归一化、HVG、PCA、邻接图和 Leiden 聚类，以主要谱系标志注释。原始计数保留用于 pseudobulk。\n\n恶性候选结合上皮/腺泡/内分泌表达谱和 infercnvpy 0.6.1 推断 CNV：hg38 refGene 坐标，101 基因窗口、步长 10。正常免疫/间质参考经谱系标志门槛及低上皮标志筛选，排除 Mixed；分为训练和留出参考，各最多 1000 细胞。候选总体 CNV 分数需超过留出参考 99% 分位与中位数加 5 MAD 的较高者，至少两条染色体窗口 RMS 高于相应留出参考 99% 分位，同聚类至少 10 个候选且模式余弦相似度 >0.3。低 CNV 导管保留 Ductal_CNVlow，不自动判为正常或恶性。\n\nGSE202051 先检查 X、raw、layers，再按已核实的表达层定义处理。所选评分层为 '+str(layer['expression_layer'])+'，数值尺度为 '+str(layer['expression_scale'])+'；用于计数分析的层为 '+str(layer.get('counts_layer') or '无已核实原始计数层')+'。计数为作者 CellBender 处理后的保留特征 UMI；215 个细胞的数值核查确认 X=ln(1+counts×10000/原始 UMI)，最大误差 4.77e-7。发表对象已做 CellBender/Scrublet 处理，采用 Level 2 Annotation，保留作者结合表达与 CNV 的细胞标签，并独立核对主要标志，不机械重复去双细胞。\n\n严格 QC：GSE205013 使用 genes >750、UMI >2500、MT <10%；GSE202051 在每患者作者保留细胞内去除最低 10% 基因数及最高 10% MT 比例范围。每个患者各细胞群分别汇总计数、平均 log-normalized 表达。基因差异使用 edgeR TMM 与稳健 quasi-likelihood 模型；评分为固定共同基因空间内的患者平均表达百分位秩。患者评分线性模型使用 HC3 稳健标准误和 t 分布区间；以治疗、取材方式及可辨识分期设计审计确定模型。BH-FDR 按队列、细胞群和比较的预设同类检验家族计算；基因与通路分开。MSigDB 2025.1.Hs 的 Hallmark/Reactome 与论文公布的 Moffitt/Hwang 基因集在治疗差异检验前冻结，覆盖率和完整成员另表保留。\n',
'## 讨论与研究边界\n\n'+interpret['limitations']+'\n\n两队列是患者间观察性比较，“已治疗”表示取样时治疗暴露，不能直接命名为耐药或持留细胞。P11 的 GEO/论文临床表矛盾保留为 Unknown，且属于肝转移，不加入原发灶主比较。P10 包含 SBRT，验证队列 CRT/CRTL/CRTN/CRTLN 等干预分别处理。两种转录组技术不构成已完成多组学；机制干预依赖仍需遗传扰动、救援与独立实验确认。\n\n'+interpret['novelty_note'],
'\n## 交付与复核入口\n\n- 下载及校验：download_integrity_report.md、01_metadata/download_manifest.tsv。\n- 患者—样本与证据：01_metadata/patient_sample_metadata.tsv、metadata_field_provenance.tsv、metadata_conflicts.tsv、expression_patient_audit.tsv。\n- 分组混杂：06_tables/patient_*_crosstab.tsv、design_matrix_audit.tsv、registered_samples_vs_patients.tsv。\n- 逐细胞结果：04_qc/*all_cell_annotations.tsv.gz；03_processed 中保留注释对象、CNV 矩阵和原始计数汇总。\n- 统计：06_tables/patient_mechanism_scores.tsv、patient_mechanism_statistics.tsv、*_edgeR.tsv.gz、cross_cohort_direction.tsv、leave_one_patient_out.tsv、background_adjusted_statistics.tsv。\n- 图件：07_figures；英文缩写见 figure_abbreviations.tsv 与 figure_captions.md。\n- 复现：scripts 与 logs/execution_ledger.jsonl 记录语法检查、运行命令、退出状态和完成标志。\n\n'+interpret['pending_items'],
'\n## 主要来源\n\n[GSE205013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205013)；[Werba 等原始研究](https://www.nature.com/articles/s41467-023-36296-4)；[GSE202051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202051)；[Hwang 等原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC10290535/)；[MSigDB 2025.1.Hs](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/)。\n']
report='\n'.join(parts)
(ROOT/'08_reports/PDAC_preanalysis_zh.md').write_text(report)
try:
    import markdown
    body=markdown.markdown(report,extensions=['tables','fenced_code'])
except ImportError:body='<pre>'+html.escape(report)+'</pre>'
style='body{max-width:1180px;margin:40px auto;padding:0 28px;font:16px/1.75 Arial,"Microsoft YaHei",sans-serif;color:#18212b}table{border-collapse:collapse;font-size:12px;display:block;overflow:auto}th,td{border:1px solid #ccd3db;padding:6px 9px}th{background:#edf1f5}h1,h2,h3{line-height:1.4}a{color:#185e97}pre{white-space:pre-wrap}'
(ROOT/'08_reports/PDAC_preanalysis_zh.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>PDAC公开数据预分析</title><style>'+style+'</style><body>'+body+'</body></html>')
print('CHINESE_REPORT_COMPLETE')
