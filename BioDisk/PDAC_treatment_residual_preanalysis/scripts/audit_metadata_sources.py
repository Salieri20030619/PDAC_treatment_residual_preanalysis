from pathlib import Path
import pandas as pd,json,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
meta=pd.read_csv(ROOT/'01_metadata/patient_sample_metadata.tsv',sep='\t',keep_default_na=False)
rows=[]
for _,r in meta.iterrows():
    for field in ['patient_id','sample_id','GSM','tissue_site','treatment_status','regimen','sampling_time','response','procedure','stage','preservation','assay']:
        source=r.source_url
        basis='GEO sample record'
        if field in ['regimen','response'] and r.dataset=='GSE205013' and r.clinical_source:
            source='https://www.nature.com/articles/s41467-023-36296-4';basis=r.clinical_source
        if field in ['regimen','stage'] and r.dataset=='GSE202051' and r.clinical_source:
            source='https://www.nature.com/articles/s41588-022-01134-8';basis=r.clinical_source+'; GEO regimen retained for Other'
        if field=='preservation' and r.dataset=='GSE205013':
            source='https://www.nature.com/articles/s41467-023-36296-4';basis='Methods fresh tissue processing; GEO extraction protocol'
        rows.append(dict(dataset=r.dataset,GSM=r.GSM,patient_id=r.patient_id,field=field,value=r[field],source_url=source,basis=basis,unknown=r[field]=='',conflict=r.metadata_conflict if field in ['treatment_status','regimen','response','sampling_time'] else ''))
pd.DataFrame(rows).to_csv(ROOT/'01_metadata/metadata_field_provenance.tsv',sep='\t',index=False)
clinical=pd.read_csv(ROOT/'01_metadata/Hwang_clinical_table1_extracted.tsv',sep='\t',keep_default_na=False)
clinical[clinical.ID.str.contains('\u2020')].to_csv(ROOT/'01_metadata/Hwang_DSP_only_excluded.tsv',sep='\t',index=False)
patients=meta[meta.patient_id!=''].drop_duplicates(['dataset','patient_id'])
counts=meta.groupby(['dataset','specimen_type'],dropna=False).agg(n_GSM=('GSM','size'),n_specimen_ID=('sample_id','nunique'),n_patient_ID=('patient_id',lambda x:len(set(x)-{''}))).reset_index()
counts.to_csv(ROOT/'06_tables/registered_samples_vs_patients.tsv',sep='\t',index=False)
missing=meta.groupby('dataset').apply(lambda x:pd.Series({c:int((x[c]=='').sum()) for c in ['patient_id','regimen','sampling_time','response','procedure','stage','preservation']}),include_groups=False)
missing.to_csv(ROOT/'06_tables/metadata_missingness.tsv',sep='\t')
report='''# 元数据来源与分组审计（表达对象映射前）

GSE205013 的 27 个 GSM 对应 27 个 Pxx 患者。原发灶切除 10 人（未治疗 6、已治疗 4），原发灶穿刺 7 人（未治疗 5、已治疗 2），肝转移穿刺 10 人。主比较以原发灶 17 人为登记候选，最终有效人数取决于细胞质量和恶性证据。

GEO 中肝转移为未治疗 9、治疗 1。论文补充表 2 将 P11 列为接受 FOLFIRINOX、疗效 PD，而 GEO 对该肝转移样本标为未治疗。不能确认测序取材与治疗的时序，故统一表将 P11 状态设为 Unknown，保留双方原始证据。该冲突不改变原发灶主比较的登记人数。P10 接受 FOLFIRINOX/SBRT；主分析和剔除 P10 的敏感性结果分开报告。FOLFIRINOX 与吉西他滨/白蛋白紫杉醇序贯患者不归作单一方案。

GSE202051 登记 43 个 PDAC 患者 ID：未治疗 18，已治疗 25。标准 CRT 为 14；CRTL 为 5；CRTN 和 CRTLN 各 1；其余方案 4。CRT 包括 FOLFIRINOX 和放疗，CRTL 另加氯沙坦；CRTN/CRTLN 包含纳武利尤单抗。技术重复 GSM 只作为重复文库，不扩充样本量。临床补充表中的 U19–U21、T26–T27 标明仅有 DSP，无单核测序，排除在本阶段单核样本量之外。

GEO 另有 6 个非恶性组织 ID 和 2 个类器官实验记录；类器官不得视作额外患者。文件说明“45 specimens”与登记的组织 ID 总数不同，将以 h5ad 实际 obs 和临床映射单独核对，当前不宣称全部可用。010nuc 与 010orgCRT 均是类器官对象；仅下载 totaldata-final-toshare 作为组织对象。数字空间实验为 NanoString GeoMx DSP，不标为 Visium。

未知临床字段为空并由 unknown_fields/metadata_missingness 标记。治疗时长保留在独立字段，不能替代测序取材距治疗结束的间隔。Hwang 表中的病理 ypT 与未治疗 T 分期不是等价基线协变量，验证主比较按真实方案分层。

队列分别来自 NYU 和 MGH 的独立原始研究，未见共同患者标识或来源复用说明；匿名 ID 无法执行遗传身份层面的重复排查。独立验证不合并细胞作统计重复。

完整逐字段证据：01_metadata/metadata_field_provenance.tsv。登记文库/标本/患者统计：06_tables/registered_samples_vs_patients.tsv。GEO 原始记录、论文全文、临床补充表及 SHA256 清单均已保存。表达对象映射完成后才可判定 METADATA_AUDIT_PASS。
'''
(ROOT/'08_reports/metadata_source_audit.md').write_text(report)
print(counts.to_string(index=False));print('METADATA_SOURCE_REVIEW_COMPLETE')
