root <- '/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis'
suppressPackageStartupMessages({library(edgeR); library(data.table)})
files <- list.files(file.path(root,'05_analysis'), pattern='_raw_pseudobulk.tsv.gz$', full.names=TRUE)
audit <- list()
for (file in files) {
  stem <- sub('_raw_pseudobulk.tsv.gz$','',basename(file))
  ds <- substr(stem,1,9)
  group <- substring(stem,11)
  if (!group %in% c('Malignant','EpiAll','MalStrictQC')) next
  dat <- as.data.frame(fread(file))
  genes <- dat[[1]]
  counts <- as.matrix(dat[,-1,drop=FALSE]); rownames(counts) <- genes
  meta <- as.data.frame(fread(file.path(root,'05_analysis',paste0(stem,'_pseudobulk_metadata.tsv'))))
  stopifnot(identical(colnames(counts),meta$patient_id),all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
  comparisons <- if (ds=='GSE205013') c('Primary','Resection','NoRT','StageAdj') else c('CRT','CRTL','AllTx')
  if (group!='Malignant') comparisons <- if (ds=='GSE205013') 'Primary' else 'CRT'
  for (comparison in comparisons) {
    sel <- rep(TRUE,nrow(meta))
    if (comparison=='Resection') sel <- meta$procedure=='Resection'
    if (comparison=='NoRT') sel <- meta$patient_id!='P10'
    if (comparison=='CRT') sel <- meta$regimen_group %in% c('Untreated','CRT')
    if (comparison=='CRTL') sel <- meta$regimen_group %in% c('Untreated','CRTL')
    md <- meta[sel,,drop=FALSE]; ct <- counts[,sel,drop=FALSE]
    md$treatment <- factor(md$treatment_status,levels=c('Untreated','Treated'))
    nn <- table(md$treatment)
    if (min(nn)<2) next
    terms <- 'treatment'
    if (ds=='GSE205013' && length(unique(md$procedure))>1) terms <- c('procedure',terms)
    if (comparison=='StageAdj' && length(unique(md$stage))>1) terms <- c('stage',terms)
    design <- model.matrix(reformulate(terms),md)
    requested_rank <- qr(design)$rank
    requested_n_parameters <- ncol(design)
    design_note <- ''
    if (comparison=='StageAdj' && qr(design)$rank<ncol(design)) {
      design <- model.matrix(~stage+treatment,md)
      design_note <- 'Procedure and stage collinear; stage replaces procedure in this sensitivity model'
    }
    key <- paste(stem,comparison,sep='_')
    ar <- data.frame(dataset=ds,cell_group=group,comparison=comparison,n_untreated=unname(nn['Untreated']),n_treated=unname(nn['Treated']),rank=qr(design)$rank,n_parameters=ncol(design),residual_df=nrow(md)-ncol(design),requested_rank=requested_rank,requested_n_parameters=requested_n_parameters,design=paste(colnames(design),collapse=';'),design_note=design_note,status='')
    if (qr(design)$rank<ncol(design) || nrow(md)-ncol(design)<3) {
      ar$status <- 'Not_estimable';audit[[key]] <- ar;next
    }
    y <- DGEList(ct)
    keep <- filterByExpr(y,design=design)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y,method='TMM')
    y <- estimateDisp(y,design=design,robust=TRUE)
    fit <- glmQLFit(y,design=design,robust=TRUE)
    test <- glmQLFTest(fit,coef=which(colnames(design)=='treatmentTreated'))
    result <- topTags(test,n=Inf,sort.by='none')$table
    result$gene <- rownames(result)
    result$dataset <- ds;result$cell_group <- group;result$comparison <- comparison
    result$n_untreated <- unname(nn['Untreated']);result$n_treated <- unname(nn['Treated'])
    result$inference <- if (min(nn)<10) 'Exploratory' else 'Patient_level_observational'
    fwrite(result,file.path(root,'06_tables',paste0(key,'_edgeR.tsv.gz')),sep='\t')
    norm <- data.frame(patient_id=md$patient_id,library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors)
    fwrite(norm,file.path(root,'05_analysis',paste0(key,'_TMM.tsv')),sep='\t')
    saveRDS(list(y=y,design=design,metadata=md,fit=fit,test=test),file.path(root,'05_analysis',paste0(key,'_edgeR.rds')))
    ar$status <- 'PASS';ar$n_genes_tested <- nrow(result);ar$n_FDR05 <- sum(result$FDR<.05);audit[[key]] <- ar
    cat('DE_COMPARISON_COMPLETE',key,'patients',nrow(md),'genes',nrow(result),'FDR05',sum(result$FDR<.05),'\n')
  }
}
if (!length(audit)) stop('No eligible raw-count comparisons')
fwrite(rbindlist(audit,fill=TRUE),file.path(root,'06_tables/pseudobulk_DE_audit.tsv'),sep='\t')
cat('PSEUDOBULK_DE_COMPLETE\n')
