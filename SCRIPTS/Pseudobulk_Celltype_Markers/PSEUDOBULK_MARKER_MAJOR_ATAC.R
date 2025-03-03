library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(rio)
library(scran)
library(scater)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(DESeq2)

outdir = 'OUT_DATA/Pseudobulk_DEGs_DARs'
dir.create(outdir, showWarnings=F)

# Load objects
seurObj = read_rds("IN_DATA/snATACseq.RDS")

seurObj$tot_acc_cicero = log2(seurObj$nCount_RNA)
meta = seurObj[[]]

# Add sample metadata
smeta = import('IN_DATA/atac_metadata.xlsx')
smeta$Species = NULL
meta = cbind(meta, smeta[match(meta$Sample, smeta$Sample_id), ])
meta$lib_batch = factor(meta$lib_batch)
seurObj@meta.data = meta

# Normalize and log transform RNA
seurObj = NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000, assay = 'RNA')

seurObj$CellType = seurObj$newannot
seurObj$CellType = gsub("L4-5_RORB_2-L4-6_RORB_1", "L4-6_RORB_3", seurObj$CellType)

# Convert to major cell types
meta = seurObj@meta.data
meta$Major = meta$CellType
meta[grepl('^L[0-9]', meta$CellType), 'Major'] = 'Excitatory'
meta[grepl('Upper|SST|PVALB|VIP|LAMP5', meta$CellType), 'Major'] = 'Inhibitory'
seurObj@meta.data = meta


# All cell types
celltypes = unique(seurObj$Major)

# Convert seurat to sce
DefaultAssay(seurObj) = 'peaks'
allSCE = as.SingleCellExperiment(seurObj)

# Pseudobulk SCE assay
allGroups = colData(allSCE)[, c("Major", "Sample_id", "human_age", "sex")]
allPseudoSCE = sumCountsAcrossCells(allSCE, allGroups)
allGroupsPseudo = colData(allPseudoSCE)[, c("Major", "Sample_id", "human_age", "sex")]
allAggMat = allPseudoSCE@assays@data$sum
colnames(allAggMat) = allGroupsPseudo[['Major']]

# Pseudobulk SCE assay
allGroupsPseudo$Major = factor(allGroupsPseudo$Major)
allGroupsPseudo$Sample_id = factor(allGroupsPseudo$Sample_id)
allGroupsPseudo$sex = factor(allGroupsPseudo$sex)

resMarkL = list()
for(i in 1:length(celltypes)){

	print(celltypes[i])

	allGroupsPseudo$ToTest = ifelse(allGroupsPseudo$Major == celltypes[i], celltypes[i], 'Others')
	allGroupsPseudo$ToTest = factor(allGroupsPseudo$ToTest, levels = c(celltypes[i], 'Others'))

	# Run DGE
	dgel = DGEList(counts = allAggMat)
	dgel = calcNormFactors(dgel)
	design = model.matrix(~ToTest + human_age + sex, data = allGroupsPseudo)
	dgel = estimateDisp(dgel, design)

	fit = glmFit(dgel,design)
	lrt = glmLRT(fit, coef = paste0('ToTestOthers'))
	res = topTags(lrt, n = nrow(allAggMat), sort.by = 'none') %>% as.data.frame

	resMarkL[[i]] = res[res$FDR < 0.05 & res$logFC < -0.3,]
	resMarkL[[i]]$CellType = celltypes[i]
	resMarkL[[i]]$Gene = rownames(resMarkL[[i]])
	
	print(i)
}

resMarkDF = do.call(rbind, resMarkL)

saveRDS(resMarkDF, paste0(outdir, '/PseudoBulk_DARs_MAJOR_MARKERS.RDS'))

# Aggregate across samples and scale for plotting
allGroups = colData(allSCE)[, c("Major")]
allPseudoSCE = sumCountsAcrossCells(allSCE, allGroups)
allGroupsPseudo = colData(allPseudoSCE)[, c("ids")]
allAggMat = allPseudoSCE@assays@data$sum
allAggMat = log2(cpm(allAggMat) + 1)
write_rds(allAggMat, paste0(outdir, '/pseudobulk_perBroadCellType_normalized_atac.RDS'))

tmp = apply(allAggMat, 1, FUN = scale)
tmp = t(tmp)
colnames(tmp) = colnames(allAggMat)
write_rds(tmp, paste0(outdir, '/pseudobulk_perBroadCellType_scaled_atac.RDS'))


