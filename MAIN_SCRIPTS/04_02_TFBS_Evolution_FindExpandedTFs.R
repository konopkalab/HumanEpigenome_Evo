rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(Seurat)
library(rio)
library(GenomicRanges)
library(rphast)
library(ape)
library(parallel)
library(Biostrings)
library(ggpubr)
library(seqinr)
library(phangorn)
library(msa)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(readr)
library(GeneOverlap)
source('Functions.R')

####
## MERGE AND PLOT OVERALL PATTERNS
####

fls = list.files(path = 'Motif_Evolution', pattern = 'motifEvo.RDS', full.names = T)
motres = lapply(fls, function(x){read_rds(x)}) %>% do.call(rbind, .)

write_rds(motres, 'Adult_AllMotifEvo.RDS')

# Plot gain vs loss
toplot = motres %>% group_by(node) %>% group_keys %>% as.data.frame
toplot$size = motres %>% group_by(node) %>% group_size
toplot$gainSize = motres %>% filter(type == 'Gain') %>% group_by(type, node) %>% group_size
toplot$lossSize = motres %>% filter(type == 'Loss') %>% group_by(type, node) %>% group_size
toplot$gainPerc = toplot$gainSize / toplot$size * 100
toplot$lossPerc = toplot$lossSize / toplot$size * 100

toplotM = melt(toplot)
toplotM = toplotM[toplotM$variable %in% c('gainPerc', 'lossPerc'),]
toplotM$variable2 = ifelse(toplotM$variable == 'gainPerc', 'Gain', 'Loss')
toplotM$node = factor(toplotM$node, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))

pdf('Adult_GainLoss_Motifs.pdf', width = 7, height = 5)
ggbarplot(toplotM, x = 'variable2', y = 'value', fill = 'variable2',
		color = 'variable2',  position = position_dodge(0.9)) +
theme_classic() + xlab('') + ylab('Percentage') +
theme(text = element_text(size=20, face = 'bold')) +
facet_wrap(~node, nrow = 1) + NoLegend() +
rotate_x_text(45)
dev.off()

# Plot number of motif changes normalized by number of substitutions
allsubs = read_rds('AdultSubstitutions.RDS')
allsubs$count = 1
allsubs2 = aggregate(count ~ nodes, data = allsubs, FUN = sum)
rownames(allsubs2) = allsubs2$nodes
allsubs2 = allsubs2[c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'),]
allsubs2$evoDist = c(6,2,8,4,13)
allsubs2$normCount = allsubs2$count / allsubs2$evoDist

motres$count = 1
motres2 = aggregate(count ~ node, data = motres, FUN = sum)
rownames(motres2) = motres2$node
motres2 = motres2[c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'),]
motres2$evoDist = c(6,2,8,4,13)
motres2$normCount1 = motres2$count / motres2$evoDist
motres2$normCount2 = motres2$count / allsubs2$count

pdf('Motif_changes_per_substitution_changes.pdf')
ggbarplot(motres2, x = 'node', y = 'normCount2', color = 'node', fill = 'node', ylim = c(0,1)) +
ylab('Motif changes / Substitution changes') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()


####
## Count Gains and Losses per node
####

#cls = read_rds('~/workdir/pr5/02_GroupCREs/ADULT/AdultCREs_SignTested.RDS')
alltfs = unique(motres$tf)

nodes = c("Human", "HC", "HCGo", "Great_Ape", "Ape")
dfL = list()
for(i in 1:length(alltfs)){
	
	tmp = motres[motres$tf == alltfs[i], ]
	tmp$count = 1
	tmp2 = aggregate(count ~ type + node, tmp, FUN = sum)
	tmp3 = dcast(tmp2, node ~ type, value.var = 'count')
	tmp3[is.na(tmp3)] = 0
	tmp3 = as.data.frame(tmp3)
	tmp3$tf = alltfs[i]

	dfL[[i]] = tmp3
	if(i%%100==0){print(i)}
}

finaldf = do.call(rbind, dfL)
finaldf$ratio = finaldf$Gain / finaldf$Loss

saveRDS(finaldf, 'Adult_AllMotifEvoCounts.RDS')

####
## Correlations of motif gain/loss between nodes
####

nodes = unique(finaldf$node)
ratiosL = lapply(nodes, function(x){finaldf[finaldf$node == x, 'ratio']})
names(ratiosL) = nodes
ratios_df = do.call(cbind, ratiosL) %>% as.data.frame

corsL = list()
for(i in 1:length(nodes)){
	corsL[[i]] = sapply(nodes, function(x){cor.test(ratios_df[,nodes[i]], ratios_df[,x])$estimate})
}
names(corsL) = nodes
cors_df = do.call(rbind, corsL)

toplot = melt(cors_df)
toplot$label = round(toplot$value, digits = 2)
toplot$Var1 = factor(toplot$Var1, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))
toplot$Var2 = gsub('.cor', '', toplot$Var2)
toplot$Var2 = factor(toplot$Var2, levels = rev(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')))

pdf('GainLoss_Cor_AcrossMotifs.pdf', width = 8, height = 6)
ggplot(toplot,aes(x=Var1,y=Var2,fill=value))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0.2, low = 'blue', high = 'red', limits=c(0,1))+
  theme_classic()+
  geom_text(aes(label = label), vjust = 0.8, colour = "black", size = 6) +
  theme(text = element_text(size=20, face = 'bold'))
dev.off()

####
## TFs that are significantly expanded in one node compared to all other nodes and overall background
####

# Test node to other nodes
nodes = c("Human", "HC", "HCGo", "Great_Ape", "Ape")
dfL2 = list()
for(j in 1:length(nodes)){

	tfs = unique(finaldf$tf)
	dfL1 = list()
	for(i in 1:length(tfs)){

		test = finaldf[finaldf$tf == tfs[i] & finaldf$node == nodes[j],]
		bcg = finaldf[finaldf$tf == tfs[i], ]

		testG = test['Gain'] %>% as.numeric
		bcgG = sum(bcg$Gain)

		testSum = sum(test[c('Gain', 'Loss')])
		bcgSum = sum(sum(bcg$Gain), sum(bcg$Loss))

		pval = prop.test(c(testG, bcgG), c(testSum, bcgSum))$p.value
		oddsR = (testG / testSum) / (bcgG / bcgSum)

		dfL1[[i]] = data.frame(tf = tfs[i], pval = pval, oddsRatio = oddsR, node = nodes[j])
	}

	df1 = do.call(rbind, dfL1)
	df1$FDR = p.adjust(df1$pval, method = 'fdr')
	dfL2[[j]] = df1
	print(j)
}

dfNode = do.call(rbind, dfL2)

# Test to all background
nodes = c("Human", "HC", "HCGo", "Great_Ape", "Ape")
bcgG = sum(finaldf$Gain)
bcgSum = sum(sum(finaldf$Gain), sum(finaldf$Loss))

dfL2 = list()
for(j in 1:length(nodes)){

	tfs = unique(finaldf$tf)
	dfL1 = list()
	for(i in 1:length(tfs)){

		test = finaldf[finaldf$tf == tfs[i] & finaldf$node == nodes[j],]
		testG = test['Gain'] %>% as.numeric
		testSum = sum(test[c('Gain', 'Loss')])

		pval = prop.test(c(testG, bcgG), c(testSum, bcgSum))$p.value
		oddsR = (testG / testSum) / (bcgG / bcgSum)

		dfL1[[i]] = data.frame(tf = tfs[i],pval = pval, oddsRatio = oddsR, node = nodes[j])
	}

	df1 = do.call(rbind, dfL1)
	df1$FDR = p.adjust(df1$pval, method = 'fdr')
	dfL2[[j]] = df1
	print(j)
}

dfBcg = do.call(rbind, dfL2)

dfNode$FDR_Bcg = dfBcg$FDR
dfNode$oddsRatioBcg = dfBcg$oddsRatio

saveRDS(dfNode, 'Adult_RawMotifAll.RDS')

enrDF = dfNode[dfNode$FDR < 0.01 & dfNode$FDR_Bcg < 0.01 & dfNode$oddsRatio > 1 & dfNode$oddsRatioBcg > 1,]
depDF = dfNode[dfNode$FDR < 0.01 & dfNode$FDR_Bcg < 0.01 & dfNode$oddsRatio < 1 & dfNode$oddsRatioBcg < 1,]

saveRDS(enrDF, 'Adult_RawMotifEnrich.RDS')
saveRDS(depDF, 'Adult_RawMotifDepleted.RDS')

####
## Filter by expression / accessibility
###

hres = enrDF[enrDF$node == 'Human', 'tf']
hcres = enrDF[enrDF$node == 'HC', 'tf']
aperes = enrDF[enrDF$node == 'Ape', 'tf']

# RNA expressed genes
rnaSeurat = read_rds('Caglayan2023_snRNAseq_Seurat_Human.RDS')

ctypes = unique(rnaSeurat$CellType)
expGnsL = list()
for(i in 1:length(ctypes)){

	tmpSeurat = subset(rnaSeurat, subset = CellType %in% ctypes[i])
	rnaMat = tmpSeurat@assays$RNA@counts
	rnaMat[rnaMat > 0] = 1
	expGnsL[[i]] = rnaMat[rowMeans(rnaMat) > 0.25,] %>% rownames
	print(i)
}

names(expGnsL) = unique(rnaSeurat$CellType)
saveRDS(expGnsL, 'Adult_expGnsL.RDS')

expGns2 = unlist(expGnsL) %>% unique
saveRDS(expGns2, 'Adult_expGns2.RDS')

# ATAC open chromatin genes
atacSeurat = read_rds('Caglayan2023_snATACseq_Seurat_Human.RDS')

ctypes = unique(atacSeurat$newannot)
accGnsL = list()
for(i in 1:length(ctypes)){

	tmpSeurat = subset(atacSeurat, subset = newannot %in% ctypes[i])
	atacMat = tmpSeurat@assays$RNA@counts
	atacMat[atacMat > 0] = 1
	accGnsL[[i]] = atacMat[rowMeans(atacMat) > 0.25,] %>% rownames
	print(i)
}

names(accGnsL) = unique(atacSeurat$newannot)
saveRDS(accGnsL, 'Adult_accGnsL.RDS')

accGns2 = unlist(accGnsL) %>% unique
saveRDS(accGns2, 'Adult_accGns2.RDS')

# Filter by expression or accessibility
hres2 = hres[hres %in% expGns2 | hres %in% accGns2] %>% sort
hcres2 = hcres[hcres %in% expGns2 | hcres %in% accGns2] %>% sort
aperes2 = aperes[aperes %in% expGns2 | aperes %in% accGns2] %>% sort

# Keep only the expressed / accessible TFs
enrDF2 = enrDF[enrDF$tf %in% Reduce(union, list(hres2, hcres2, aperes2)),]
depDF2 = depDF[depDF$tf %in% Reduce(union, list(hres2, hcres2, aperes2)),]

saveRDS(enrDF2, 'Adult_RawMotifEnrich_ExpressedAccessible.RDS')
saveRDS(depDF2, 'Adult_RawMotifDepleted_ExpressedAccessible.RDS')

# Plot number of TFs
tmp = enrDF2
tmp$node = factor(tmp$node, levels = nodes)
toplot = table(tmp$node) %>% as.data.frame

pdf('Adult_number_of_significantly_expanded_motifs.pdf')
ggbarplot(toplot, x = 'Var1', y = 'Freq', color = 'gray', fill = 'gray') +
ylab('Number of significantly expanded TFs') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

tmp = depDF2
tmp$node = factor(tmp$node, levels = nodes)
toplot = table(tmp$node) %>% as.data.frame

pdf('Adult_number_of_significantly_depleted_motifs.pdf')
ggbarplot(toplot, x = 'Var1', y = 'Freq', color = 'gray', fill = 'gray') +
ylab('Number of significantly expanded TFs') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

# Plot enrichment stats of all TFs

# HUMAN #
enrDF2_human = enrDF2[enrDF2$node == 'Human', 'tf']
toplot = dfNode[dfNode$tf %in% enrDF2_human,]
toplot$log10FDR = -log10(toplot$FDR)
toplot$is_sign = ifelse(toplot$FDR_Bcg < 0.05 & toplot$FDR < 0.05, '**',
			ifelse(!(toplot$FDR_Bcg < 0.05) & toplot$FDR < 0.05, '*', ''))

toplot$node = factor(toplot$node, levels = rev(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')))

pdf('Adult_Human_AllDivTFs_EnrichStats.pdf', width = 20, height = 6)
ggscatter(toplot, x = 'tf', y = 'node', color = 'oddsRatio', size = 'log10FDR') +
  scale_size_continuous(range = c(4,15)) +
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  geom_text(aes(label = is_sign), vjust = 0.8, colour = "darkgreen", size = 12 ) +
  #facet_wrap(~tf) +
  rotate_x_text(45)
dev.off()

# APE #
enrDF2_human = enrDF2[enrDF2$node == 'Ape', 'tf']
toplot = dfNode[dfNode$tf %in% enrDF2_human,]
toplot$log10FDR = -log10(toplot$FDR)

toplot$is_sign = ifelse(toplot$FDR_Bcg < 0.05 & toplot$FDR < 0.05, '**',
			ifelse(!(toplot$FDR_Bcg < 0.05) & toplot$FDR < 0.05, '*', ''))

toplot$node = factor(toplot$node, levels = rev(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')))

pdf('Adult_ape_AllDivTFs_EnrichStats.pdf', width = 20, height = 6)
ggscatter(toplot, x = 'tf', y = 'node', color = 'oddsRatio', size = 'log10FDR') +
  scale_size_continuous(range = c(8,15)) +
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  geom_text(aes(label = is_sign), vjust = 0.8, colour = "darkgreen", size = 12 ) +
  #facet_wrap(~tf) +
  rotate_x_text(45)
dev.off()

####
## Self Enrichment
###
library(GeneOverlap)
bcg = length(tfs)

# Prepare data
tmpL = split(enrDF2, enrDF2$node)
expL = lapply(tmpL, function(x){x$tf})
names(expL) = names(tmpL)
expL = expL[c('Human', 'HC', 'Ape')]

resgom = newGOM(expL, expL, genome.size = bcg)
pvalmat = getMatrix(resgom, name="pval")
pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
pvalmelt = melt(pvalmat, value.name = 'pval')
pvalmelt$log10_FDR = -log10(pvalmelt$pval)
pvalmelt = pvalmelt[pvalmelt$Var1 != pvalmelt$Var2,]

oddsmat = getMatrix(resgom, name="odds.ratio")
oddsmelt = melt(oddsmat, value.name = 'odds.ratio')
oddsmelt = oddsmelt[oddsmelt$Var1 != oddsmelt$Var2,]

jacmat = getMatrix(resgom, name="Jaccard")
colnames(jacmat) = rownames(jacmat)
jacmelt = melt(jacmat, value.name = 'jaccard')
jacmelt = jacmelt[jacmelt$Var1 != jacmelt$Var2,]

resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio, JAC = jacmelt$jaccard)
resdf$log10_round = round(resdf$log10_FDR, digits = 2)
resdf$OR_round = round(resdf$OR, digits = 2)
resdf$JAC_round = round(resdf$JAC, digits = 2)
resdf$fdr_Sc = formatC(resdf$pval, format = "e", digits = 2)

pdf('Adult_Similarity_of_significantly_expanded_motifs.pdf', width = 7, height = 5)
ggscatter(resdf, x = 'Var1', y = 'Var2', color = 'log10_FDR', size = 'OR') +
  scale_size_continuous(range = c(3,10)) +
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 2, low = 'blue', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(45)
dev.off()


####
## Find TFs with significant difference of expansion between Human - HC
###

tmpTFs = enrDF2[enrDF2$node %in% c('Human', 'HC'), 'tf'] %>% unique
dfL = list()
for(i in 1:length(tmpTFs)){

	test = finaldf[finaldf$tf == tfs[i] & finaldf$node == 'Human',]
	bcg = finaldf[finaldf$tf == tfs[i] & finaldf$node == 'HC',]

	testG = test['Gain'] %>% as.numeric
	bcgG = bcg['Gain'] %>% as.numeric

	testSum = sum(test[c('Gain', 'Loss')])
	bcgSum = sum(bcg[c('Gain', 'Loss')])

	pval = prop.test(c(testG, bcgG), c(testSum, bcgSum))$p.value
	oddsR = (testG / testSum) / (bcgG / bcgSum)

	dfL[[i]] = data.frame(tf = tmpTFs[i], pval = pval, oddsRatio = oddsR)
}

diffHomin = do.call(rbind, dfL)
diffHomin$FDR = p.adjust(diffHomin$pval, method = 'fdr')

diffHomin[diffHomin$FDR < 0.05,] # No significant difference


####
## BARPLOTS
####

# Load motif evolution results
motres = read_rds('Adult_AllMotifEvo.RDS')
adult_enrDF2 = read_rds('Adult_RawMotifEnrich_ExpressedAccessible.RDS')
fetal_enrDF2 = read_rds('Adult_RawMotifEnrich_ExpressedAccessible.RDS')
finaldf = read_rds('Adult_AllMotifEvoCounts.RDS')
finaldf$node = factor(finaldf$node, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))

adult_enrDF2$dev = 'Adult'
fetal_enrDF2$dev = 'Fetal'
enrDF2 = rbind(adult_enrDF2, fetal_enrDF2)

# TFs to plot
tfs = enrDF2[enrDF2$node == 'Human', 'tf'] %>% unique
toplot = finaldf[finaldf$tf %in% tfs,]

# Background gain/loss ratio
bcg = sum(finaldf$Gain) / sum(finaldf$Loss)

adult_depDF2 = readRDS('Adult_RawMotifDepleted_ExpressedAccessible.RDS')
adult_enrDF2$id = paste0(adult_enrDF2$tf, '_', adult_enrDF2$node)
adult_depDF2$id = paste0(adult_depDF2$tf, '_', adult_depDF2$node)

toplot$id = paste0(toplot$tf, '_', toplot$node)
toplot$is_sign = ifelse(toplot$id %in% adult_enrDF2$id , 'Sign',
		ifelse(toplot$id %in% adult_depDF2$id, 'Sign_Dep', 'NS'))

pdf('AllHumanExpanded_TFs_ADULT_GainLoss.pdf', width = 18, height = 14)
ggbarplot(toplot, x = 'node', y = 'ratio', fill = 'lightblue',
		color = 'is_sign', ylim = c(0,2.2), size = 1, palette = c('lightblue', 'red', 'black')) +
theme_classic() +
ylab('Gain / Loss') + xlab('') +
theme(text = element_text(size=25, face = 'bold'), plot.title = element_text(hjust = 0.5)) +
facet_wrap(~tf) +
geom_hline(yintercept = bcg, linetype = 'dashed', color = 'red') +
rotate_x_text(45)
dev.off()

####
## Cell Type Enrichments
###

# Load motif evolution results
motres = read_rds('Adult_AllMotifEvo.RDS')
adult_enrDF2 = readRDS('Adult_RawMotifEnrich_ExpressedAccessible.RDS')
fetal_enrDF2 = readRDS('Adult_RawMotifEnrich_ExpressedAccessible.RDS')

adult_enrDF2$dev = 'Adult'
fetal_enrDF2$dev = 'Fetal'
enrDF2 = rbind(adult_enrDF2, fetal_enrDF2)

# Peak to gene links
gene_link = read_rds('DATASETS/Ma_LinkPeaksToGenes_FINAL.RDS')
gene_link$peak2 = sub('-', ':', gene_link$peak)

## HUMAN ##

# The node to run cell type enrichments
the_node = 'Ape'
tfs = enrDF2[enrDF2$node == the_node, 'tf'] %>% unique

# Motifs
tfs = enrDF2[enrDF2$node == the_node, 'tf'] %>% unique

# Cell type markers
ctmarks = read_rds('DATASETS/Adult_PseudoBulk_DEGs_MAJOR_MARKERS.RDS')
ctmarks = ctmarks[ctmarks$CellType != 'Others',]
tmpL = split(ctmarks, ctmarks$CellType)
ctL = lapply(tmpL, function(x){x[, 'Gene']})
ctL = ctL[c('Excitatory', 'Inhibitory', 'OPC', 'MOL', 'Astrocyte', 'Microglia')]

# Find gene linkes to tf target gains / losses
linksL = list()
links_saveL = list()
for(i in 1:length(tfs)){

	motres_gain = motres[motres$tf %in% tfs[i] & motres$node %in% the_node & motres$type == 'Gain',]
	motres_loss = motres[motres$tf %in% tfs[i] & motres$node %in% the_node & motres$type == 'Loss',]

	gene_link_gain = gene_link[gene_link$peak2 %in% motres_gain$CRE, 'gene'] %>% unique
	gene_link_loss = gene_link[gene_link$peak2 %in% motres_loss$CRE, 'gene'] %>% unique

	vars = c(rep('Gain', length(gene_link_gain)), rep('Loss', length(gene_link_loss)))
	linksL[[i]] = data.frame(vals = c(gene_link_gain, gene_link_loss), vars = vars, tf = tfs[i])

	# Save all gained/lost targets and their links
	q1 = gene_link[gene_link$peak2 %in% motres_gain$CRE, ] %>%
			add_column(type = 'gain') %>% add_column(tf = tfs[i])

	q2 = gene_link[gene_link$peak2 %in% motres_loss$CRE, ] %>%
			add_column(type = 'loss') %>% add_column(tf = tfs[i])

	links_saveL[[i]] = rbind(q1, q2)

}

linksdf = do.call(rbind, linksL)
links_savedf = do.call(rbind, links_saveL)

# Combine MEF2 related TFs
linksdf$tf = gsub('MEF2.*', 'MEF2', linksdf$tf)
links_savedf$tf = gsub('MEF2.*', 'MEF2', links_savedf$tf)

# Combine HD-LIM TFs
linksdf$tf = gsub('LHX6|LMX1A|LMX1B', 'HD-LIM', linksdf$tf)
links_savedf$tf = gsub('LHX6|LMX1A|LMX1B', 'HD-LIM', links_savedf$tf)

linksdf$id = paste0(linksdf$tf, '_', linksdf$vars)

linksdfL = split(linksdf, linksdf$id)
nms = names(linksdfL)
linksdfL = lapply(linksdfL, function(x){x$vals})
names(linksdfL) = nms

# Save
write_rds(links_savedf, paste0(the_node, '_Motif_All_GeneLinks.RDS'))

# Background
alAggMat = read_rds('DATASETS/Adult_RNA_pseudobulk_perBroadCellType_normalized.RDS')
bcg_1 = rownames(alAggMat)
bcg_2 = unique(gene_link$gene)
bcg = length(union(bcg_1, bcg_2))

# Enrichment in cell type marker genes
toplot = geneOvEnr(linksdfL, ctL, bcg, plot = F, label_type = 'symbol', hg = 10, wd = 14)
toplot$sign_label = ifelse(toplot$FDR < 0.05 & toplot$value > 1.3, '*', '')
toplot$value = ifelse(toplot$value > 1.5, 1.5, toplot$value)

pdf(paste0('ADULT_', the_node, '_AllDivTFs_CellTypeEnrich.pdf'), height =14, width = 14)
ggscatter(toplot, x = 'Var2', y = 'type', color = 'value', size = 'log10FDR') +
  scale_size_continuous(range = c(6,12)) +
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red', limits = c(0.5,1.5)) +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  geom_text(aes(label = sign_label), vjust = 0.8, colour = "darkgreen", size = 12 ) +
  facet_wrap(~tf, ncol = 5) +
  rotate_x_text(90)
dev.off()

# Highlight target genes
test = links_savedf[links_savedf$tf == 'CLOCK' & links_savedf$type == 'gain', 'gene'] %>% table %>% sort %>% tail(20)
link_bcg = table(gene_link$gene)
sort(test/link_bcg[names(test)], decreasing = T)

## APE ##

# The node to run cell type enrichments
the_node = 'Ape'

# Motifs
tfs = enrDF2[enrDF2$node == the_node, 'tf']

# Cell type markers
ctmarks = read_rds('DATASETS/Adult_PseudoBulk_DEGs_MAJOR_MARKERS.RDS')
ctmarks = ctmarks[ctmarks$CellType != 'Others',]
tmpL = split(ctmarks, ctmarks$CellType)
ctL = lapply(tmpL, function(x){x[, 'Gene']})
ctL = ctL[c('Excitatory', 'Inhibitory', 'OPC', 'MOL', 'Astrocyte', 'Microglia')]

# Find gene linkes to tf target gains / losses
linksL = list()
links_saveL = list()
for(i in 1:length(tfs)){

	motres_gain = motres[motres$tf %in% tfs[i] & motres$node %in% the_node & motres$type == 'Gain',]
	motres_loss = motres[motres$tf %in% tfs[i] & motres$node %in% the_node & motres$type == 'Loss',]

	gene_link_gain = gene_link[gene_link$peak2 %in% motres_gain$CRE, 'gene'] %>% unique
	gene_link_loss = gene_link[gene_link$peak2 %in% motres_loss$CRE, 'gene'] %>% unique

	vars = c(rep('Gain', length(gene_link_gain)), rep('Loss', length(gene_link_loss)))
	linksL[[i]] = data.frame(vals = c(gene_link_gain, gene_link_loss), vars = vars, tf = tfs[i])

	# Save all gained/lost targets and their links
	q1 = gene_link[gene_link$peak2 %in% motres_gain$CRE, ] %>%
			add_column(type = 'gain') %>% add_column(tf = tfs[i])

	q2 = gene_link[gene_link$peak2 %in% motres_loss$CRE, ] %>%
			add_column(type = 'loss') %>% add_column(tf = tfs[i])

	links_saveL[[i]] = rbind(q1, q2)

}

linksdf = do.call(rbind, linksL)
links_savedf = do.call(rbind, links_saveL)

# Combine MEF2 TFs
linksdf$tf = gsub('MEF2[A-Z]', 'MEF2', linksdf$tf)
links_savedf$tf = gsub('MEF2[A-Z]', 'MEF2', links_savedf$tf)

linksdf$id = paste0(linksdf$tf, '_', linksdf$vars)

linksdfL = split(linksdf, linksdf$id)
nms = names(linksdfL)
linksdfL = lapply(linksdfL, function(x){x$vals})
names(linksdfL) = nms

# Save
write_rds(links_savedf, paste0(the_node, '_Motif_All_GeneLinks.RDS'))

# Background
alAggMat = read_rds('DATASETS/Adult_RNA_pseudobulk_perBroadCellType_normalized.RDS')
bcg_1 = rownames(alAggMat)
bcg_2 = unique(gene_link$gene)
bcg = length(union(bcg_1, bcg_2))

# Enrichment in cell type marker genes
toplot = geneOvEnr(linksdfL, ctL, bcg, plot = F, label_type = 'symbol', hg = 10, wd = 14)
toplot$sign_label = ifelse(toplot$FDR < 0.05 & toplot$value > 1.3, '*', '')
toplot$value = ifelse(toplot$value > 1.5, 1.5, toplot$value)

pdf(paste0('ADULT_', the_node, '_AllDivTFs_CellTypeEnrich.pdf'), height = 11, width = 14)
ggscatter(toplot, x = 'Var2', y = 'type', color = 'value', size = 'log10FDR') +
  scale_size_continuous(range = c(4,12)) +
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red', limits = c(0.5,1.5)) +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  geom_text(aes(label = sign_label), vjust = 0.8, colour = "darkgreen", size = 12 ) +
  facet_wrap(~tf, ncol = 5) +
  rotate_x_text(90)
dev.off()

# Highlight target genes
test = links_savedf[links_savedf$tf == 'MEF2' & links_savedf$type == 'gain', 'gene'] %>% table %>% sort %>% tail(20)
link_bcg = table(gene_link$gene)
sort(test/link_bcg[names(test)], decreasing = T)


