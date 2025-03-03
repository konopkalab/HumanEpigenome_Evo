require(plyr)
require(dplyr)
require(tidyverse)
require(tidyr)
require(ggplot2)
require(reshape2)
require(data.table)
require(Seurat)
require(rio)
require(GenomicRanges)
require(ggpubr)
require(TFBSTools)
require(JASPAR2020)
require(motifmatchr)
require(readr)
require(GeneOverlap)

args = commandArgs(trailingOnly = TRUE)

for (arg in args) {
	split_arg <- strsplit(arg, "=")[[1]]
	var_name <- split_arg[1]
	var_value <- split_arg[2]
	if(grepl(',', var_value)){
		var_value = strsplit(var_value, ',') %>% unlist()
	}
	print(var_name)
	print(var_value)
	assign(var_name, var_value)
}



if(new_or_ms == 'MS'){
	motres = readRDS('IN_DATA/AllMotifEvo.RDS')
} else if(new_or_ms == 'NEW') {
	fls = list.files(path = outdir, pattern = 'motifEvo.RDS', full.names = T)
	motres = lapply(fls, function(x){read_rds(x)}) %>% do.call(rbind, .)
	write_rds(motres, paste0(outdir, '/AllMotifEvo.RDS'))
} else{
	stop('Error: new_or_ms should be either NEW or MS')
}


plt_pref = paste0(outdir, '/', plt_pref)
dir.create(outdir, showWarnings = F)

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

pdf(paste0(plt_pref, '_GainLoss_Motifs.pdf'), width = 7, height = 5)
print( ggbarplot(toplotM, x = 'variable2', y = 'value', fill = 'variable2',
		color = 'variable2',  position = position_dodge(0.9)) +
theme_classic() + xlab('') + ylab('Percentage') +
theme(text = element_text(size=20, face = 'bold')) +
facet_wrap(~node, nrow = 1) + NoLegend() +
rotate_x_text(45) )
dev.off()


####
## Count Gains and Losses per node
####

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

finaldf = rbindlist(dfL, fill = T) %>% as.data.frame
finaldf$ratio = finaldf$Gain / finaldf$Loss

saveRDS(finaldf, paste0(outdir, '/MotifEvoCounts.RDS'))

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

pdf(paste0(plt_pref, '_GainLoss_Cor_AcrossMotifs_Adult.pdf'), width = 8, height = 6)
print( ggplot(toplot,aes(x=Var1,y=Var2,fill=value))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0.2, low = 'blue', high = 'red', limits=c(0,1))+
  theme_classic()+
  geom_text(aes(label = label), vjust = 0.8, colour = "black", size = 6) +
  theme(text = element_text(size=20, face = 'bold')) )
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

saveRDS(dfNode, paste0(outdir, '/Motif_Test_AllStats.RDS'))

enrDF = dfNode[dfNode$FDR < 0.01 & dfNode$FDR_Bcg < 0.01 & dfNode$oddsRatio > 1 & dfNode$oddsRatioBcg > 1,]
depDF = dfNode[dfNode$FDR < 0.01 & dfNode$FDR_Bcg < 0.01 & dfNode$oddsRatio < 1 & dfNode$oddsRatioBcg < 1,]

saveRDS(enrDF, paste0(outdir, '/Motif_Test_Enriched.RDS'))
saveRDS(depDF, paste0(outdir, '/Motif_Test_Depleted.RDS'))

####
## Filter by expression / accessibility
###

hres = enrDF[enrDF$node == 'Human', 'tf']
hcres = enrDF[enrDF$node == 'HC', 'tf']
aperes = enrDF[enrDF$node == 'Ape', 'tf']

# RNA expressed genes
rnaSeurat = read_rds(snrnaseq)

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
expGns2 = unlist(expGnsL) %>% unique
saveRDS(expGns2, paste0(outdir, '/expressed_genes.RDS'))

# ATAC open chromatin genes
atacSeurat = read_rds(snATACseq)

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
accGns2 = unlist(accGnsL) %>% unique
saveRDS(accGns2, paste0(outdir, '/accessible_genes.RDS'))

# Filter by expression or accessibility
accGns2 = readRDS(paste0(outdir, '/accessible_genes.RDS'))
expGns2 = readRDS(paste0(outdir, '/expressed_genes.RDS'))

hres2 = hres[hres %in% expGns2 | hres %in% accGns2] %>% sort
hcres2 = hcres[hcres %in% expGns2 | hcres %in% accGns2] %>% sort
aperes2 = aperes[aperes %in% expGns2 | aperes %in% accGns2] %>% sort

# Keep only the expressed / accessible TFs
enrDF2 = enrDF[enrDF$tf %in% Reduce(union, list(hres2, hcres2, aperes2)),]
depDF2 = depDF[depDF$tf %in% Reduce(union, list(hres2, hcres2, aperes2)),]
saveRDS(enrDF2, paste0(outdir, '/Motif_Test_Enriched_ExpressedAccessible.RDS'))
saveRDS(depDF2, paste0(outdir, '/Motif_Test_Depleted_ExpressedAccessible.RDS'))


# Plot number of TFs
tmp = enrDF2
tmp$node = factor(tmp$node, levels = nodes)
toplot = table(tmp$node) %>% as.data.frame

pdf(paste0(plt_pref, '_number_of_significantly_expanded_motifs.pdf'))
print( ggbarplot(toplot, x = 'Var1', y = 'Freq', color = 'gray', fill = 'gray') +
ylab('Number of significantly expanded TFs') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend() )
dev.off()

tmp = depDF2
tmp$node = factor(tmp$node, levels = nodes)
toplot = table(tmp$node) %>% as.data.frame

pdf(paste0(plt_pref, '_number_of_significantly_depleted_motifs.pdf'))
print( ggbarplot(toplot, x = 'Var1', y = 'Freq', color = 'gray', fill = 'gray') +
ylab('Number of significantly expanded TFs') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend() )
dev.off()



# Plot enrichment stats of all TFs

# HUMAN #
enrDF2_human = enrDF2[enrDF2$node == 'Human', 'tf']
toplot = dfNode[dfNode$tf %in% enrDF2_human,]
toplot$log10FDR = -log10(toplot$FDR)
toplot$is_sign = ifelse(toplot$FDR_Bcg < 0.05 & toplot$FDR < 0.05, '**',
			ifelse(!(toplot$FDR_Bcg < 0.05) & toplot$FDR < 0.05, '*', ''))

toplot$node = factor(toplot$node, levels = rev(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')))

pdf(paste0(plt_pref, '_Human_AllDivTFs_EnrichStats.pdf'), width = 20, height = 6)
print( ggscatter(toplot, x = 'tf', y = 'node', color = 'oddsRatio', size = 'log10FDR') +
  scale_size_continuous(range = c(4,15)) +
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  geom_text(aes(label = is_sign), vjust = 0.8, colour = "darkgreen", size = 12 ) +
  rotate_x_text(45) )
dev.off()

# APE #
enrDF2_human = enrDF2[enrDF2$node == 'Ape', 'tf']
toplot = dfNode[dfNode$tf %in% enrDF2_human,]
toplot$log10FDR = -log10(toplot$FDR)

toplot$is_sign = ifelse(toplot$FDR_Bcg < 0.05 & toplot$FDR < 0.05, '**',
			ifelse(!(toplot$FDR_Bcg < 0.05) & toplot$FDR < 0.05, '*', ''))

toplot$node = factor(toplot$node, levels = rev(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')))

pdf(paste0(plt_pref, '_ape_AllDivTFs_EnrichStats.pdf'), width = 20, height = 6)
print( ggscatter(toplot, x = 'tf', y = 'node', color = 'oddsRatio', size = 'log10FDR') +
  scale_size_continuous(range = c(8,15)) +
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  geom_text(aes(label = is_sign), vjust = 0.8, colour = "darkgreen", size = 12 ) +
  rotate_x_text(45) )
dev.off()


