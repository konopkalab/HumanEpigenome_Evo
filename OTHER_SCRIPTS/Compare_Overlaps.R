rm(list = ls())
library('bedr')
library(ape)
library(dplyr)
library(parallel)
library(Biostrings)
library(ggpubr)
library(Seurat)
library(reshape2)
library(pheatmap)
source('Functions.R')

# Load datasets
chars = rio::import('HAR_CORTEX_CRE_FINAL_DF.xlsx')
chars = chars[!duplicated(chars$CRE),]
hars = rio::import('Doan_Combined_HARs_hg38.xlsx')
hvarsAll = readRDS('AdultCREs_SignTested_SD_1_FC_ConsAdded.RDS')
f_hvarsAll = readRDS('FetalCREs_SignTested_SD_1_FC_ConsAdded.RDS')

# Sort regions
chars = bedr.sort.region(chars, verbose = F)
hars = bedr.sort.region(hars, verbose = F)

# Overlaps -- HAR, ADULT
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')
nodes2 = paste0(nodes, 'Sign')
dfL = list()
for(i in 1:length(nodes)){

	# Prepare lineage divergent CREs
	cres = hvarsAll[ hvarsAll[[ nodes2[i] ]] == 'Sign', 'CRE']
	cres_bed = strToBed(cres, sep = ':|-')
	cres_bed = bedr.sort.region(cres_bed, verbose = F)

	# Overlap hars and divergent cres. Find the ratio of hars in divergent cres
	ov_hars = bedr(input = list(a = hars, b = cres_bed), method = "intersect", params = c("-loj"), verbose = F)
	rat = sum(ov_hars$V2 != '-1') / nrow(ov_hars)

	dfL[[i]] = data.frame(node = nodes[i], ratio = rat, type = 'HAR', dev = 'ADULT')
	print(i)
}

df1 = do.call(rbind, dfL)


# Overlaps -- CORTICAL HAR, ADULT
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')
nodes2 = paste0(nodes, 'Sign')
dfL = list()
for(i in 1:length(nodes)){

	# Prepare lineage divergent CREs
	cres = hvarsAll[ hvarsAll[[ nodes2[i] ]] == 'Sign', 'CRE']
	cres_bed = strToBed(cres, sep = ':|-')
	cres_bed = bedr.sort.region(cres_bed, verbose = F)

	# Overlap hars and divergent cres. Find the ratio of hars in divergent cres
	ov_hars = bedr(input = list(a = chars, b = cres_bed), method = "intersect", params = c("-loj"), verbose = F)
	rat = sum(ov_hars$V2 != '-1') / nrow(ov_hars)

	dfL[[i]] = data.frame(node = nodes[i], ratio = rat, type = 'CORTICAL_HAR', dev = 'ADULT')
	print(i)
}

df2 = do.call(rbind, dfL)


# Overlaps -- HAR, FETAL
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')
nodes2 = paste0(nodes, 'Sign')
dfL = list()
for(i in 1:length(nodes)){

	# Prepare lineage divergent CREs
	cres = f_hvarsAll[ f_hvarsAll[[ nodes2[i] ]] == 'Sign', 'CRE']
	cres_bed = strToBed(cres, sep = ':|-')
	cres_bed = bedr.sort.region(cres_bed, verbose = F)

	# Overlap hars and divergent cres. Find the ratio of hars in divergent cres
	ov_hars = bedr(input = list(a = hars, b = cres_bed), method = "intersect", params = c("-loj"), verbose = F)
	rat = sum(ov_hars$V2 != '-1') / nrow(ov_hars)

	dfL[[i]] = data.frame(node = nodes[i], ratio = rat, type = 'HAR')
	print(i)
}

df3 = do.call(rbind, dfL)


# Overlaps -- CORTICAL HAR, ADULT
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')
nodes2 = paste0(nodes, 'Sign')
dfL = list()
for(i in 1:length(nodes)){

	# Prepare lineage divergent CREs
	cres = f_hvarsAll[ f_hvarsAll[[ nodes2[i] ]] == 'Sign', 'CRE']
	cres_bed = strToBed(cres, sep = ':|-')
	cres_bed = bedr.sort.region(cres_bed, verbose = F)

	# Overlap hars and divergent cres. Find the ratio of hars in divergent cres
	ov_hars = bedr(input = list(a = chars, b = cres_bed), method = "intersect", params = c("-loj"), verbose = F)
	rat = sum(ov_hars$V2 != '-1') / nrow(ov_hars)

	dfL[[i]] = data.frame(node = nodes[i], ratio = rat, type = 'CORTICAL_HAR', dev = 'ADULT')
	print(i)
}

df4 = do.call(rbind, dfL)


# Plot overlaps
toplot = do.call(rbind, list(df1,df2))

pdf('HAR_Overlap_Comparison.pdf', width = 10, height = 5)
ggbarplot(toplot, x = 'node', y = 'ratio', fill = 'node', color = 'node') +
ylab('Ratio') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
facet_wrap(~type) +
NoLegend()
dev.off()







