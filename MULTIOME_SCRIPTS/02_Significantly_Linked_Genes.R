rm(list = ls())
library(Signac)
library(Seurat)
library(Matrix)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(readr)
library(tgutil)
library(tidyverse)
library(ggpubr)
source("utility_functions.R")

# Read data
pks_to_gns = read_rds('MULTIOME_FETAL/PKS_TO_GNS/Trevino_LinkPeaksToGenes_FINAL.RDS')
pks_to_gns$FDR = p.adjust(pks_to_gns$pval, method = 'fdr')

# FDR cutoff of 0.05 does not filter additional links in this dataset
pks_to_gns$peak2 = sub('-', ':', pks_to_gns$peak)
pks_to_gns$gene = factor(pks_to_gns$gene)
gns_bcg = pks_to_gns$gene %>% unique

hvarsAll = read_rds('02_GroupCREs/FETAL/FetalCREs_SignTested_SD_1_FC_ConsAdded.RDS')

# All CREs
cre_all = hvarsAll$CRE
gns_cre_all = pks_to_gns[pks_to_gns$peak2 %in% cre_all, 'gene'] %>% table

####
## Human divergent CREs
####

cre_test = hvarsAll[hvarsAll$HumanSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$HumanSign == 'Sign')
randL = list()
for(i in 1:1000){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', i))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	randL[[i]] = gns_cre_rand
	if(i%%100 == 0){print(i)}
}

rand_df = do.call(cbind, randL)
rand_df = rand_df[names(gns_cre_test), ]

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
resL = list()
for(i in 1:length(gnstotest)){

	#if(gns_cre_test[gnstotest[i]] < 3){next}
	gn = gnstotest[i]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[gn,] %>% as.numeric %>% median

	pval = sum(rand_df[gn,] > obs) / 1000
	resL[[i]] = data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
	if(i%%100 == 0){print(i)}
}

resdf = do.call(rbind, resL)
saveRDS(resdf, 'MULTIOME_FETAL/Linked_Gene_Stats_Human.RDS')


####
## HC divergent CREs
####

cre_test = hvarsAll[hvarsAll$HCSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$HCSign == 'Sign')
randL = list()
for(i in 1:1000){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', i))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	randL[[i]] = gns_cre_rand
	if(i%%100 == 0){print(i)}
}

rand_df = do.call(cbind, randL)
rand_df = rand_df[names(gns_cre_test), ]

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
resL = list()
for(i in 1:length(gnstotest)){

	#if(gns_cre_test[gnstotest[i]] < 3){next}
	gn = gnstotest[i]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[gn,] %>% as.numeric %>% median

	pval = sum(rand_df[gn,] > obs) / 1000
	resL[[i]] = data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
	if(i%%100 == 0){print(i)}
}

resdf = do.call(rbind, resL)
saveRDS(resdf, 'MULTIOME_FETAL/Linked_Gene_Stats_HC.RDS')

####
## HCGo divergent CREs
####

cre_test = hvarsAll[hvarsAll$HCGoSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$HCGoSign == 'Sign')
randL = list()
for(i in 1:1000){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', i))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	randL[[i]] = gns_cre_rand
	if(i%%100 == 0){print(i)}
}

rand_df = do.call(cbind, randL)
rand_df = rand_df[names(gns_cre_test), ]

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
resL = list()
for(i in 1:length(gnstotest)){

	#if(gns_cre_test[gnstotest[i]] < 3){next}
	gn = gnstotest[i]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[gn,] %>% as.numeric %>% median

	pval = sum(rand_df[gn,] > obs) / 1000
	resL[[i]] = data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
	if(i%%100 == 0){print(i)}
}

resdf = do.call(rbind, resL)
saveRDS(resdf, 'MULTIOME_FETAL/Linked_Gene_Stats_HCGo.RDS')

####
## Great_Ape divergent CREs
####

cre_test = hvarsAll[hvarsAll$Great_ApeSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$Great_ApeSign == 'Sign')
randL = list()
for(i in 1:1000){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', i))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	randL[[i]] = gns_cre_rand
	if(i%%100 == 0){print(i)}
}

rand_df = do.call(cbind, randL)
rand_df = rand_df[names(gns_cre_test), ]

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
resL = list()
for(i in 1:length(gnstotest)){

	#if(gns_cre_test[gnstotest[i]] < 3){next}
	gn = gnstotest[i]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[gn,] %>% as.numeric %>% median

	pval = sum(rand_df[gn,] > obs) / 1000
	resL[[i]] = data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
	if(i%%100 == 0){print(i)}
}

resdf = do.call(rbind, resL)
saveRDS(resdf, 'MULTIOME_FETAL/Linked_Gene_Stats_Great_Ape.RDS')

####
## Ape divergent CREs
####

cre_test = hvarsAll[hvarsAll$ApeSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$ApeSign == 'Sign')
randL = list()
for(i in 1:1000){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', i))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	randL[[i]] = gns_cre_rand
	if(i%%100 == 0){print(i)}
}

rand_df = do.call(cbind, randL)
rand_df = rand_df[names(gns_cre_test), ]

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
resL = list()
for(i in 1:length(gnstotest)){

	#if(gns_cre_test[gnstotest[i]] < 3){next}
	gn = gnstotest[i]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[gn,] %>% as.numeric %>% median

	pval = sum(rand_df[gn,] > obs) / 1000
	resL[[i]] = data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
	if(i%%100 == 0){print(i)}
}

resdf = do.call(rbind, resL)
saveRDS(resdf, 'MULTIOME_FETAL/Linked_Gene_Stats_Ape.RDS')


####
## Read all results
####

human_resdf = read_rds('MULTIOME_FETAL/Linked_Gene_Stats_Human.RDS') %>% add_column(Divergence = 'Human')
hc_resdf = read_rds('MULTIOME_FETAL/Linked_Gene_Stats_HC.RDS') %>% add_column(Divergence = 'HC')
hcgo_resdf = read_rds('MULTIOME_FETAL/Linked_Gene_Stats_HCGo.RDS') %>% add_column(Divergence = 'HCGo')
greatApe_resdf = read_rds('MULTIOME_FETAL/Linked_Gene_Stats_Great_Ape.RDS') %>% add_column(Divergence = 'Great_Ape')
ape_resdf = read_rds('MULTIOME_FETAL/Linked_Gene_Stats_Ape.RDS') %>% add_column(Divergence = 'Ape')

# Multiple testing correction
human_resdf$FDR = p.adjust(human_resdf$pval, method = 'fdr')
hc_resdf$FDR = p.adjust(hc_resdf$pval, method = 'fdr')
hcgo_resdf$FDR = p.adjust(hcgo_resdf$pval, method = 'fdr')
greatApe_resdf$FDR = p.adjust(greatApe_resdf$pval, method = 'fdr')
ape_resdf$FDR = p.adjust(ape_resdf$pval, method = 'fdr')


# Divergent genes. Only the ones with at least N peaks linked to it among all peaks.
bcg_pks_thresh = 5
obs_pks_thresh = 2

human_div = human_resdf[human_resdf$pval < 0.05 &
			human_resdf$all_peaks >= bcg_pks_thresh & 
			human_resdf$obs_peaks >= obs_pks_thresh, ]

hc_div = hc_resdf[hc_resdf$pval < 0.05 &
		hc_resdf$all_peaks >= bcg_pks_thresh &
		hc_resdf$obs_peaks >= obs_pks_thresh, ]

hcgo_div = hcgo_resdf[hcgo_resdf$pval < 0.05 &
			hcgo_resdf$all_peaks >= bcg_pks_thresh &
			hcgo_resdf$obs_peaks >= obs_pks_thresh, ]

gape_div = greatApe_resdf[greatApe_resdf$pval < 0.05 &
			greatApe_resdf$all_peaks >= bcg_pks_thresh &
			greatApe_resdf$obs_peaks >= obs_pks_thresh, ]

ape_div = ape_resdf[ape_resdf$pval< 0.05 &
			ape_resdf$all_peaks >= bcg_pks_thresh &
			ape_resdf$obs_peaks >= obs_pks_thresh, ]


all_res = do.call(rbind, list(human_resdf, hc_resdf, hcgo_resdf, greatApe_resdf, ape_resdf))
all_div = do.call(rbind, list(human_div, hc_div, hcgo_div, gape_div, ape_div))

write_rds(all_res, 'MULTIOME_FETAL/Linked_Gene_Final_Results_All.RDS')
write_rds(all_div, 'MULTIOME_FETAL/Linked_Gene_Final_Results_Sign_Nominal.RDS')


# Keep only significant ones after FDR correction.
all_div = all_div[all_div$FDR < 0.05, ]
write_rds(all_div, 'MULTIOME_FETAL/Linked_Gene_Final_Results_Sign_FDR.RDS')

# Plot number of genes per category
all_div = read_rds('MULTIOME_FETAL/Linked_Gene_Final_Results_Sign_FDR.RDS')

groups = unique(all_div[['Divergence']])
vals = table(all_div[['Divergence']])[groups] %>% as.numeric
toplot = data.frame(vars = groups, vals = vals)
colors = c('blue', 'lightblue', 'purple', 'pink', 'red')

pdf('Fetal_number_of_divergent_genes.pdf')
ggbarplot(toplot, x = 'vars', y = 'vals', color = 'vars', fill = 'vars', palette = colors) +
	ylab('Number of significantly linked genes\n(FDR<0.05)') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	NoLegend()
dev.off()


# All genes tested for DEGs
alldeg_res = readRDS('PSEUDOBULK_DEGs_ALL.RDS')
all_gns_degres = alldeg_res$Gene %>% unique

# All genes tested for peak-link
seq_evo_bcg = read_rds('pr_gns.RDS') %>% unlist %>% as.character

####
## HS-DEG ENRICHMENT
####

# Species-specific differentially expressed genes
hsdegs = alldeg_res[alldeg_res$Evolution == 'Human_Specific', 'Gene'] %>% unique
#mvshcdegs = alldeg_res[alldeg_res$Evolution == 'MvsHC', 'Gene'] %>% unique

human_div = all_div[all_div$Divergence == 'Human', 'gene']
hc_div = all_div[all_div$Divergence == 'HC', 'gene']
hcgo_div = all_div[all_div$Divergence == 'HCGo', 'gene']
gape_div = all_div[all_div$Divergence == 'Great_Ape', 'gene']
ape_div = all_div[all_div$Divergence == 'Ape', 'gene']


# Observed
human_obs_ov = sum(human_div %in% hsdegs) / length(human_div)
hc_obs_ov = sum(hc_div %in% hsdegs) / length(hc_div)
hcgo_obs_ov = sum(hcgo_div %in% hsdegs) / length(hcgo_div)
gape_obs_ov = sum(gape_div %in% hsdegs) / length(gape_div)
ape_obs_ov = sum(ape_div %in% hsdegs) / length(ape_div)

# Calculate p value compared to background
human_rand_ovs = sapply(1:1000, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% human_div) / length(human_div)})
human_pval = sum(human_rand_ovs > human_obs_ov) / 1000

hc_rand_ovs = sapply(1:1000, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% hc_div) / length(hc_div)})
hc_pval = sum(hc_rand_ovs > hc_obs_ov) / 1000

hcgo_rand_ovs = sapply(1:1000, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% hcgo_div) / length(hcgo_div)})
hcgo_pval = sum(hcgo_rand_ovs > hcgo_obs_ov) / 1000

gape_rand_ovs = sapply(1:1000, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% gape_div) / length(gape_div)})
gape_pval = sum(gape_rand_ovs > gape_obs_ov) / 1000

ape_rand_ovs = sapply(1:1000, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% ape_div) / length(ape_div)})
ape_pval = sum(ape_rand_ovs > ape_obs_ov) / 1000


rands = rbind(human_rand_ovs, hc_rand_ovs, hcgo_rand_ovs, gape_rand_ovs, ape_rand_ovs) %>% as.data.frame
colnames(rands) = paste0('Randomized_', 1:ncol(rands))
rands$Observed = c(human_obs_ov, hc_obs_ov, hcgo_obs_ov, gape_obs_ov, gape_obs_ov)
rands$ids = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')

toplot = reshape2::melt(rands)
toplot$ids = factor(toplot$ids, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))

pdf('HSDEG_SequenceEvo_Enrichment_FETAL.pdf')
ggplot(toplot, aes(x = ids, y = value)) +
geom_boxplot(outlier.shape = NA, color = 'blue') +
theme_classic() +
theme(text = element_text(size=20)) +
xlab('') + ylab('Overlap_of_peakLinkedGenes\nand_HS-DEGs') +
rotate_x_text(90) +
geom_point(data = toplot[toplot$variable == 'Observed',], color = 'red', size = 3) +
theme(strip.text.x = element_text(angle = 90))
dev.off()

# Export p-value table
toexport = toplot[toplot$variable == 'Observed',]
toexport$pval = c(human_pval, hc_pval, hcgo_pval, gape_pval, ape_pval)

library(gridExtra)
pdf('HSDEG_SequenceEvo_Enrichment_TABLE_FETAL.pdf')
grid.table(toexport)
dev.off()

write_rds(toplot, 'MULTIOME_FETAL/HSDEG_SequenceEvo_Enrichment.RDS')



####
## BARPLOT RAGs
####

# Add background number of CRE links for all genes
all_res = read_rds('MULTIOME_FETAL/Linked_Gene_Final_Results_All.RDS')
all_div = read_rds('MULTIOME_FETAL/Linked_Gene_Final_Results_Sign_FDR.RDS')

all_res$id = paste0(all_res$gene, '_', all_res$Divergence)
all_div$id = paste0(all_div$gene, '_', all_div$Divergence)

toplot = all_res
toplot$is_sign = ifelse(toplot$id %in% all_div$id, 'Sign', 'NS')
#toplot$is_sign = ifelse(toplot$FDR < 0.05, 'Sign', 'NS')
toplot$is_sign = factor(toplot$is_sign, levels = c('NS', 'Sign'))

tmp = toplot %>% group_by(gene) %>% summarize(obs_peaks = round(mean(bcg_peaks_mean))) %>% as.data.frame
tmp$is_sign = 'NS'
tmp$Divergence = 'Background'

toplot = toplot[, c('gene', 'obs_peaks', 'is_sign', 'Divergence')]
toplot = rbind(toplot, tmp)
toplot$id = paste0(toplot$gene, '_', toplot$Divergence)
toplot = toplot[!(duplicated(toplot$id)),]

write_rds(toplot, 'MULTIOME_FETAL/Linked_Gene_Final_Results_All_FORPLOTTING.RDS')


# Plot selected genes
toplot = read_rds('MULTIOME_FETAL/Linked_Gene_Final_Results_All_FORPLOTTING.RDS')
gns = c('TFEB')
toplot_sub = toplot[toplot$gene %in% gns,]

pdf('TFEB_RAG_Fetal.pdf', height = 5, width = 5)
ggbarplot(toplot_sub, x = 'Divergence', y = 'obs_peaks', color = 'is_sign', palette = c('lightblue', 'red'), fill = 'lightblue', size = 2) +
	ylab('Number of linked\n and divergent CREs') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	facet_wrap(~gene, nrow = 2) +
	NoLegend()
dev.off()



####
## ADULT-FETAL RAG OVERLAP
####

adult_all_div = read_rds('MULTIOME_ADULT/Linked_Gene_Final_Results_Sign_FDR.RDS')
fetal_all_div = read_rds('MULTIOME_FETAL/Linked_Gene_Final_Results_Sign_FDR.RDS')

nodes = unique(adult_all_div$Divergence)

for(i in 1:length(nodes)){

	q1 = adult_all_div[adult_all_div$Divergence == nodes[i], 'gene']
	q2 = fetal_all_div[fetal_all_div$Divergence == nodes[i], 'gene']

	gnL = list(q1, q2)
	names(gnL) = c(paste0("Adult\n", nodes[i], "\nRAGs"), paste0("Fetal\n", nodes[i], "\nRAGs"))
	fn = paste0('Fetal_Adult_RAG_Overlap_', nodes[i])
	plotVenn(gnL, fn = fn, margins = 0.1)
}



####
## RAG OVERLAP ACROSS LINEAGES
####

# Upset plot
library(UpSetR)
library(ComplexUpset)

# ADULT #

# List of gene sets
gnL = split(adult_all_div, adult_all_div$Divergence)
gnL2 = lapply(gnL, function(x){x$gene})
names(gnL2) = names(gnL)

# Boolean matrix of genes per set
allgns = Reduce(union, gnL2)
logmat = sapply(names(gnL2), function(x){allgns %in% gnL2[[x]]}) %>% as.data.frame
rownames(logmat) = allgns

logmat = logmat[,rev(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))]
sps = colnames(logmat)

# Display all non-zero intersections
pdf('ADULT_RAG_Overlap_UpsetPlot.pdf', width = 8)
ComplexUpset::upset(logmat, sps, name='sps', width_ratio=0.1, sort_sets = F,
set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90))),
themes=upset_modify_themes(
        list(
            'intersections_matrix'=theme(text=element_text(size=30)),
            'Intersection size'=theme(text=element_text(size=30)),
            'overall_sizes'=theme(text=element_text(size=30)),
            'default'=theme(text=element_text(size=30))
        )))
dev.off()



# FETAL #

# List of gene sets
gnL = split(fetal_all_div, fetal_all_div$Divergence)
gnL2 = lapply(gnL, function(x){x$gene})
names(gnL2) = names(gnL)

# Boolean matrix of genes per set
allgns = Reduce(union, gnL2)
logmat = sapply(names(gnL2), function(x){allgns %in% gnL2[[x]]}) %>% as.data.frame
rownames(logmat) = allgns

logmat = logmat[,rev(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))]
sps = colnames(logmat)

# Display all non-zero intersections
pdf('FETAL_RAG_Overlap_UpsetPlot.pdf', width = 12)
ComplexUpset::upset(logmat, sps, name='sps', width_ratio=0.1, sort_sets = F,
set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90))),
themes=upset_modify_themes(
        list(
            'intersections_matrix'=theme(text=element_text(size=30)),
            'Intersection size'=theme(text=element_text(size=30)),
            'overall_sizes'=theme(text=element_text(size=30)),
            'default'=theme(text=element_text(size=30))
        )))
dev.off()





