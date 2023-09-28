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
source('Functions.R')

####
## MOTIF EVOLUTION AND CHROMATIN ACCESSIBILITY ENRICHMENT
####

# Load DAR results
darsAll = readRDS('DATASETS/Adult_SpeciesSpecific_DARs_snATACseq.RDS')

# Species specific dars
darsH = darsAll[darsAll$Evolution == 'Human_Specific', 'Gene'] %>% unique
darsC = darsAll[darsAll$Evolution == 'Chimp_Specific', 'Gene'] %>% unique
darsM = darsAll[darsAll$Evolution == 'MvsHC', 'Gene'] %>% unique
tmp = unique(darsAll$Gene)

# CREs with no accessibility difference across species
darsCons = setdiff(tmp, unique(c(darsH, darsC, darsM))) %>% unique
pbmarkL = list(darsH, darsCons)
names(pbmarkL) = c('HumanAcc', 'ConservedAcc')

# TF enrichment results
motevo = readRDS('Adult_RawMotifEnrich_ExpressedAccessible.RDS')
human_tfs = motevo[motevo$node == 'Human', 'tf']
hc_tfs = motevo[motevo$node == 'HC', 'tf']
ape_tfs = motevo[motevo$node == 'Ape', 'tf']

# Load motif evolution results
motres = readRDS('Adult_AllMotifEvo.RDS')
motres$type = factor(motres$type)

# Create the TF List -- HUMAN-HOMININ TFs
tf = human_tfs
human_tf = motres[motres$tf %in% tf & motres$node == 'Human', 'CRE'] %>% unique
hc_tf = motres[motres$tf %in% tf & motres$node == 'HC', 'CRE'] %>% unique
hcgo_tf = motres[motres$tf %in% tf & motres$node == 'HCGo', 'CRE'] %>% unique
gape_tf = motres[motres$tf %in% tf & motres$node == 'Great_Ape', 'CRE'] %>% unique
ape_tf = motres[motres$tf %in% tf & motres$node == 'Ape', 'CRE'] %>% unique

human_tf2 = setdiff(human_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, ape_tf)))
hc_tf2 = setdiff(hc_tf, Reduce(union, list(human_tf, hcgo_tf, gape_tf, ape_tf)))
hcgo_tf2 = setdiff(hcgo_tf, Reduce(union, list(hc_tf, human_tf, gape_tf, ape_tf)))
gape_tf2 = setdiff(gape_tf, Reduce(union, list(hc_tf, hcgo_tf, human_tf, ape_tf)))
ape_tf2 = setdiff(ape_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, human_tf)))

tfL = list(human_tf2, hc_tf2, hcgo_tf2, gape_tf2, ape_tf2)
names(tfL) = paste0('HumanExpandedTFs_', c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))
human_tfL = tfL

# Create the TF List -- HOMININ TFs
tf = hc_tfs
human_tf = motres[motres$tf %in% tf & motres$node == 'Human', 'CRE'] %>% unique
hc_tf = motres[motres$tf %in% tf & motres$node == 'HC', 'CRE'] %>% unique
hcgo_tf = motres[motres$tf %in% tf & motres$node == 'HCGo', 'CRE'] %>% unique
gape_tf = motres[motres$tf %in% tf & motres$node == 'Great_Ape', 'CRE'] %>% unique
ape_tf = motres[motres$tf %in% tf & motres$node == 'Ape', 'CRE'] %>% unique

human_tf2 = setdiff(human_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, ape_tf)))
hc_tf2 = setdiff(hc_tf, Reduce(union, list(human_tf, hcgo_tf, gape_tf, ape_tf)))
hcgo_tf2 = setdiff(hcgo_tf, Reduce(union, list(hc_tf, human_tf, gape_tf, ape_tf)))
gape_tf2 = setdiff(gape_tf, Reduce(union, list(hc_tf, hcgo_tf, human_tf, ape_tf)))
ape_tf2 = setdiff(ape_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, human_tf)))

tfL = list(human_tf2, hc_tf2, hcgo_tf2, gape_tf2, ape_tf2)
names(tfL) = paste0('HomininExpandedTFs_', c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))
hc_tfL = tfL

# Create the TF List -- Ape TFs
tf = ape_tfs
human_tf = motres[motres$tf %in% tf & motres$node == 'Human', 'CRE'] %>% unique
hc_tf = motres[motres$tf %in% tf & motres$node == 'HC', 'CRE'] %>% unique
hcgo_tf = motres[motres$tf %in% tf & motres$node == 'HCGo', 'CRE'] %>% unique
gape_tf = motres[motres$tf %in% tf & motres$node == 'Great_Ape', 'CRE'] %>% unique
ape_tf = motres[motres$tf %in% tf & motres$node == 'Ape', 'CRE'] %>% unique

human_tf2 = setdiff(human_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, ape_tf)))
hc_tf2 = setdiff(hc_tf, Reduce(union, list(human_tf, hcgo_tf, gape_tf, ape_tf)))
hcgo_tf2 = setdiff(hcgo_tf, Reduce(union, list(hc_tf, human_tf, gape_tf, ape_tf)))
gape_tf2 = setdiff(gape_tf, Reduce(union, list(hc_tf, hcgo_tf, human_tf, ape_tf)))
ape_tf2 = setdiff(ape_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, human_tf)))

tfL = list(human_tf2, hc_tf2, hcgo_tf2, gape_tf2, ape_tf2)
names(tfL) = paste0('ApeExpandedTFs_', c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))
ape_tfL = tfL

# Create the TF List -- Other TFs
tf = setdiff(unique(motres$tf), Reduce(union, list(human_tfs, hc_tfs, ape_tfs)))
human_tf = motres[motres$tf %in% tf & motres$node == 'Human', 'CRE'] %>% unique
hc_tf = motres[motres$tf %in% tf & motres$node == 'HC', 'CRE'] %>% unique
hcgo_tf = motres[motres$tf %in% tf & motres$node == 'HCGo', 'CRE'] %>% unique
gape_tf = motres[motres$tf %in% tf & motres$node == 'Great_Ape', 'CRE'] %>% unique
ape_tf = motres[motres$tf %in% tf & motres$node == 'Ape', 'CRE'] %>% unique

human_tf2 = setdiff(human_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, ape_tf)))
hc_tf2 = setdiff(hc_tf, Reduce(union, list(human_tf, hcgo_tf, gape_tf, ape_tf)))
hcgo_tf2 = setdiff(hcgo_tf, Reduce(union, list(hc_tf, human_tf, gape_tf, ape_tf)))
gape_tf2 = setdiff(gape_tf, Reduce(union, list(hc_tf, hcgo_tf, human_tf, ape_tf)))
ape_tf2 = setdiff(ape_tf, Reduce(union, list(hc_tf, hcgo_tf, gape_tf, human_tf)))

tfL = list(human_tf2, hc_tf2, hcgo_tf2, gape_tf2, ape_tf2)
names(tfL) = paste0('OtherTFs_', c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))
cons_tfL = tfL

# Run enrichments
evoL = c(human_tfL, hc_tfL, ape_tfL, cons_tfL)
bcg = length(motres$CRE %>% unique)
toplot = geneOvEnrDep(gnL1 = evoL, gnL2 = pbmarkL, bcg, plot = F, hg = 6, wd = 20)
toplot$var1 = gsub('Great_Ape', 'GreatApe', toplot$var1)
toplot$lineage = gsub('.*_', '', toplot$var1)
toplot$tfGroup = gsub('_.*', '', toplot$var1)
toplot$tfGroup = factor(toplot$tfGroup, levels = c('HumanExpandedTFs', 'HomininExpandedTFs', 'ApeExpandedTFs', 'OtherTFs'))

pdf('Motif_Accessibility_Association_GAINLOSS_IN_ALL_LINEAGES.pdf', height = 6, width = 20)
ggscatter(toplot, x = 'lineage', y = 'var2', color = 'OddsRatio', size = 'log10FDR') +
		geom_label(data = toplot, aes(label = sign_label), color="black",
			label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
		labs(x='', y='') +
		scale_size_continuous(range = c(10,25)) +
		scale_color_gradient2(midpoint = 1, low = 'lightblue', high = 'red') +
		theme_classic() +
		facet_wrap(~tfGroup) +
		theme(text = element_text(size=20, face = 'bold')) +
		geom_text(aes(label = sign_label_2), vjust = 1.3, colour = "darkgreen", fontface = 'bold', size = 20 ) +
		rotate_x_text(45)
dev.off()

# Only keep human changes and rerun
evoL = c(human_tfL['HumanExpandedTFs_Human'], ape_tfL['ApeExpandedTFs_Human'], cons_tfL['OtherTFs_Human'])
names(evoL) = c('HumanExpandedTFs', 'ApeExpandedTFs', 
'OtherTFs')

bcg = length(motres$CRE %>% unique)
toplot = geneOvEnrDep(gnL1 = evoL, gnL2 = pbmarkL, bcg, plot = F, hg = 6, wd = 20)
toplot$var1 = gsub('Great_Ape', 'GreatApe', toplot$var1)

toplot$var1 = factor(toplot$var1, levels = c('HumanExpandedTFs', 'ApeExpandedTFs', 'OtherTFs'))
toplot$OddsRatio = ifelse(toplot$OddsRatio < 0.8, 0.8, toplot$OddsRatio)

pdf('Motif_Accessibility_Association_GAINLOSS_IN_HUMAN.pdf', height = 8, width = 10)
ggscatter(toplot, x = 'var1', y = 'var2', color = 'OddsRatio', size = 'log10FDR') +
		geom_label(data = toplot, aes(label = sign_label), color="black",
			label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
		labs(x='', y='') +
		scale_size_continuous(range = c(10,25)) +
		scale_color_gradient2(midpoint = 1.1, low = 'lightblue', high = 'red', limits = c(0.8, 2)) +
		theme_classic() +
		theme(text = element_text(size=20, face = 'bold')) +
		geom_text(aes(label = sign_label_2), vjust = 1.3, colour = "darkgreen", fontface = 'bold', size = 20 ) +
		rotate_x_text(45)
dev.off()


