rm(list = ls())
library(dplyr)
library(parallel)
library(Biostrings)
library(ggpubr)
library(Seurat)
library(reshape2)
source('Functions.R')

# Load data
adultSubs = readRDS('02_GroupCREs/ADULT/AdultSubstitutions.RDS')
fetalSubs = readRDS('02_GroupCREs/FETAL/FetalSubstitutions.RDS')
harSubs = readRDS('02_GroupCREs/HAR/HARSubstitutions.RDS')
hvarsAll = readRDS('02_GroupCREs/ADULT/AdultCREs_SignTested_SD_1_FC_ConsAdded.RDS')
f_hvarsAll = readRDS('02_GroupCREs/FETAL/FetalCREs_SignTested_SD_1_FC_ConsAdded.RDS')

# Conversion rates -- HAR
q1 = harSubs[harSubs$nodes == 'Human',]
q1$change = paste0(q1[[ 'Human' ]], q1[[ 'HC' ]])

to_at = c('AC', 'AG', 'TC', 'TG')
to_gc = c('CA', 'GA', 'CT', 'GT')
at_at = c('AT', 'TA')
gc_gc = c('GC', 'CG')

r1 = sum(q1$change %in% to_at) / nrow(q1)
r2 = sum(q1$change %in% to_gc) / nrow(q1)
r3 = sum(q1$change %in% at_at) / nrow(q1)
r4 = sum(q1$change %in% gc_gc) / nrow(q1)

df1 = data.frame(vars = c('TO_AT', 'TO_GC', 'AT_AT', 'GC_GC'), vals = c(r1,r2,r3,r4))
df1$node = 'HAR'
df1$type = 'HAR'

# Conversion rates -- Adult
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Anthropoids')
nodes2 = paste0(nodes, 'Sign')
dfL = list()
for(i in 1:5){

	cres = hvarsAll[ hvarsAll[[ nodes2[i] ]] == 'Sign', 'CRE']
	q1 = adultSubs[adultSubs$nodes == nodes[i] & adultSubs$CRE %in% cres,]
	q1$change = paste0(q1[[ nodes[i] ]], q1[[ nodes[i+1] ]])

	to_at = c('AC', 'AG', 'TC', 'TG')
	to_gc = c('CA', 'GA', 'CT', 'GT')
	at_at = c('AT', 'TA')
	gc_gc = c('GC', 'CG')

	r1 = sum(q1$change %in% to_at) / nrow(q1)
	r2 = sum(q1$change %in% to_gc) / nrow(q1)
	r3 = sum(q1$change %in% at_at) / nrow(q1)
	r4 = sum(q1$change %in% gc_gc) / nrow(q1)

	dfL[[i]] = data.frame(vars = c('TO_AT', 'TO_GC', 'AT_AT', 'GC_GC'), vals = c(r1,r2,r3,r4))
	dfL[[i]]$node = nodes[i]
}

df2 = do.call(rbind, dfL)
df2$type = 'Adult'

# Conversion rates -- Fetal
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Anthropoids')
nodes2 = paste0(nodes, 'Sign')
dfL = list()
for(i in 1:5){toplot$

	cres = f_hvarsAll[ f_hvarsAll[[ nodes2[i] ]] == 'Sign', 'CRE']
	q1 = fetalSubs[fetalSubs$nodes == nodes[i] & fetalSubs$CRE %in% cres,]
	q1$change = paste0(q1[[ nodes[i] ]], q1[[ nodes[i+1] ]])

	to_at = c('AC', 'AG', 'TC', 'TG')
	to_gc = c('CA', 'GA', 'CT', 'GT')
	at_at = c('AT', 'TA')
	gc_gc = c('GC', 'CG')

	r1 = sum(q1$change %in% to_at) / nrow(q1)
	r2 = sum(q1$change %in% to_gc) / nrow(q1)
	r3 = sum(q1$change %in% at_at) / nrow(q1)
	r4 = sum(q1$change %in% gc_gc) / nrow(q1)

	dfL[[i]] = data.frame(vars = c('TO_AT', 'TO_GC', 'AT_AT', 'GC_GC'), vals = c(r1,r2,r3,r4))
	dfL[[i]]$node = nodes[i]
}

df3 = do.call(rbind, dfL)
df3$type = 'Fetal'

# Combine and plot
toplot = do.call(rbind, list(df1, df2, df3))
toplot$id = paste0(toplot$node, toplot$type)
toplot = toplot[toplot$vars == 'TO_GC',]

pdf('GC_Conversion_Rates_Comparison.pdf', width = 10)
ggbarplot(toplot, x = 'id', y = 'vals', fill = 'type', color = 'type') +
ylab('Ratio') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90)
dev.off()




