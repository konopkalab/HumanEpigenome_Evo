rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(reshape2)
library(data.table)
library(ggpubr)
library(Seurat)

# BARPLOT TOTSNPs
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons')
dfL = list()
for(i in 1:length(nodes)){

	fls = list.files(path = paste0('TOP20K_EXPAND_50kb/ADULT_', nodes[i]), pattern = '.*M', full.names = T)
	totsnps = sapply(fls, function(x){read.table(x) %>% as.numeric}) %>% sum
	
	dfL[[i]] = data.frame(totsnps = totsnps, node = nodes[i])
}

toplot = do.call(rbind, dfL)

pdf('TOTSNPs_ADULT.pdf', width = 6, height = 5)
ggbarplot(toplot, y = 'totsnps', x = 'node', color = 'node', fill = 'node') +
rotate_x_text(90) + theme(text=element_text(size=20)) + xlab('') + ylab('Total number of SNPs') +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20), legend.pos = 'right') +
NoLegend()
dev.off()


