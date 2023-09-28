rm(list = ls())
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
library(ggrepel)

# Set variables
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons')
type = 'TOP20K_EXPAND_50kb'
dev = 'ADULT'

# Combine regression p-values
final_list = list()
for(j in 1:length(nodes)){

	subdir = paste0(dev, '_', nodes[j])
	resf = list.files(path = paste0('~/workdir/pr5/05_LDSC/', type, '/RESULTS/', subdir), pattern = 'cell', full = T)

	resl = list()
	for(i in 1:length(resf)){
		res = read.table(resf[i], header = T)
		fulln = gsub(".*/", '', resf[i]) %>% gsub('\\..*', '', .)
		res$Categ = gsub('.*_', '', fulln)
		res$dev = gsub('_.*', '', fulln)
		res$node = nodes[j]
		resl[[i]] = res
	}

	final_list[[j]] = do.call(rbind, resl)
	final_list[[j]]$node = nodes[j]
}

finaldf = do.call(rbind, final_list)
finaldf$Categ = gsub('AUT', 'ASD', finaldf$Categ)
finaldf = finaldf[!(finaldf$Categ %in% c('HEIGHT', 'OCD', 'ANXIETY')),]

# Reshape for plotting
finaldf2 = split(finaldf, finaldf$Categ) %>% do.call(rbind, .)
finaldf2$FDR = split(finaldf, finaldf$Categ) %>% lapply(., function(x){p.adjust(x$Coefficient_P_value, method = 'BH')}) %>% unlist
finaldf2$log10FDR = -log10(finaldf2$FDR)
finaldf2$is_sign = ifelse(finaldf2$FDR < 0.01, '**',
			ifelse(finaldf2$FDR < 0.1, '*', ''))


finaldf2$CoefScaled = finaldf2$Coefficient * 1e+9
finaldf2$Categ = factor(finaldf2$Categ, levels = c('AD', 'ADHD', 'ASD', 'BP', 'SCZ', 'MDD', 'INT', 'COG', 'OST', 'CAD'))
finaldf2$node = factor(finaldf2$node, levels = rev(nodes))

toplot = finaldf2
toplot$CoefScaled = ifelse(toplot$CoefScaled > 20, 20,
			ifelse(toplot$CoefScaled < -20, -20, toplot$CoefScaled))

pdf(paste0('LDSC_', dev, '_', type, '.pdf'), width = 10, height = 6)
ggscatter(toplot, x = 'Categ', y = 'node', color = 'CoefScaled', size = 'log10FDR') +
rotate_x_text(90) + theme(text=element_text(size=20)) + xlab('') + ylab('') +
geom_text(aes(label = is_sign), vjust = 0.8, colour = "darkgreen", size = 12 ) +
scale_size_continuous(range = c(6,20), limits = c(0,12)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20), legend.pos = 'right') +
scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0, limits = c(-20,20))
dev.off()


saveRDS(finaldf2, paste0(type, '/RESULTS/', dev, '_All_Stats.RDS'))





