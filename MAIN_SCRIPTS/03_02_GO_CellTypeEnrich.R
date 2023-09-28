rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
library(Seurat)
library(dplyr)
library(tidyverse)
library(rio)
library(GenomicRanges)
library(data.table)
library(Signac)
library(pheatmap)
source('Functions.R')


# Load datasets
gores = readRDS('GO_allFinalDF.RDS')
cre_gene_goterm = readRDS('GO_Sign_Assoc_Gns_Pks.RDS')

# Terms to plot
sign_go = c("adenosine deaminase activity" , "lipoprotein particle binding", "methionine metabolic process")
nodes = c('Human', 'Human', 'Human')

####
## CELL TYPE ACCESSIBILITY HEATMAP
####

mat = readRDS('DATASETS/ATAC_pseudobulk_perBroadCellType_scaled.RDS')

aggdfL = list()
for(i in 1:length(sign_go)){

	div_cres = cre_gene_goterm[cre_gene_goterm$name == sign_go[i] & cre_gene_goterm$Divergence == nodes[i], c('seqnames', 'start', 'end')]
	div_cres = paste0(div_cres[,1], ':', div_cres[,2], '-', div_cres[,3])

	toplot = mat[div_cres,] %>% melt
	ctypesOrder = c("Excitatory", "Inhibitory", "OPC", "MOL", "Astrocyte", "Microglia")
	toplot$Var2 = factor(toplot$Var2, levels = ctypesOrder)

	aggdf = aggregate(value~Var2, data = toplot, FUN = mean)
	aggdf$name = sign_go[i]
	aggdfL[[i]] = aggdf
}

toplot = do.call(rbind, aggdfL)


# Subsets only
subsets = sign_go
toplot1 = toplot[toplot$name %in% subsets,]
toplot1$name = lapply(toplot1$name, function(x){wrapper(x, width = 20)}) %>% unlist
toplot1$name = factor(toplot1$name, levels = rev(unique(toplot1$name)))

pdf('ADULT_ALL_GO_CELLTYPE.pdf', width = 10, height = 6)
ggplot(toplot1, aes(x=Var2, y=name, fill = value))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic() + rotate_x_text(45) +
  theme(text = element_text(size=30, face = 'bold'), legend.pos = 'top')
dev.off()


####
## CELL TYPE EXPRESSION HEATMAP
####

mat = read_rds('DATASETS/RNA_pseudobulk_perBroadCellType_scaled.RDS')
all_div = read_rds('DATASETS/Linked_Gene_Final_Results_Sign_FDR.RDS')

pltL = list()
div_rag_gnsL = list()
for(i in 1:length(sign_go)){

	div_gns = cre_gene_goterm[cre_gene_goterm$name == sign_go[i] & cre_gene_goterm$Divergence == nodes[i], 'gene_symbol']
	div_gns = unique(div_gns)

	rags = all_div[all_div$Divergence == nodes[i], 'gene']
	div_rag_gnsL[[i]] = intersect(div_gns, rags)
	print(div_rag_gnsL[[i]])

	cmngns = intersect(rownames(mat), div_gns)
	if(length(cmngns) <= 1){next}

	toplot = mat[cmngns,]
	toplot = toplot[!(is.nan(toplot[,1])),]
	pltL[[i]] = pheatmap(t(toplot), fontsize = 18)
}


ggexport(pltL[[1]], filename = paste0(sign_go[1], '_GENES_ADULT_GO_CELLTYPE.pdf'), height = 4, width = 6)
ggexport(pltL[[2]], filename = paste0(sign_go[2], '_GENES_ADULT_GO_CELLTYPE.pdf'), height = 4, width = 10)
ggexport(pltL[[3]], filename = paste0(sign_go[3], '_GENES_ADULT_GO_CELLTYPE.pdf'), height = 4, width = 6)


####
## BARPLOT OF ALL GO ASSOCIATED GENES WITH # OF LINKED CRES
####

# Load peak-gene linkage dataset for plotting
toplot = read_rds('DATASETS/Linked_Gene_Final_Results_All_FORPLOTTING.RDS')

cre_linksL = list()
for(i in 1:length(sign_go)){

	# All genes per category
	allgo_gns = cre_gene_goterm[cre_gene_goterm$name == sign_go[i] & cre_gene_goterm$Divergence == 'Human', 'gene_symbol']
	allgo_gns = unique(allgo_gns)

	# Number of CREs linked to each gene
	cre_links = toplot[toplot$gene %in% allgo_gns & toplot$Divergence == nodes[i], ]

	# Remove the genes with no CRE link
	cre_links = cre_links[cre_links$obs_peaks != 0,]
	cre_links$go = sign_go[i]

	cre_linksL[[i]] = cre_links

}

toplot_links = do.call(rbind, cre_linksL)
levs = toplot_links[order(toplot_links$obs_peaks, decreasing = T), 'gene']
toplot_links$gene = factor(toplot_links$gene, levels = levs)

# Plot term 1
toplot_links_1 = toplot_links[toplot_links$go %in% sign_go[1],]
pdf(paste0(sign_go[1], '_Adult_GO_AllGenes.pdf'), height = 5, width = 5)
ggbarplot(toplot_links_1, x = 'gene', y = 'obs_peaks', color = 'is_sign', palette = c('lightblue', 'red'), fill = 'lightblue', size = 2) +
	ylab('Number of linked\n and divergent CREs') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20, face = 'italic'),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	NoLegend()
dev.off()

# Plot term 2
toplot_links_2 = toplot_links[toplot_links$go %in% sign_go[2],]
pdf(paste0(sign_go[2], '_Adult_GO_AllGenes.pdf'), height = 5, width = 10)
ggbarplot(toplot_links_2, x = 'gene', y = 'obs_peaks', color = 'is_sign', palette = c('lightblue', 'red'), fill = 'lightblue', size = 2) +
	ylab('Number of linked\n and divergent CREs') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20, face = 'italic'),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	NoLegend()
dev.off()


# Plot term 3
toplot_links_3 = toplot_links[toplot_links$go %in% sign_go[3],]
pdf(paste0(sign_go[3], '_Adult_GO_AllGenes.pdf'), height = 5, width = 5)
ggbarplot(toplot_links_3, x = 'gene', y = 'obs_peaks', color = 'is_sign', palette = c('lightblue', 'red'), fill = 'lightblue', size = 2) +
	ylab('Number of linked\n and divergent CREs') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20, face = 'italic'),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	NoLegend()
dev.off()

####
## BARPLOT OF GO ASSOCIATED RAGS
####

# Load peak-gene linkage dataset for plotting
toplot = read_rds('DATASETS/Linked_Gene_Final_Results_All_FORPLOTTING.RDS')

gnsL = list()
for(i in 1:length(sign_go)){

	# Genes per category
	allgo_gns = cre_gene_goterm[cre_gene_goterm$name == sign_go[i] & cre_gene_goterm$Divergence == 'Human', 'gene_symbol']
	allgo_gns = unique(allgo_gns)

	# Intersect with RAGs
	rags = toplot[toplot$Divergence == nodes[i] & toplot$is_sign == 'Sign', 'gene']
	gnstoplot = intersect(allgo_gns, rags)

	gnsL[[i]] = gnstoplot

}

gnstoplot = unique(unlist(gnsL))
toplot_sub = toplot[toplot$gene %in% gnstoplot,]

pdf('Adult_GO_RAGs.pdf', height = 5, width = 8)
ggbarplot(toplot_sub, x = 'Divergence', y = 'obs_peaks', color = 'is_sign', palette = c('lightblue', 'red'), fill = 'lightblue', size = 2) +
	ylab('Number of linked\n and divergent CREs') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	facet_wrap(~gene, nrow = 1) +
	NoLegend()
dev.off()


