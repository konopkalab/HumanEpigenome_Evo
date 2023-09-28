rm(list = ls())
library(rphast)
library(ape)
library(dplyr)
library(parallel)
library(Biostrings)
library(ggpubr)
library(Seurat)
library(reshape2)
library(pheatmap)
source('Functions.R')

####
## EVO GROUP - CELL TYPE MARKER ENRICHMENT
####

# Load evolutionary groups
evoGroups = readRDS('AdultCREs_SignTested_SD_1_FC_ConsAdded.RDS')
human_cres = evoGroups[evoGroups$HumanSign == 'Sign', 'CRE']
hc_cres = evoGroups[evoGroups$HCSign == 'Sign', 'CRE']
hcgo_cres = evoGroups[evoGroups$HCGoSign == 'Sign', 'CRE']
gape_cres = evoGroups[evoGroups$Great_ApeSign == 'Sign', 'CRE']
ape_cres = evoGroups[evoGroups$ApeSign == 'Sign', 'CRE']
con_cres = evoGroups[evoGroups$consSign == 'Sign', 'CRE']

evoL = list(human_cres, hc_cres, hcgo_cres, gape_cres, ape_cres, con_cres)
names(evoL) = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Cons')

# Load cell type markers
ctmarks = readRDS('DATASETS/Adult_PseudoBulk_DARs_MAJOR_MARKERS.RDS')
ctmarks = ctmarks[ctmarks$logFC < -1,]
ctmarksL = split(ctmarks, ctmarks$CellType)
ctmarksL = lapply(ctmarksL, function(x){x$Gene})

# Run enrichment
matForDAR = readRDS('DATASETS/Adult_ATAC_pseudobulk_perBroadCellType_normalized.RDS')
bcg = nrow(matForDAR)

# Fisher's exact test enrichment and depletion
toplot = geneOvEnrDep(gnL1 = evoL, gnL2 = ctmarksL, bcg, plot = F, hg = 10, wd = 14)
toplot$sign_label = formatC(toplot$FDR, format = "e", digits = 0)

pdf('Adult_EvoGroup_CellTypeEnrichment.pdf', width = 10, height = 6)
ggscatter(toplot, x = 'var1', y = 'var2', color = 'OddsRatio', size = 'log10FDR') +
geom_label(data = toplot, aes(label = sign_label), color="black",
	label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
labs(x="", y="") +
scale_size_continuous(range = c(5,25)) +
scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red', limits = c(0.4,2)) +
theme_classic() +
theme(text = element_text(size=20, face = 'bold')) +
geom_text(aes(label = sign_label_2), vjust = 1.2, colour = "darkgreen", fontface = 'bold', size = 20 ) +
rotate_x_text(45)
dev.off()

# Load evolutionary groups
evoGroups = readRDS('FetalCREs_SignTested_SD_1_FC_ConsAdded.RDS')
human_cres = evoGroups[evoGroups$HumanSign == 'Sign', 'CRE']
hc_cres = evoGroups[evoGroups$HCSign == 'Sign', 'CRE']
hcgo_cres = evoGroups[evoGroups$HCGoSign == 'Sign', 'CRE']
gape_cres = evoGroups[evoGroups$Great_ApeSign == 'Sign', 'CRE']
ape_cres = evoGroups[evoGroups$ApeSign == 'Sign', 'CRE']
con_cres = evoGroups[evoGroups$consSign == 'Sign', 'CRE']

evoL = list(human_cres, hc_cres, hcgo_cres, gape_cres, ape_cres, con_cres)
names(evoL) = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Cons')

# Load cell type markers
ctmarks = readRDS('DATASETS/Fetal_PseudoBulk_DARs_MAJOR_MARKERS.RDS')
ctmarks = ctmarks[ctmarks$logFC < -1,]
ctmarks = ctmarks[ctmarks$CellType != 'Others',]
ctmarksL = split(ctmarks, ctmarks$CellType)
ctmarksL = lapply(ctmarksL, function(x){x$Gene})

# Run enrichment
matForDAR = readRDS('Fetal_ATAC_pseudobulk_perBroadCellType_normalized.RDS')
bcg = nrow(matForDAR)

# Fisher's exact test enrichment and depletion
toplot = geneOvEnrDep(gnL1 = evoL, gnL2 = ctmarksL, bcg, plot = F, hg = 10, wd = 14, fn = 'tmp')
toplot$sign_label = formatC(toplot$FDR, format = "e", digits = 0)

pdf('Fetal_EvoGroup_CellTypeEnrichment.pdf', width = 10, height = 6)
ggscatter(toplot, x = 'var1', y = 'var2', color = 'OddsRatio', size = 'log10FDR') +
geom_label(data = toplot, aes(label = sign_label), color="black",
	label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
labs(x="", y="") +
scale_size_continuous(range = c(10,25)) +
scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red', limits = c(0.3,1.4)) +
theme_classic() +
theme(text = element_text(size=20, face = 'bold')) +
geom_text(aes(label = sign_label_2), vjust = 1.2, colour = "darkgreen", fontface = 'bold', size = 20 ) +
rotate_x_text(45)
dev.off()


