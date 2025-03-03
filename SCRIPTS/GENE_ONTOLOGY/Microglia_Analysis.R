library(dplyr)
library(parallel)
library(Biostrings)
library(ggpubr)
library(Seurat)
library(reshape2)
library(readr)
library(pheatmap)
library(gridExtra)


plt_pref = 'Fetal_Microglia'

####
## FETAL MICROGLIA RAGs GO
####

# Load fetal cell type markers
ctmarks = readRDS('IN_DATA/PseudoBulk_DARs_MAJOR_MARKERS.RDS')
ctmarks = ctmarks[ctmarks$logFC < -1,]
micMarks = ctmarks[ctmarks$CellType == 'Microglia', 'Gene']

# Load evolutionary groups
evoGroups = readRDS('IN_DATA/Fetal_CREs_SignTested_ConsAdded.RDS')
human_cres = evoGroups[evoGroups$HumanSign == 'Sign', 'CRE']
hc_cres = evoGroups[evoGroups$HCSign == 'Sign', 'CRE']
hcgo_cres = evoGroups[evoGroups$HCGoSign == 'Sign', 'CRE']
gape_cres = evoGroups[evoGroups$Great_ApeSign == 'Sign', 'CRE']
ape_cres = evoGroups[evoGroups$ApeSign == 'Sign', 'CRE']
all_cres = Reduce(union, list(human_cres, hc_cres, hcgo_cres, gape_cres, ape_cres))

# Fetal peaks to genes
pks_to_gns = read_rds('IN_DATA/Fetal_LinksPeaksToGenes_AllStats.RDS')
pks_to_gns$peak2 = sub('-', ':', pks_to_gns$peak)

humanMicGns = pks_to_gns[pks_to_gns$peak2 %in% intersect(human_cres, micMarks), 'gene'] %>% unique
hcMicGns = pks_to_gns[pks_to_gns$peak2 %in% intersect(hc_cres, micMarks), 'gene'] %>% unique
hcgoMicGns = pks_to_gns[pks_to_gns$peak2 %in% intersect(hcgo_cres, micMarks), 'gene'] %>% unique
gapeMicGns = pks_to_gns[pks_to_gns$peak2 %in% intersect(gape_cres, micMarks), 'gene'] %>% unique
apeMicGns = pks_to_gns[pks_to_gns$peak2 %in% intersect(ape_cres, micMarks), 'gene'] %>% unique
allMicGns = pks_to_gns[pks_to_gns$peak2 %in% micMarks, 'gene'] %>% unique
allMicGns = setdiff(allMicGns, Reduce(union, list(humanMicGns, hcMicGns, hcgoMicGns, gapeMicGns, apeMicGns)))

####
## PLOT IN DEVELOPING MICROGLIA DATASET (POPOVA 2021)
####

seurObj = read_rds('IN_DATA/Popova_Microglia_snRNAseq.rds')
seurObj = NormalizeData(seurObj)

# Rename clusters as described in the paper
meta = seurObj_cortex@meta.data
meta$idents = Idents(seurObj_cortex)
meta$newannot = 'Unknown'
meta[meta$idents == 1, 'newannot'] = 'Axon tract\nassociated'
meta[meta$idents == 2, 'newannot'] = 'Homeostatic'
meta[meta$idents == 3, 'newannot'] = 'Ex vivo\nactivated'
meta[meta$idents == 4, 'newannot'] = 'Cytokine\nassociated'
meta[meta$idents == 5, 'newannot'] = 'Perivascular'
seurObj_cortex@meta.data = meta

# Remove ex vivo activated since it is a tissue handling artifact
seurObj_cortex = subset(seurObj_cortex, subset = newannot == 'Ex vivo\nactivated', invert = T)

# Rerun 2d embedding after the subset
seurObj_cortex = RunUMAP(seurObj_cortex, dims = 1:10, reduction = 'harmonySCT')

Idents(seurObj_cortex) = seurObj_cortex$newannot
pdf(paste0('PLOTS/', plt_pref, '_Microglia_UMAP_FetalBrain.pdf'))
print( DimPlot(seurObj_cortex, label = T, label.size = 7, raster = T) +
NoLegend() )
dev.off()

# Recreate seurat object for RNA assay
seurObj_cortex2 = CreateSeuratObject(counts = seurObj_cortex@assays$SCT@counts, meta.data = seurObj_cortex@meta.data)
seurObj_cortex2 = NormalizeData(seurObj_cortex2)

# Loop through the subtypes 
ctypes = unique(seurObj_cortex2$newannot)
fc_thresh = 0.25
dfL = list()
for(i in 1:length(ctypes)){

	# Human
	humanGns = intersect(humanMicGns, rownames(seurObj_cortex2))
	humanres = FindMarkers(seurObj_cortex2, ident.1 = ctypes[i], group.by = 'newannot', features = humanGns, only.pos = T)
	humanres = humanres[humanres$p_val_adj < 0.05 & humanres$avg_log2FC > fc_thresh,]
	humanratio = nrow(humanres) / length(humanGns)

	# HC
	hcGns = intersect(hcMicGns, rownames(seurObj_cortex2))
	hcres = FindMarkers(seurObj_cortex2, ident.1 = ctypes[i], group.by = 'newannot', features = hcGns, only.pos = T)
	hcres = hcres[hcres$p_val_adj < 0.05 & hcres$avg_log2FC > fc_thresh,]
	hcratio = nrow(hcres) / length(hcGns)

	# HCGo
	hcgoGns = intersect(hcgoMicGns, rownames(seurObj_cortex2))
	hcgores = FindMarkers(seurObj_cortex2, ident.1 = ctypes[i], group.by = 'newannot', features = hcgoGns, only.pos = T)
	hcgores = hcgores[hcgores$p_val_adj < 0.05 & hcgores$avg_log2FC > fc_thresh,]
	hcgoratio = nrow(hcgores) / length(hcgoGns)

	# Great ape
	gapeGns = intersect(gapeMicGns, rownames(seurObj_cortex2))
	gaperes = FindMarkers(seurObj_cortex2, ident.1 = ctypes[i], group.by = 'newannot', features = gapeGns, only.pos = T)
	gaperes = gaperes[gaperes$p_val_adj < 0.05 & gaperes$avg_log2FC > fc_thresh,]
	gaperatio = nrow(gaperes) / length(gapeGns)

	# Ape
	apeGns = intersect(apeMicGns, rownames(seurObj_cortex2))
	aperes = FindMarkers(seurObj_cortex2, ident.1 = ctypes[i], group.by = 'newannot', features = apeGns, only.pos = T)
	aperes = aperes[aperes$p_val_adj < 0.05 & aperes$avg_log2FC > fc_thresh,]
	aperatio = nrow(aperes) / length(apeGns)

	# All
	allGns = intersect(allMicGns, rownames(seurObj_cortex2))
	allres = FindMarkers(seurObj_cortex2, ident.1 = ctypes[i], group.by = 'newannot', features = allGns, only.pos = T)
	allres = allres[allres$p_val_adj < 0.05 & allres$avg_log2FC > fc_thresh,]
	allratio = nrow(allres) / length(allGns)

	dfL[[i]] = data.frame(degN = c(nrow(humanres), nrow(hcres), nrow(hcgores), nrow(gaperes), nrow(aperes), nrow(allres)),
				ratios = c(humanratio, hcratio, hcgoratio, gaperatio, aperatio, allratio),
				divergence = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'All'),
				ctype = ctypes[i])
}

df = do.call(rbind, dfL)
df$divergence = factor(df$divergence, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'All'))

pdf(paste0('PLOTS/', plt_pref, '_DivergentMicrogliaMarkers_CellTypeAssociation.pdf'), width = 16, height = 5)
print( ggbarplot(df, x = 'divergence', y = 'ratios', fill = 'divergence',
		color = 'divergence', position = position_dodge(0.9)) +
  theme_classic() + xlab('') + ylab('Ratio of subtype markers within\nmicroglia marker group') +
  theme(text = element_text(size=20, face = 'bold'), legend.pos = 'right') +
  facet_wrap(~ctype, scales = 'free', nrow = 1) +
  rotate_x_text(45) +
  NoLegend() )
dev.off()

# Calculate divergent vs non-divergent p-val with chi-square
statL = list()
for(i in 1:length(dfL)){

	qw = dfL[[i]]
	n1 = qw[qw$divergence != 'All', 'degN']
	n2 = qw[qw$divergence == 'All', 'degN']

	s1 = sum(n1/qw[qw$divergence != 'All', 'ratios'])
	s2 = n2/qw[qw$divergence == 'All', 'ratios']
	n1 = sum(n1)

	res = prop.test(c(n1,n2), c(s1,s2))

	statL[[i]] = data.frame(ctype = qw[1, 'ctype'], pval = res$p.val, odds_ratio = res$estimate[1] / res$estimate[2])
}

toexport = do.call(rbind, statL)

pdf(paste0('PLOTS/', plt_pref, '_MICROGLIA_SUBTYPE_EVO.pdf'))
print( grid.table(toexport) )
dev.off()

# All cytokine associated microglia markers
cytoMarks = FindMarkers(seurObj_cortex2, ident.1 = "Cytokine\nassociated", group.by = 'newannot', features = rownames(seurObj_cortex2), only.pos = T)
cytoMarkGns = cytoMarks[cytoMarks$p_val_adj < 0.05 & cytoMarks$avg_log2FC > 0.25,] %>% rownames

# All RDGs
all_div = read_rds('Fetal_Linked_Gene_Final_Results_Sign.RDS')
all_div = all_div[all_div$obs_peaks >= 5,]

# Cytokine-associated marker RDGs
cytoMarkRDGs = all_div[all_div$gene %in% cytoMarkGns, 'gene'] %>% unique
divMicAssocGns = Reduce(union, list(humanMicGns, hcMicGns, hcgoMicGns, gapeMicGns, apeMicGns))
divMicAssocRDGs = intersect(divMicAssocGns, cytoMarkRDGs) %>% unique

divMicAssocRDGs = cytoMarks[divMicAssocRDGs,] %>% .[order(.[, 'avg_log2FC'], decreasing = T),] %>% rownames

seurObj_cortex2$newannot = factor(seurObj_cortex2$newannot, rev(c("Cytokine\nassociated", "Axon tract\nassociated", "Homeostatic", "Perivascular")))

pdf(paste0(plt_pref, '_Microglia_CytokineMarker_RDGs.pdf'), height = 5, width = 6)
print( DotPlot(seurObj_cortex2, features = divMicAssocRDGs, group.by = 'newannot', dot.scale = 10) +
	theme(axis.text.x = element_text(size=20, face = 'italic'),
		axis.text.y = element_text(size=20, face = 'bold')) +
	xlab('') +
	rotate_x_text(45) )
dev.off()

