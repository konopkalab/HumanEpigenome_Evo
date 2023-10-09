rm(list = ls())
library(dplyr)
library(parallel)
library(Biostrings)
library(ggpubr)
library(Seurat)
library(reshape2)
library(readr)
library(pheatmap)
source('utility_functions.R')

####
## FETAL MICROGLIA RAGs GO
####

# Load fetal cell type markers
ctmarks = readRDS('CELLTYPE_MARKERS/FETAL/PseudoBulk_DARs_MAJOR_MARKERS.RDS')
ctmarks = ctmarks[ctmarks$logFC < -1,]
micMarks = ctmarks[ctmarks$CellType == 'Microglia',]

# Load evolutionary groups
evoGroups = readRDS('02_GroupCREs/FETAL/FetalCREs_SignTested_SD_1_FC_ConsAdded.RDS')
human_cres = evoGroups[evoGroups$HumanSign == 'Sign', 'CRE']
hc_cres = evoGroups[evoGroups$HCSign == 'Sign', 'CRE']
hcgo_cres = evoGroups[evoGroups$HCGoSign == 'Sign', 'CRE']
gape_cres = evoGroups[evoGroups$Great_ApeSign == 'Sign', 'CRE']
ape_cres = evoGroups[evoGroups$ApeSign == 'Sign', 'CRE']
all_cres = Reduce(union, list(human_cres, hc_cres, hcgo_cres, gape_cres, ape_cres))

# Load regulation accelerated genes
all_rags = read_rds('MULTIOME_FETAL/Linked_Gene_Final_Results_Sign_FDR.RDS')
pks_to_gns = read_rds('MULTIOME_FETAL/PKS_TO_GNS/Trevino_LinkPeaksToGenes_FINAL.RDS')
pks_to_gns$peak2 = sub('-', ':', pks_to_gns$peak)

# Evolutionarily divergent fetal microglia markers
all_cres_mic = intersect(micMarks$Gene, all_cres)

# Genes associated with these markers
micLinkedGns = pks_to_gns[pks_to_gns$peak2 %in% all_cres_mic, 'gene'] %>% unique

# Divergent genes
micLinkedRAGs = intersect(micLinkedGns, unique(all_rags$gene))

# Run gene ontology enrichment
rnadf = GOenrich(micLinkedRAGs, micLinkedGns, pCut = 1, qCut = 1, species = 'human')
write_rds(rnadf, 'MULTIOME_FETAL/GO_Enrich_MicrogliaLinkedRAGs.RDS')

# Extract significant results
rnadf_sign = rnadf[rnadf$p.adjust < 0.05, ]
rnadf_sign_fetal = rnadf_sign


# PLOT TERM GENES
ind = 2
term_name = rnadf_sign_fetal$Description[ind]

gogns = strsplit(rnadf_sign_fetal$geneID[ind], '/') %>% unlist
all_ragsL = split(all_rags, all_rags$Divergence)
all_ragsL = lapply(all_ragsL, function(x){x$gene})
dfL = list()
for(i in 1:length(all_ragsL)){
	df = sapply(gogns, function(x){sum(x %in% all_ragsL[[i]])}) %>% as.data.frame
	colnames(df) = 'Presence'
	df$node = names(all_ragsL)[i]
	df$gene = rownames(df)
	dfL[[i]] = df
}

toplot = do.call(rbind, dfL)
ord_gns = toplot[toplot$Presence == 1, 'gene'] %>% table %>% sort %>% names
toplot$gene = factor(toplot$gene, levels = ord_gns)
toplot$node = factor(toplot$node, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))

pdf(paste0(term_name, '_Fetal_Genes.pdf'), width = 7, height = 7)
ggscatter(toplot, y = 'gene', x = 'node', size = 'Presence') +
  scale_size_continuous(range = c(0,15)) +
  labs(x="", y="") +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold.italic')) +
  scale_x_discrete(position = "top") +
  rotate_x_text(-45)
dev.off()


#seurObj = read_rds('RNASEQ_FETAL/seurObj_clustered.RDS')
#seurObj = NormalizeData(seurObj)
#seurObj = subset(seurObj, subset = Major == 'Others', invert = T)
#DotPlot(seurObj, features = c(gogns), group.by = 'Major')


####
## PLOT IN DEVELOPING MICROGLIA DATASET (POPOVA 2021)
####

seurObj = read_rds('EXTERNAL_DATASETS/POPOVA_2021/biccn_cb.rds')
seurObj = NormalizeData(seurObj)

human_mic_rags = c('DOCK8', 'PLEC', 'CCL3', 'CCL4', 'CCL3L3', 'ADAM8', 'CRK')
seurObj_cortex = subset(seurObj, subset = Area == 'Cortex')
seurObj_cortex = subset(seurObj_cortex, subset = orig.ident != 'cerebellum_31')

# Rerun 2d embedding after the subset
seurObj_cortex = RunUMAP(seurObj_cortex, dims = 1:10, reduction = 'harmonySCT')

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

seurObj_cortex2 = seurObj_cortex
Idents(seurObj_cortex2) = seurObj_cortex2$newannot
pdf('Microglia_UMAP_FetalCortex.pdf')
DimPlot(seurObj_cortex2, label = T, label.size = 7, raster = T) +
NoLegend()
dev.off()

pdf('Microglia_HumanRAGs_FetalCortex.pdf', height = 6, width = 10)
DotPlot(seurObj_cortex2, features = human_mic_rags, group.by = 'newannot') +
	theme(axis.text.x = element_text(size=20, face = 'italic'),
		axis.text.y = element_text(size=20, face = 'bold')) +
	rotate_x_text(45)
dev.off()



####
## ADULT MICROGLIA RAGs GO
####

# Load adult cell type markers
ctmarks = readRDS('CELLTYPE_MARKERS/ADULT/PseudoBulk_DARs_MAJOR_MARKERS.RDS')
ctmarks = ctmarks[ctmarks$logFC < -1,]
micMarks = ctmarks[ctmarks$CellType == 'Microglia',]

# Load evolutionary groups
evoGroups = readRDS('02_GroupCREs/ADULT/AdultCREs_SignTested_SD_1_FC_ConsAdded.RDS')
human_cres = evoGroups[evoGroups$HumanSign == 'Sign', 'CRE']
hc_cres = evoGroups[evoGroups$HCSign == 'Sign', 'CRE']
hcgo_cres = evoGroups[evoGroups$HCGoSign == 'Sign', 'CRE']
gape_cres = evoGroups[evoGroups$Great_ApeSign == 'Sign', 'CRE']
ape_cres = evoGroups[evoGroups$ApeSign == 'Sign', 'CRE']
all_cres = Reduce(union, list(human_cres, hc_cres, hcgo_cres, gape_cres, ape_cres))

# Load regulation accelerated genes
all_rags = read_rds('MULTIOME_ADULT/Linked_Gene_Final_Results_Sign_FDR.RDS')
pks_to_gns = read_rds('MULTIOME_ADULT/PKS_TO_GNS/Ma_LinkPeaksToGenes_FINAL.RDS')
pks_to_gns$peak2 = sub('-', ':', pks_to_gns$peak)

# Evolutionarily divergent adult microglia markers
all_cres_mic = intersect(micMarks$Gene, all_cres)

# Genes associated with these markers
micLinkedGns = pks_to_gns[pks_to_gns$peak2 %in% all_cres_mic, 'gene'] %>% unique

# Divergent genes
micLinkedRAGs = intersect(micLinkedGns, unique(all_rags$gene))

# Run gene ontology enrichment
rnadf = GOenrich(micLinkedRAGs, micLinkedGns, pCut = 1, qCut = 1, species = 'human')
rnadf_adult = rnadf


####
## PLOT
####

# Combine both
rnadf_adult$type = 'Adult'
rnadf_sign_fetal$type = 'Fetal'
toplot = rbind(rnadf_adult, rnadf_sign_fetal)

# Selected GO terms
terms = rnadf_sign_fetal$Description
terms = terms[terms != 'intermediate filament cytoskeleton organization']

# Extract selected GO terms
plotenr = toplot[toplot$Description %in% terms,]

# For plotting
plotenr$log10FDR = -log10(plotenr$qvalue)
plotenr$is_sign = ifelse(plotenr$qvalue < 0.05, '*', '')
plotenr$gonameY = sapply(plotenr$Description, function(x){wrapper(x, width = 20)})

pdf('Microglia_All_RAGs_Gene_Ontology.pdf', width = 6, height = 5)
ggplot(plotenr, aes(x = type, y = gonameY)) +
  geom_point(aes(size = log10FDR, fill = OddsRatio), shape = 21, colour = "black") +
  scale_fill_gradient2(midpoint = 1.3, low = 'blue', high = 'red') +
  theme_classic() + xlab('') + ylab('') +
  theme(text = element_text(size=20, face = 'bold'), legend.pos = 'right') +
  rotate_x_text(45) +
  geom_text(aes(label = is_sign), vjust = 0.8, colour = "blue", size = 12 ) +
  scale_size_continuous(range = c(7,20))
dev.off()














