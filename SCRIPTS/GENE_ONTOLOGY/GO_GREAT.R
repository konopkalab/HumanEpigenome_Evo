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
library(rGREAT)
library(pheatmap)
library(readr)
library(rrvgo)
library("org.Hs.eg.db")

####
## Prepare for Gene Ontology Enrichments
####

# Set outdir
outdir = 'OUT_DATA/GREAT'

# Set prefix
pref = 'Adult'

# Read data
hvarsAll = readRDS('OUT_DATA/Evo_groups/CREs_SignTested_ConsAdded.RDS')

# Keep the ones linked to genes
pks_to_gns = read_rds('OUT_DATA/RDG/Multiome_randgns_links_All.RDS')
pks_to_gns$peak2 = sub('-', ':', pks_to_gns$peak)
hvarsAll_Linked = hvarsAll[hvarsAll$CRE %in% pks_to_gns$peak2,]

# Extract background regions
bcgCREs = hvarsAll_Linked$CRE
bcg_bed = bcgCREs %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(chr = as.character(V1), start = as.numeric(as.character(V2)),
			end = as.numeric(as.character(V3))) %>%
			makeGRangesFromDataFrame

# Remove ungrouped CREs
hvars = hvarsAll_Linked[!(hvarsAll_Linked$HumanSign == 'NS' & hvarsAll_Linked$HCSign == 'NS' &
		hvarsAll_Linked$HCGoSign == 'NS' & hvarsAll_Linked$Great_ApeSign == 'NS' &
		hvarsAll_Linked$ApeSign == 'NS' & hvarsAll_Linked$consSign == 'NS'), ]


####
## Find Gene Ontology Enrichments
####

groups = paste0(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons'), 'Sign')
gns_pks_all_L = list()
tableAll_L = list()
for(i in 1:length(groups)){

	name = gsub('Sign', '', groups[i])

	# Extract regions to test
	testCREs = hvars[hvars[groups[i]] == 'Sign', 'CRE']
	testCREs = unique(testCREs)
	test_bed = testCREs %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
			do.call(rbind, .) %>% as.data.frame() %>%
			mutate(chr = as.character(V1), start = as.numeric(as.character(V2)),
				end = as.numeric(as.character(V3))) %>%
			makeGRangesFromDataFrame

	# Run enrichment locally
	resMF = great(test_bed, "GO:MF", biomart_dataset = 'hsapiens_gene_ensembl', background = bcg_bed)
	resBP = great(test_bed, "GO:BP", biomart_dataset = 'hsapiens_gene_ensembl', background = bcg_bed)
	resCC = great(test_bed, "GO:CC", biomart_dataset = 'hsapiens_gene_ensembl', background = bcg_bed)

	resAll = c(resMF, resBP, resCC)
	names(resAll) = c('MF', 'BP', 'CC')
	saveRDS(resAll, paste0(outdir, '/GO_resAll_', name, '.RDS'))

	# Get enrichment tables
	table_MF = getEnrichmentTable(resMF) %>% add_column(type = 'MF') %>% add_column(Divergence = name)
	table_BP = getEnrichmentTable(resBP) %>% add_column(type = 'BP') %>% add_column(Divergence = name)
	table_CC = getEnrichmentTable(resCC) %>% add_column(type = 'CC') %>% add_column(Divergence = name)

	tableAll = do.call(rbind, list(table_MF, table_BP, table_CC)) %>% add_column(Divergence = name)
	tableAll_L[[i]] = tableAll

	# Find genes for significant terms
	enriched_names = tableAll[tableAll$p_adjust < 0.05 & tableAll$fold_enrichment > 1.3, 'description']
	enriched_ids = tableAll[tableAll$p_adjust < 0.05 & tableAll$fold_enrichment > 1.3, 'id']

	if(length(enriched_ids) == 0){
		next
	} else{

		gns_pks_L = list()
		for(j in 1:length(enriched_ids)){

			gotype = tableAll[tableAll$id == enriched_ids[j], 'type']
			tmpres = getRegionGeneAssociations(resAll[[gotype]], term_id = enriched_ids[j])
			tmpres$annotated_genes = sapply(tmpres$annotated_genes, function(x){x[1]})

			tmpres$gene_symbol = mapIds(org.Hs.eg.db, keys = as.character(tmpres$annotated_genes),
						keytype = "ENSEMBL", column="SYMBOL")
			tmpres = as.data.frame(tmpres)
			tmpres$id = enriched_ids[j]
			tmpres$name = enriched_names[j]

			gns_pks_L[[j]] = tmpres
			if(j%%10==0){print(j)}
		}

		gns_pks_df = do.call(rbind, gns_pks_L)
		gns_pks_df$seqnames = paste0('chr', gns_pks_df$seqnames)
		gns_pks_df$Divergence = name

		gns_pks_all_L[[i]] = gns_pks_df

		print(i)
	}
}


allres = do.call(rbind, tableAll_L)
saveRDS(allres, paste0(outdir, '/GO_allFinalDF.RDS'))

go_assoc_gns_pks = do.call(rbind, gns_pks_all_L)
saveRDS(go_assoc_gns_pks, paste0(outdir, '/GO_Sign_Assoc_Gns_Pks.RDS'))

####
## Find Significant Enrichments
####

# Read datasets
all_div = read_rds('OUT_DATA/RDG/Linked_Gene_Final_Results_Sign.RDS')
all_res = read_rds('OUT_DATA/RDG/Linked_Gene_Final_Results_All.RDS')
gores = readRDS(paste0(outdir, '/GO_allFinalDF.RDS'))
cre_gene_goterm = readRDS(paste0(outdir, '/GO_Sign_Assoc_Gns_Pks.RDS'))

# Find significant GO enrichments
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons')
sign_goL = list()
for(j in 1:length(nodes)){

	# Divergent genes and related stats
	sp_all = all_res[all_res$Divergence == nodes[j], ]
	sp_div = all_div[all_div$Divergence == nodes[j], 'gene']

	# Enriched go terms
	sp_gonames = gores[gores$p_adjust < 0.05 & gores$fold_enrichment > 1.3 & gores$Divergence == nodes[j], 'description']

	if(length(sp_gonames) == 0){next}

	# Find number of divergent genes and cres per go term enrichment
	totgnsL = list()
	totcreL = list()
	for(i in 1:length(sp_gonames)){
		gns = cre_gene_goterm[cre_gene_goterm$name == sp_gonames[i], 'gene_symbol'] %>% unique
		gns_inters = intersect(gns, sp_div)
		totcreL[[i]] = sp_all[sp_all$gene %in% gns_inters, 'obs_peaks'] %>% sum
		totgnsL[[i]] = length(gns_inters)
	}

	# Keep only the ones with N divergent genes and M divergent CREs.
	sign_goL[[j]] = sp_gonames[unlist(totgnsL) >= 1 & unlist(totcreL) >= 10]

	print(nodes[j])
	if(sum(unlist(totgnsL)) > 0){print(summary(unlist(totgnsL)))}
	if(sum(unlist(totcreL)) > 0){print(summary(unlist(totcreL)))}
}

names(sign_goL) = nodes

# Save final results
signs = unlist(sign_goL) %>% unique
gores_final = gores[gores$description %in% signs,]
gores_final$Divergence.1 = NULL
saveRDS(gores_final, paste0(outdir, '/GO_allFinalDF_Filtered.RDS'))


####
## Plot Number of Significant Enrichments
####

toplot = sapply(sign_goL, function(x){length(x)}) %>% as.data.frame
colnames(toplot)[1] = 'gonumber'
toplot$nodes = rownames(toplot)

pdf(paste0('PLOTS/', pref, '_Adult_GO_Enrich_TotNumber.pdf'))
print( ggbarplot(toplot, x = 'nodes', y = 'gonumber',
	color = 'nodes', fill = 'nodes') +
	ylab('Number of significant GO terms') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	NoLegend() )
dev.off()


####
## Plot Significant Enrichments with Revigo
####

gotypes = c('BP', 'MF', 'CC')
node = 'Human'

for(i in 1:length(gotypes)){

	goids = gores[gores$description %in% sign_goL[[node]] &
			gores$Divergence == node &
			gores$type == gotypes[i], 'id']

	goscores = gores[gores$description %in% sign_goL[[node]] &
			gores$Divergence == node &
			gores$type == gotypes[i], 'p_adjust'] %>% -log10(.)
	names(goscores) = goids

	simMatrix = calculateSimMatrix(goids,
		                        orgdb="org.Hs.eg.db",
		                        ont=gotypes[i],
		                        method="Rel")


	reducedTerms = reduceSimMatrix(simMatrix,
		                        goscores,
		                        threshold=0.7,
		                        orgdb="org.Hs.eg.db")

	pdf(paste0('PLOTS/', pref, '_Revigo_', node, '_', gotypes[i], '.pdf'))
	print( treemapPlot(reducedTerms) )
	dev.off()
}



####
## Plot Significant Enrichments
####

# Load adult Great results
fetal_gores = readRDS(paste0(outdir, '/Fetal_GO_allFinalDF.RDS'))
sign_go = unlist(sign_goL)

# Plot the significant terms
pltL = list()
for(i in 1:length(sign_go)){

	goname = sign_go[i]

	# Adult
	plotenr = gores[gores$description ==  goname,]
	plotenr$log10FDR = -log10(plotenr$p_adjust)
	plotenr$Divergence = factor(plotenr$Divergence, c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons'))
	plotenr$is_sign = ifelse(plotenr$p_adjust < 0.05, '*', '')
	plotenr$dev = 'Adult'
	plotenr_adult = plotenr

	# Fetal
	plotenr = fetal_gores[fetal_gores$description ==  goname,]
	plotenr$log10FDR = -log10(plotenr$p_adjust)
	plotenr$Divergence = factor(plotenr$Divergence, c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons'))
	plotenr$is_sign = ifelse(plotenr$p_adjust < 0.05, '*', '')
	plotenr$dev = 'Fetal'
	plotenr_fetal = plotenr

	# Combine
	plotenr = rbind(plotenr_fetal, plotenr_adult)

	wrapper <- function(x, ...) 
	{
	  paste(strwrap(x, ...), collapse = "\n")
	}

	plotenr$gonameY = sapply(plotenr$description, function(x){wrapper(x, width = 20)})

	pltL[[i]] = ggplot(plotenr, aes(x = Divergence, y = gonameY)) +
	  geom_point(aes(size = log10FDR, fill = fold_enrichment), shape = 21, colour = "black") +
	  scale_fill_gradient2(midpoint = 1, low = 'blue', high = 'red') +
	  theme_classic() + xlab('') + ylab('') +
	  theme(text = element_text(size=20, face = 'bold'), legend.pos = c('right'))+
	  rotate_x_text(45) +
	  facet_wrap(~dev) +
	  geom_text(aes(label = is_sign), vjust = 0.8, colour = "blue", size = 12 ) +
	  scale_size_continuous(range = c(1,10))
}

ggexport(pltL, filename = paste0('PLOTS/', pref, '_GO_Enrich_ALL.pdf'), width = 10, height = 5)

# Plot the significant terms

# RDGs
human_rdg = all_div[all_div$Divergence == 'Human', 'gene']

sign_goH = sign_go[grepl('Human', names(sign_go))]
pltL = list()
toplotL = list()
for(i in 1:length(sign_goH)){

	goname = sign_goH[i]

	# Plot the genes and number of linked CREs
	gns = cre_gene_goterm[cre_gene_goterm$name == goname, 'gene_symbol'] %>% unique
	gns_rdg = intersect(gns, human_rdg)
	toplot = all_res[all_res$gene %in% gns,]
	toplot = toplot[toplot$Divergence == 'Human',]
	toplot = toplot[order(toplot$obs_peaks, decreasing = T),]
	toplot = toplot[toplot$obs_peaks > 0,]

	toplot$is_rdg = ifelse(toplot$gene %in% gns_rdg, 'RDG', 'Non-RDG')

	wrapper <- function(x, ...) 
	{
	  paste(strwrap(x, ...), collapse = "\n")
	}

	toplot$gonameY = wrapper(goname, width = 20)

	pdf(paste0('PLOTS/', goname, '_', pref, '_GO_assoc_genes.pdf'), width = 10)
	print( ggbarplot(toplot, x = 'gene', y = 'obs_peaks',
				color = 'is_rdg', fill = 'is_rdg', palette = c('gray', 'blue')) +
			ylab('Number of linked divergent CREs') + xlab('') +
			theme(text = element_text(size=20)) +
			theme(axis.text.x = element_text(size=20),
					axis.text.y = element_text(size=20),
					axis.title = element_text(size=20)) +
			rotate_x_text(45) )
	dev.off()
}

