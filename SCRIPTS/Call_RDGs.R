require(Signac)
require(Seurat)
require(Matrix)
require(EnsDb.Hsapiens.v86)
require(BSgenome.Hsapiens.UCSC.hg38)
require(dplyr)
require(readr)
require(tgutil)
require(tidyverse)
require(ggpubr)
require(parallel)
require(gridExtra)

args = commandArgs(trailingOnly = TRUE)

for (arg in args) {
	split_arg <- strsplit(arg, "=")[[1]]
	var_name <- split_arg[1]
	var_value <- split_arg[2]
	if(grepl(',', var_value)){
		var_value = strsplit(var_value, ',') %>% unlist()
	}
	print(var_name)
	print(var_value)
	assign(var_name, var_value)
}


# Read data
if(new_or_ms == 'MS'){
	pks_to_gns = read_rds('IN_DATA/Link_Peaks_To_Genes_Result.RDS')
} else if(new_or_ms == 'NEW') {
	fls = list.files(path = outdir, pattern = 'Multiome_randgns_links', full.names = T)
	pks_to_gns = lapply(fls, function(x){read_rds(x)}) %>% do.call(rbind, .)
} else{
	stop('Error: new_or_ms should be either NEW or MS')
}


nthreads = as.numeric(nthreads)
obs_peaksN = as.numeric(obs_peaksN)
all_peaksN = as.numeric(all_peaksN)
fdr = as.numeric(fdr)
dir.create(outdir, showWarnings = F)

# Further FDR cutoff of 0.05 although it usually does not filter any more links
pks_to_gns$FDR = p.adjust(pks_to_gns$pval, method = 'fdr')
pks_to_gns = pks_to_gns[pks_to_gns$FDR < 0.05,]

pks_to_gns$peak2 = sub('-', ':', pks_to_gns$peak)
pks_to_gns$gene = factor(pks_to_gns$gene)
gns_bcg = pks_to_gns$gene %>% unique


# All CREs
hvarsAll = read_rds('OUT_DATA/Evo_groups/CREs_SignTested_ConsAdded.RDS')
cre_all = hvarsAll$CRE
gns_cre_all = pks_to_gns[pks_to_gns$peak2 %in% cre_all, 'gene'] %>% table

####
## Human divergent CREs
####

number_of_perm = 10000
cre_test = hvarsAll[hvarsAll$HumanSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$HumanSign == 'Sign')

randL = mclapply(1:number_of_perm, mc.cores = nthreads, function(x){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', x))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	if(x%%100 == 0){print(x)}
	gns_cre_rand
})


rand_df = do.call(cbind, randL)

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
rand_df = data.table::as.data.table(rand_df)

resL = mclapply(1:length(gnstotest), mc.cores = nthreads, function(x){
	gn = gnstotest[x]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[x,] %>% as.numeric %>% median

	pval = sum(rand_df[x,] > obs) / number_of_perm
	if(x%%100 == 0){print(x)}
	data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
})

resdf = do.call(rbind, resL)
saveRDS(resdf, paste0(outdir, '/Linked_Gene_Stats_Human_perm10000_noPreFilter.RDS'))

####
## HC divergent CREs
####

cre_test = hvarsAll[hvarsAll$HCSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$HCSign == 'Sign')

randL = mclapply(1:number_of_perm, mc.cores = nthreads, function(x){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', x))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	if(x%%100 == 0){print(x)}
	gns_cre_rand
})


rand_df = do.call(cbind, randL)

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
rand_df = data.table::as.data.table(rand_df)

resL = mclapply(1:length(gnstotest), mc.cores = nthreads, function(x){
	gn = gnstotest[x]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[x,] %>% as.numeric %>% median

	pval = sum(rand_df[x,] > obs) / number_of_perm
	if(x%%100 == 0){print(x)}
	data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
})

resdf = do.call(rbind, resL)
saveRDS(resdf, paste0(outdir, '/Linked_Gene_Stats_HC_perm10000_noPreFilter.RDS'))

####
## HCGo divergent CREs
####

cre_test = hvarsAll[hvarsAll$HCGoSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$HCGoSign == 'Sign')
randL = mclapply(1:number_of_perm, mc.cores = nthreads, function(x){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', x))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	if(x%%100 == 0){print(x)}
	gns_cre_rand
})


rand_df = do.call(cbind, randL)

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
rand_df = data.table::as.data.table(rand_df)

resL = mclapply(1:length(gnstotest), mc.cores = nthreads, function(x){
	gn = gnstotest[x]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[x,] %>% as.numeric %>% median

	pval = sum(rand_df[x,] > obs) / number_of_perm
	if(x%%100 == 0){print(x)}
	data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
})

resdf = do.call(rbind, resL)
saveRDS(resdf, paste0(outdir, '/Linked_Gene_Stats_HCGo_perm10000_noPreFilter.RDS'))

####
## Great_Ape divergent CREs
####

cre_test = hvarsAll[hvarsAll$Great_ApeSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$Great_ApeSign == 'Sign')
randL = mclapply(1:number_of_perm, mc.cores = nthreads, function(x){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', x))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	if(x%%100 == 0){print(x)}
	gns_cre_rand
})


rand_df = do.call(cbind, randL)

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
rand_df = data.table::as.data.table(rand_df)

resL = mclapply(1:length(gnstotest), mc.cores = nthreads, function(x){
	gn = gnstotest[x]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[x,] %>% as.numeric %>% median

	pval = sum(rand_df[x,] > obs) / number_of_perm
	if(x%%100 == 0){print(x)}
	data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
})

resdf = do.call(rbind, resL)
saveRDS(resdf, paste0(outdir, '/Linked_Gene_Stats_GreatApe_perm10000_noPreFilter.RDS'))

####
## Ape divergent CREs
####

cre_test = hvarsAll[hvarsAll$ApeSign == 'Sign', 'CRE']
gns_cre_test = pks_to_gns[pks_to_gns$peak2 %in% cre_test, 'gene'] %>% table

# Gene links for random CREs
totdiv = sum(hvarsAll$ApeSign == 'Sign')
randL = mclapply(1:number_of_perm, mc.cores = nthreads, function(x){

	cre_rand = hvarsAll$CRE %>% sample(., totdiv)
	gns_cre_rand = pks_to_gns[pks_to_gns$peak2 %in% cre_rand, 'gene'] %>% table() %>% as.data.frame
	colnames(gns_cre_rand) = c('gene', paste0('Randomized_', x))
	rownames(gns_cre_rand) = gns_cre_rand$gene
	gns_cre_rand$gene = NULL
	if(x%%100 == 0){print(x)}
	gns_cre_rand
})


rand_df = do.call(cbind, randL)

# Genes that have significantly more links in the divergent group than the background
gnstotest = names(gns_cre_test)
rand_df = data.table::as.data.table(rand_df)

resL = mclapply(1:length(gnstotest), mc.cores = nthreads, function(x){
	gn = gnstotest[x]
	obs = gns_cre_test[gn]
	all = gns_cre_all[gn]
	rand_mean = rand_df[x,] %>% as.numeric %>% median

	pval = sum(rand_df[x,] > obs) / number_of_perm
	if(x%%100 == 0){print(x)}
	data.frame(gene = gn, pval = pval,
				obs_peaks = obs,
				bcg_peaks_mean = rand_mean,
				all_peaks = all)
})

resdf = do.call(rbind, resL)
saveRDS(resdf, paste0(outdir, '/Linked_Gene_Stats_Ape_perm10000_noPreFilter.RDS'))

####
## Read all results
####

human_resdf = read_rds(paste0(outdir, '/Linked_Gene_Stats_Human_perm10000_noPreFilter.RDS')) %>% add_column(Divergence = 'Human')
hc_resdf = read_rds(paste0(outdir, '/Linked_Gene_Stats_HC_perm10000_noPreFilter.RDS')) %>% add_column(Divergence = 'HC')
hcgo_resdf = read_rds(paste0(outdir, '/Linked_Gene_Stats_HCGo_perm10000_noPreFilter.RDS')) %>% add_column(Divergence = 'HCGo')
greatApe_resdf = read_rds(paste0(outdir, '/Linked_Gene_Stats_GreatApe_perm10000_noPreFilter.RDS')) %>% add_column(Divergence = 'Great_Ape')
ape_resdf = read_rds(paste0(outdir, '/Linked_Gene_Stats_Ape_perm10000_noPreFilter.RDS')) %>% add_column(Divergence = 'Ape')

# Add filter
human_resdf2 = human_resdf[human_resdf$obs_peaks >= obs_peaksN & human_resdf$all_peaks >= all_peaksN, ]
hc_resdf2 = hc_resdf[hc_resdf$obs_peaks >= obs_peaksN & hc_resdf$all_peaks >= all_peaksN, ]
hcgo_resdf2 = hcgo_resdf[hcgo_resdf$obs_peaks >= obs_peaksN & hcgo_resdf$all_peaks >= all_peaksN, ]
greatApe_resdf2 = greatApe_resdf[greatApe_resdf$obs_peaks >= obs_peaksN & greatApe_resdf$all_peaks >= all_peaksN, ]
ape_resdf2 = ape_resdf[ape_resdf$obs_peaks >= obs_peaksN & ape_resdf$all_peaks >= all_peaksN, ]

# Multiple testing correction
human_resdf2$FDR = p.adjust(human_resdf2$pval, method = 'fdr')
hc_resdf2$FDR = p.adjust(hc_resdf2$pval, method = 'fdr')
hcgo_resdf2$FDR = p.adjust(hcgo_resdf2$pval, method = 'fdr')
greatApe_resdf2$FDR = p.adjust(greatApe_resdf2$pval, method = 'fdr')
ape_resdf2$FDR = p.adjust(ape_resdf2$pval, method = 'fdr')


# Divergent genes. Only the ones with at least N peaks linked to it among all peaks.
fdr_cutoff = fdr

human_div = human_resdf2[human_resdf2$FDR < fdr_cutoff, ]

hc_div = hc_resdf2[hc_resdf2$FDR < fdr_cutoff, ]

hcgo_div = hcgo_resdf2[hcgo_resdf2$FDR < fdr_cutoff, ]

gape_div = greatApe_resdf2[greatApe_resdf2$FDR < fdr_cutoff, ]

ape_div = ape_resdf2[ape_resdf2$FDR< fdr_cutoff, ]


all_res = do.call(rbind, list(human_resdf, hc_resdf, hcgo_resdf, greatApe_resdf, ape_resdf))
all_div = do.call(rbind, list(human_div, hc_div, hcgo_div, gape_div, ape_div))

write_rds(all_res, paste0(outdir, '/Linked_Gene_Final_Results_All.RDS'))
write_rds(all_div, paste0(outdir, '/Linked_Gene_Final_Results_Sign.RDS'))


####
## HS-DEG ENRICHMENT
####

all_res = read_rds(paste0(outdir, '/Linked_Gene_Final_Results_All.RDS'))
all_div = read_rds(paste0(outdir, '/Linked_Gene_Final_Results_Sign.RDS'))

# All genes tested for DEGs
alldeg_res = readRDS('IN_DATA/PSEUDOBULK_DEGs_ALL.RDS')
all_gns_degres = alldeg_res$Gene %>% unique

# All genes tested for peak-link
seq_evo_bcg = read_rds('IN_DATA/pr_gns.RDS') %>% unlist %>% as.character

# Species-specific differentially expressed genes
hsdegs = alldeg_res[alldeg_res$Evolution == 'Human_Specific', 'Gene'] %>% unique

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
nperm = 10000

human_rand_ovs = sapply(1:nperm, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% human_div) / length(human_div)})
human_pval = sum(human_rand_ovs > human_obs_ov) / nperm

hc_rand_ovs = sapply(1:nperm, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% hc_div) / length(hc_div)})
hc_pval = sum(hc_rand_ovs > hc_obs_ov) / nperm

hcgo_rand_ovs = sapply(1:nperm, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% hcgo_div) / length(hcgo_div)})
hcgo_pval = sum(hcgo_rand_ovs > hcgo_obs_ov) / nperm

gape_rand_ovs = sapply(1:nperm, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% gape_div) / length(gape_div)})
gape_pval = sum(gape_rand_ovs > gape_obs_ov) / nperm

ape_rand_ovs = sapply(1:nperm, function(x){sum(sample(all_gns_degres, length(hsdegs)) %in% ape_div) / length(ape_div)})
ape_pval = sum(ape_rand_ovs > ape_obs_ov) / nperm

rands = rbind(human_rand_ovs, hc_rand_ovs, hcgo_rand_ovs, gape_rand_ovs, ape_rand_ovs) %>% as.data.frame
colnames(rands) = paste0('Randomized_', 1:ncol(rands))
rands$Observed = c(human_obs_ov, hc_obs_ov, hcgo_obs_ov, gape_obs_ov, ape_obs_ov)
rands$ids = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')

toplot = reshape2::melt(rands)
toplot$ids = factor(toplot$ids, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))

pdf(plt_pref, '_HSDEG_RDG_Enrichment.pdf'))
print( ggplot(toplot, aes(x = ids, y = value)) +
geom_boxplot(outlier.shape = NA, color = c('blue', 'lightblue', 'purple', 'pink', 'peachpuff'),lwd=1.5) +
theme_classic() +
theme(text = element_text(size=20)) +
xlab('') + ylab('Overlap_of_peakLinkedGenes\nand_HS-DEGs') +
rotate_x_text(90) +
geom_point(data = toplot[toplot$variable == 'Observed',], color = 'red', size = 3) +
theme(strip.text.x = element_text(angle = 90)) +
ylim(c(0,0.4)) )
dev.off()

# Export p-value table
toexport = toplot[toplot$variable == 'Observed',]
toexport$pval = c(human_pval, hc_pval, hcgo_pval, gape_pval, ape_pval)
toexport$FDR = p.adjust(toexport$pval, method = 'BH')

pdf(paste0(plt_pref, '_HSDEG_RDG_Enrichment_TABLE.pdf'))
print( grid.table(toexport) )
dev.off()

####
## PLOT RDGs
####

# Plot number of genes per category
groups = unique(all_div[['Divergence']])
vals = table(all_div[['Divergence']])[groups] %>% as.numeric
toplot = data.frame(vars = groups, vals = vals)
colors = c('blue', 'lightblue', 'purple', 'pink', 'red')

pdf(paste0(plt_pref, '_Number_of_divergent_genes.pdf'))
print( ggbarplot(toplot, x = 'vars', y = 'vals', color = 'vars', fill = 'vars', palette = colors) +
	ylab('Number of significantly linked genes') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	NoLegend() )
dev.off()

# Per lineage
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')
for(i in 1:length(nodes)){

	gns = all_div[all_div$Divergence == nodes[i], 'gene'] %>% unique
	toplot = all_res[all_res$gene %in% gns,]
	toplot$ratToAll = toplot$obs_peaks / toplot$all_peaks
	toplot$ratToAll = ifelse(toplot$ratToAll > 0.75, 0.75, toplot$ratToAll)

	pdf(paste0(nodes[i], '_RAGs_Adult.pdf'), height = 15)
	print( ggscatter(toplot, x = 'Divergence', y = 'gene', color = 'ratToAll', size = 'obs_peaks') +
		ylab('Number of significantly linked genes') + xlab('') +
		scale_size_continuous(range = c(2,10)) +
		theme(text = element_text(size=20), legend.pos = 'right') +
		theme(axis.text.x = element_text(size=20),
				axis.text.y = element_text(size=20, face = c('italic')),
				axis.title = element_text(size=20)) +
		scale_color_gradient2(midpoint = 0.1, low = 'blue', mid = 'white',
				high = 'red', limits = c(0,0.75)) +
		rotate_x_text(45) )
	dev.off()

}
