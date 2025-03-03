require(rphast)
require(ape)
require(dplyr)
require(parallel)
require(Biostrings)
require(ggpubr)
require(Seurat)
require(reshape2)
require(pheatmap)
require(readr)

args = commandArgs(trailingOnly = TRUE)

for (arg in args) {
	split_arg <- strsplit(arg, "=")[[1]]
	var_name <- split_arg[1]
	var_value <- split_arg[2]
	print(var_name)
	print(var_value)
	assign(var_name, var_value)
}

# Nonconfigurable variables
branch_lengths = c(6,2,8,4,13)
names(branch_lengths) = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')
nodes = c(names(branch_lengths), 'Anthropoids')
plt_pref = paste0('PLOTS/', plt_pref)

dir.create(outdir, showWarnings = F)
dir.create('PLOTS', showWarnings = F)

####
## CONVERSION RATES BETWEEN AT AND GC
####

dfSubs = read_rds(paste0(subsdir, '/AncestralSubstitutions_ALL.RDS'))

dfL = list()
for(i in 1:5){

	q1 = dfSubs[dfSubs$nodes == nodes[i],]
	q1$change = paste0(q1[[ nodes[i] ]], q1[[ nodes[i+1] ]])
	table(q1$change) %>% prop.table

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

toplot = do.call(rbind, dfL)

pdf(paste0(plt_pref, '_substitutions_conversion.pdf'))
print( ggbarplot(toplot, x = 'node', y = 'vals', fill = 'vars', color = 'vars') +
ylab('Ratio') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) )
dev.off()


####
## CREATE CREs-Lineages substitution matrix
####

# Count the substitutions in CREs
allcounts = dfSubs %>% group_by(CRE, nodes) %>% group_keys %>% as.data.frame
allcounts$sizes = dfSubs %>% group_by(CRE, nodes) %>% group_size

countsBED = allcounts$CRE %>% strsplit(., ':|-') %>% do.call(rbind, .) %>% as.data.frame
countsBED$V2 = as.numeric(countsBED$V2)
countsBED$V3 = as.numeric(countsBED$V3)
allcounts$length = countsBED$V3 - countsBED$V2
allcounts$ratio = allcounts$sizes / allcounts$length
allcounts$nodes = factor(allcounts$nodes, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))

saveRDS(allcounts, paste0(outdir, '/Substitutions_Counted_NotNormalized.RDS'))

# Plot total numbers
allcountsAgg = aggregate(sizes~nodes, data = allcounts, FUN = sum)

pdf(paste0(plt_pref, '_Brain_total_substitutions_raw.pdf'))
print( ggbarplot(allcountsAgg, x = 'nodes', y = 'sizes', fill = 'nodes', color = 'nodes') +
ylab('Total number of substitutions') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend() )
dev.off()

allcountsAgg$normFactors = branch_lengths
allcountsAgg$normSizes = allcountsAgg$sizes / allcountsAgg$normFactors

pdf(paste0(plt_pref, '_Brain_total_substitutions_normalized.pdf'))
print( ggbarplot(allcountsAgg, x = 'nodes', y = 'normSizes', fill = 'nodes', color = 'nodes') +
ylab('Per MY number of substitutions') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend() )
dev.off()


# Adjust the ratios by branch length (Locke et al., 2011)
allcounts2 = allcounts
allcounts2$normRatio = allcounts2$ratio
allcounts2[allcounts2$nodes == 'Human', 'normRatio'] = allcounts2[allcounts2$nodes == 'Human', 'ratio'] / 6
allcounts2[allcounts2$nodes == 'HC', 'normRatio'] = allcounts2[allcounts2$nodes == 'HC', 'ratio'] / 2
allcounts2[allcounts2$nodes == 'HCGo', 'normRatio'] = allcounts2[allcounts2$nodes == 'HCGo', 'ratio'] / 8
allcounts2[allcounts2$nodes == 'Great_Ape', 'normRatio'] = allcounts2[allcounts2$nodes == 'Great_Ape', 'ratio'] / 4
allcounts2[allcounts2$nodes == 'Ape', 'normRatio'] = allcounts2[allcounts2$nodes == 'Ape', 'ratio'] / 13

saveRDS(allcounts2, paste0(outdir, '/Substitutions_Counted_BranchNormalized.RDS'))


# Reshape to long format for grouping based on substitutions across branches
allcountsLong = dcast(allcounts2, CRE ~ nodes, value.var = 'normRatio')
allcountsLong[is.na(allcountsLong)] = 0
rownames(allcountsLong) = allcountsLong$CRE
allcountsLong$CRE = NULL
allcountsLong = as.matrix(allcountsLong)

# Find average number of substitutions per 1kb per million years
allcountsLong = allcountsLong * 1000
allcountsLong = as.data.frame(allcountsLong)

saveRDS(allcountsLong, paste0(outdir, '/Substitutions_Counted_BranchNormalized_LongFormat.RDS'))

# Plot normalized substitutions levels per node
toplot = melt(allcountsLong)
pdf(paste0(plt_pref, '_substitution_counts_branch_normalized.pdf'), width = 4)
print( ggboxplot(toplot, x = 'variable', y = 'value', outlier.shape = NA, color = 'variable', ylim = c(0, 20)) +
ylab('Number of substitutions per MY per KB') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend() )
dev.off()

# Plot it without branch correction
tmpLong = dcast(allcounts2, CRE ~ nodes, value.var = 'ratio')
tmpLong[is.na(tmpLong)] = 0
rownames(tmpLong) = tmpLong$CRE
tmpLong$CRE = NULL
tmpLong = as.matrix(tmpLong)
tmpLongDF = tmpLong * 1000
tmpLongDF = as.data.frame(tmpLongDF)

toplot = melt(tmpLongDF)
pdf(paste0(plt_pref, '_substitution_counts_raw.pdf'), width = 4)
ggboxplot(toplot, x = 'variable', y = 'value', outlier.shape = NA, color = 'variable', ylim = c(0, 20)) +
ylab('Number of substitutions per KB') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

####
## Find CREs with significantly more substitutions in certain lineages
####

# Find the proportions across lineages
gres = unique(dfSubs$CRE)
obsL = list()
obsL = mclapply(1:length(gres), mc.cores = 30, function(x){
	obs = dfSubs[dfSubs$CRE == gres[x], 'nodes'] %>% table %>% prop.table
	missing = nodes[!(nodes %in% names(obs))]
	obs = setNames(c(obs, rep(0, length(missing))), c(names(obs), missing))
	if(x%%1000 == 0){print(x)}
	obs = obs[nodes]
	t(obs) %>% as.data.frame %>% mutate(CRE = gres[x])
})

obsDF = do.call(rbind, obsL)
saveRDS(obsDF, paste0(outdir, '/obsDF.RDS'))

# Set the background
subsize = min(nrow(dfSubs)/2, 10000) # Thus will always be 10000 for the full dataset
bcg = lapply(1:1000, function(x){dfSubs[sample(1:nrow(dfSubs), subsize), 'nodes'] %>% table %>% prop.table}) %>% do.call(rbind, .)

# Statistically compare the proportions across lineages
bcg = bcg[,nodes[nodes != 'Anthropoids']]
pvalsL = list()
for(i in 1:nrow(obsDF)){

	hpval = sum(bcg[, 'Human'] >= obsDF[i, 'Human']) / nrow(bcg)
	hcpval = sum(bcg[, 'HC'] >= obsDF[i, 'HC']) / nrow(bcg)
	hcgopval = sum(bcg[, 'HCGo'] >= obsDF[i, 'HCGo']) / nrow(bcg)
	gapepval = sum(bcg[, 'Great_Ape'] >= obsDF[i, 'Great_Ape']) / nrow(bcg)
	apepval = sum(bcg[, 'Ape'] >= obsDF[i, 'Ape']) / nrow(bcg)

	hFC = obsDF[i, 'Human'] / median(bcg[, 'Human'])
	hcFC = obsDF[i, 'HC'] / median(bcg[, 'HC'])
	hcgoFC = obsDF[i, 'HCGo'] / median(bcg[, 'HCGo'])
	gapeFC = obsDF[i, 'Great_Ape'] / median(bcg[, 'Great_Ape'])
	apeFC = obsDF[i, 'Ape'] / median(bcg[, 'Ape'])

	pvalsL[[i]] = data.frame(CRE = obsDF[i, 'CRE'], hpval = hpval,
			hcpval = hcpval, hcgopval = hcgopval,
			gapepval = gapepval, apepval = apepval, hFC = hFC,
			hcFC = hcFC, hcgoFC = hcgoFC,
			gapeFC = gapeFC, apeFC = apeFC)

	if(i%%1000 == 0){print(i)}
}

pvalsDF = do.call(rbind, pvalsL)

saveRDS(pvalsDF, paste0(outdir, '/pvalsDF.RDS'))

# Combine the p-values and the normalized substitution rates
allcountsLong = readRDS(paste0(outdir, '/Substitutions_Counted_BranchNormalized_LongFormat.RDS'))
allcountsLong = allcountsLong[pvalsDF$CRE,] # for reordering, number of rows are the same
colnames(allcountsLong) = paste0(colnames(allcountsLong), '_norm')

# Add raw substitutions counts as well
allcounts2 = readRDS(paste0(outdir, '/Substitutions_Counted_BranchNormalized.RDS'))

rawLong = dcast(allcounts2, CRE ~ nodes, value.var = 'sizes')
rawLong[is.na(rawLong)] = 0
rownames(rawLong) = rawLong$CRE
rawLong$CRE = NULL
rawLong = as.matrix(rawLong)
rawLong = rawLong[pvalsDF$CRE,] # for reordering, number of rows are the same
colnames(rawLong) = paste0(colnames(rawLong), '_raw')

####
## Assign evolutionary divergence per CRE
####

finalDF = cbind(pvalsDF, allcountsLong, rawLong)
avHuman = finalDF[,'Human_norm'] %>% scale
avHC = finalDF[,'HC_norm'] %>% scale
avHCGo = finalDF[,'HCGo_norm'] %>% scale
avGreatApe = finalDF[,'Great_Ape_norm'] %>% scale
avApe = finalDF[,'Ape_norm'] %>% scale

# Greater than user-defined standard deviation more divergent than the rest of the CREs per node
# Significantly more divergent than the other nodes

finalDF$HumanSign = 'NS'
finalDF[p.adjust(finalDF$hpval, method = 'BH') < fdr_cutoff & finalDF$hFC > fc_cutoff &
	avHuman > sd_cutoff & finalDF$Human_raw >= 2, 'HumanSign'] = 'Sign'

finalDF$HCSign = 'NS'
finalDF[p.adjust(finalDF$hcpval, method = 'BH') < fdr_cutoff & finalDF$hcFC > fc_cutoff &
	avHC > sd_cutoff & finalDF$HC_raw >= 2, 'HCSign'] = 'Sign'

finalDF$HCGoSign = 'NS'
finalDF[p.adjust(finalDF$hcgopval, method = 'BH') < fdr_cutoff & finalDF$hcgoFC > fc_cutoff &
	avHCGo > sd_cutoff & finalDF$HCGo_raw >= 2, 'HCGoSign'] = 'Sign'

finalDF$Great_ApeSign = 'NS'
finalDF[p.adjust(finalDF$gapepval, method = 'BH') < fdr_cutoff & finalDF$gapeFC > fc_cutoff &
	avGreatApe > sd_cutoff & finalDF$Great_Ape_raw >= 2, 'Great_ApeSign'] = 'Sign'

finalDF$ApeSign = 'NS'
finalDF[p.adjust(finalDF$apepval, method = 'BH') < fdr_cutoff & finalDF$apeFC > fc_cutoff &
	avApe > sd_cutoff & finalDF$Ape_raw >= 2, 'ApeSign'] = 'Sign'

saveRDS(finalDF, paste0(outdir, '/CREs_SignTested.RDS'))


####
## Identify conserved CREs and plot groups
####

humdivCRE = finalDF[finalDF$HumanSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
hcdivCRE = finalDF[finalDF$HCSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
hcgodivCRE = finalDF[finalDF$HCGoSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
gapedivCRE = finalDF[finalDF$Great_ApeSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
apedivCRE = finalDF[finalDF$ApeSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique

cons = finalDF[finalDF$Human_norm < mean(finalDF$Human_norm) &
	finalDF$HC_norm < mean(finalDF$HC_norm) &
	finalDF$HCGo_norm < mean(finalDF$HCGo_norm) &
	finalDF$Great_Ape_norm < mean(finalDF$Great_Ape_norm) &
	finalDF$Ape_norm < mean(finalDF$Ape_norm), 'CRE'] %>% gsub(':|-', '_', .) %>% unique

cons = setdiff(cons, Reduce(union, list(humdivCRE, hcdivCRE, hcgodivCRE, gapedivCRE, apedivCRE)))

finalDF$consSign = ifelse(gsub(':|-', '_', finalDF$CRE) %in% cons, 'Sign', 'NS')
saveRDS(finalDF, paste0(outdir, '/CREs_SignTested_ConsAdded.RDS'))


# Plot sizes of the groups
consCRE = finalDF[finalDF$consSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
alldivs = Reduce(union, list(humdivCRE, hcdivCRE, hcgodivCRE, gapedivCRE, apedivCRE, consCRE))
toplot = data.frame(vars = c('Non-Divergent', 'Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Conserved'),
			vals = c(nrow(finalDF) - length(alldivs), length(humdivCRE), length(hcdivCRE), length(hcgodivCRE), length(gapedivCRE), length(apedivCRE), length(consCRE)))

pdf(paste0(plt_pref, '_groupsNumbers.pdf'))
print( ggbarplot(toplot, x = 'vars', y = 'vals', fill = 'vars', color = 'vars') +
ylab('Total number of CREs') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(45) +
NoLegend() )
dev.off()


# Plot normalized substitution amounts per group
library(ggrastr)
groups = paste0(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons'), 'Sign')
fcs = c('hFC', 'hcFC', 'hcgoFC', 'gapeFC', 'apeFC')

for(i in 1:length(groups)){

	toplot = finalDF[finalDF[groups[i]] == 'Sign',]
	toplot = melt(toplot[,c('Human_norm', 'HC_norm', 'HCGo_norm', 'Great_Ape_norm', 'Ape_norm')])

	pdf(paste0(plt_pref, groups[i], '_SubstitutionsGrouped.pdf'), width = 4)
	print( ggboxplot(toplot, x = 'variable', y = 'value', outlier.shape = NA,
			color = 'variable', ylim = c(0,6)) +
	ylab('Substitutions per MY per KB') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(90) +
	NoLegend() )
	dev.off()

}


####
## Similarity between groups
####

greL = list(humdivCRE, hcdivCRE, hcgodivCRE, gapedivCRE, apedivCRE, consCRE)
names(greL) = c('HumanSubs', 'HCSubs', 'HCGoSubs', 'GreatApeSubs', 'ApeSubs', 'consSubs')

# Upset plot
library(UpSetR)

interList = list(c("ApeSubs"), c("GreatApeSubs"), c("HCGoSubs"), c("HCSubs"), c("HumanSubs"),
					c("HumanSubs","HCSubs"), c("HumanSubs","HCGoSubs"), c("HumanSubs","GreatApeSubs"), c("HumanSubs","ApeSubs"),
					c("HCSubs","HCGoSubs"), c("HCSubs","GreatApeSubs"), c("HCSubs","ApeSubs"),
					c("HCGoSubs","GreatApeSubs"), c("HCGoSubs","ApeSubs"),
					c("GreatApeSubs","ApeSubs"))

pdf(paste0(plt_pref, '_divergence_Nodes_Overlap_Upset.pdf'), width = 10)
print( UpSetR::upset(fromList(greL), names(greL), order.by = "degree", keep.order = T,  intersections = interList, nsets = 5,
	number.angles = 30, point.size = 3.5, line.size = 2, 
	mainbar.y.label = "Number of CREs", sets.x.label = "Total CREs per group",
	text.scale = c(2, 1.5, 1.5, 1.5, 2, 1.5),
	decreasing = c(F,F)) )
dev.off()
