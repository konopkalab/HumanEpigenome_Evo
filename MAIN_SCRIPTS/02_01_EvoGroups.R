rm(list = ls())
library(rphast)
library(ape)
library(dplyr)
library(parallel)
library(Biostrings)
library(ggpubr)
library(seqinr)
library(phangorn)
library(msa)
source('Functions.R')

####
## Combine data
####

# Read and combine all substitutions
adultSubs = list.files(path = 'Ancestral_Substitutions', pattern = 'AncestralSubstitutions.RDS', full.names = T) %>% lapply(., function(x){readRDS(x)}) %>% do.call(rbind, .)

saveRDS(adultSubs, 'AdultSubstitutions.RDS')

####
## CONVERSION RATES BETWEEN AT AND GC
####

nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Anthropoids')
dfL = list()
for(i in 1:5){

	q1 = adultSubs[adultSubs$nodes == nodes[i],]
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

pdf('Adult_substitutions_conversion.pdf')
ggbarplot(toplot, x = 'node', y = 'vals', fill = 'vars', color = 'vars') +
ylab('Ratio') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90)
dev.off()


####
## CREATE CREs-Lineages substitution matrix
####


# Count the substitutions in CREs
allcounts = adultSubs %>% group_by(CRE, nodes) %>% group_keys %>% as.data.frame
allcounts$sizes = adultSubs %>% group_by(CRE, nodes) %>% group_size

countsBED = allcounts$CRE %>% strsplit(., ':|-') %>% do.call(rbind, .) %>% as.data.frame
countsBED$V2 = as.numeric(countsBED$V2)
countsBED$V3 = as.numeric(countsBED$V3)
allcounts$length = countsBED$V3 - countsBED$V2
allcounts$ratio = allcounts$sizes / allcounts$length
allcounts$nodes = factor(allcounts$nodes, levels = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'))

saveRDS(allcounts, 'AdultSubstitutions_Counted_NotNormalized.RDS')

# Plot total numbers
allcountsAgg = aggregate(sizes~nodes, data = allcounts, FUN = sum)

pdf('Adult_Brain_total_substitutions_raw.pdf')
ggbarplot(allcountsAgg, x = 'nodes', y = 'sizes', fill = 'nodes', color = 'nodes') +
ylab('Total number of substitutions') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

allcountsAgg$normFactors = c(6,2,8,4,13)
allcountsAgg$normSizes = allcountsAgg$sizes / allcountsAgg$normFactors

pdf('Adult_Brain_total_substitutions_normalized.pdf')
ggbarplot(allcountsAgg, x = 'nodes', y = 'normSizes', fill = 'nodes', color = 'nodes') +
ylab('Per MY number of substitutions') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()


# Adjust the ratios by branch length (Locke et al., 2011)
allcounts2 = allcounts
allcounts2$normRatio = allcounts2$ratio
allcounts2[allcounts2$nodes == 'Human', 'normRatio'] = allcounts2[allcounts2$nodes == 'Human', 'ratio'] / 6
allcounts2[allcounts2$nodes == 'HC', 'normRatio'] = allcounts2[allcounts2$nodes == 'HC', 'ratio'] / 2
allcounts2[allcounts2$nodes == 'HCGo', 'normRatio'] = allcounts2[allcounts2$nodes == 'HCGo', 'ratio'] / 8
allcounts2[allcounts2$nodes == 'Great_Ape', 'normRatio'] = allcounts2[allcounts2$nodes == 'Great_Ape', 'ratio'] / 4
allcounts2[allcounts2$nodes == 'Ape', 'normRatio'] = allcounts2[allcounts2$nodes == 'Ape', 'ratio'] / 13

saveRDS(allcounts2, 'AdultSubstitutions_Counted_BranchNormalized.RDS')

# Reshape to long format for grouping based on substitutions across branches
allcountsLong = dcast(allcounts2, CRE ~ nodes, value.var = 'normRatio')
allcountsLong[is.na(allcountsLong)] = 0
rownames(allcountsLong) = allcountsLong$CRE
allcountsLong$CRE = NULL
allcountsLong = as.matrix(allcountsLong)

# Find average number of substitutions per 1kb per million years
allcountsLong = allcountsLong * 1000
allcountsLong = as.data.frame(allcountsLong)

saveRDS(allcountsLong, 'AdultSubstitutions_Counted_BranchNormalized_LongFormat.RDS')

# Plot histogram of normalized substitutions
pdf('Number_of_subs_perMY_perKB_ADULT.pdf')
hist(melt(allcountsLong)$value, breaks = 1000, xlim = c(0,3))
dev.off()


toplot = melt(allcountsLong)
pdf('ADULT_Substitution_counts_branch_normalized.pdf', width = 4)
ggboxplot(toplot, x = 'variable', y = 'value', outlier.shape = NA, color = 'variable', ylim = c(0, 20)) +
ylab('Number of substitutions per MY per KB') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

# Plot without branch correction
tmpLong = dcast(allcounts2, CRE ~ nodes, value.var = 'ratio')
tmpLong[is.na(tmpLong)] = 0
rownames(tmpLong) = tmpLong$CRE
tmpLong$CRE = NULL
tmpLong = as.matrix(tmpLong)
tmpLongDF = tmpLong * 1000
tmpLongDF = as.data.frame(tmpLongDF)

toplot = melt(tmpLongDF)
pdf('ADULT_Substitution_counts_raw.pdf', width = 4)
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

# Set the background
bcg = lapply(1:1000, function(x){adultSubs[sample(1:nrow(adultSubs), 10000), 'nodes'] %>% table %>% prop.table}) %>% do.call(rbind, .)

# Find the proportions across lineages
cres = unique(adultSubs$CRE)
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape')
obsL = list()
obsL = mclapply(1:length(cres), mc.cores = 30, function(x){
	obs = adultSubs[adultSubs$CRE == cres[x], 'nodes'] %>% table %>% prop.table
	missing = nodes[!(nodes %in% names(obs))]
	obs = setNames(c(obs, rep(0, length(missing))), c(names(obs), missing))
	if(x%%1000 == 0){print(x)}
	obs = obs[nodes]
	t(obs) %>% as.data.frame %>% mutate(CRE = cres[x])
})

obsDF = do.call(rbind, obsL)
saveRDS(obsDF, 'obsDF.RDS')

# Statistically compare the proportions across lineages
bcg = bcg[,nodes]
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

saveRDS(pvalsDF, 'pvalsDF.RDS')

# Combine the p-values and the normalized substitution rates
allcountsLong = readRDS('AdultSubstitutions_Counted_BranchNormalized_LongFormat.RDS')
allcountsLong = allcountsLong[pvalsDF$CRE,] # for reordering, number of rows are the same
colnames(allcountsLong) = paste0(colnames(allcountsLong), '_norm')

# Add raw substitutions counts as well
allcounts2 = readRDS('AdultSubstitutions_Counted_BranchNormalized.RDS')

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

# Greater than 1 standard deviation more divergent than the rest of the CREs per node
# Significantly more divergent than the other nodes (FDR < 0.05)
thresh = 1

finalDF$HumanSign = 'NS'
finalDF[p.adjust(finalDF$hpval, method = 'BH') < 0.05 & finalDF$hFC > 1.5 &
	avHuman > thresh & finalDF$Human_raw >= 2, 'HumanSign'] = 'Sign'

finalDF$HCSign = 'NS'
finalDF[p.adjust(finalDF$hcpval, method = 'BH') < 0.05 & finalDF$hcFC > 1.5 &
	avHC > thresh & finalDF$HC_raw >= 2, 'HCSign'] = 'Sign'

finalDF$HCGoSign = 'NS'
finalDF[p.adjust(finalDF$hcgopval, method = 'BH') < 0.05 & finalDF$hcgoFC > 1.5 &
	avHCGo > thresh & finalDF$HCGo_raw >= 2, 'HCGoSign'] = 'Sign'

finalDF$Great_ApeSign = 'NS'
finalDF[p.adjust(finalDF$gapepval, method = 'BH') < 0.05 & finalDF$gapeFC > 1.5 &
	avGreatApe > thresh & finalDF$Great_Ape_raw >= 2, 'Great_ApeSign'] = 'Sign'

finalDF$ApeSign = 'NS'
finalDF[p.adjust(finalDF$apepval, method = 'BH') < 0.05 & finalDF$apeFC > 1.5 &
	avApe > thresh & finalDF$Ape_raw >= 2, 'ApeSign'] = 'Sign'


saveRDS(finalDF, 'AdultCREs_SignTested_SD_1_FC.RDS')


####
## Identify conserved CREs and plot groups
####

adultSub = readRDS('AdultCREs_SignTested_SD_1_FC.RDS')

humdivCRE = adultSub[adultSub$HumanSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
hcdivCRE = adultSub[adultSub$HCSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
hcgodivCRE = adultSub[adultSub$HCGoSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
gapedivCRE = adultSub[adultSub$Great_ApeSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
apedivCRE = adultSub[adultSub$ApeSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique

cons = finalDF[finalDF$Human_norm < mean(finalDF$Human_norm) &
	finalDF$HC_norm < mean(finalDF$HC_norm) &
	finalDF$HCGo_norm < mean(finalDF$HCGo_norm) &
	finalDF$Great_Ape_norm < mean(finalDF$Great_Ape_norm) &
	finalDF$Ape_norm < mean(finalDF$Ape_norm), 'CRE'] %>% gsub(':|-', '_', .) %>% unique

cons = setdiff(cons, Reduce(union, list(humdivCRE, hcdivCRE, hcgodivCRE, gapedivCRE, apedivCRE)))

finalDF$consSign = ifelse(gsub(':|-', '_', finalDF$CRE) %in% cons, 'Sign', 'NS')
saveRDS(finalDF, 'AdultCREs_SignTested_SD_1_FC_ConsAdded.RDS')

# Plot sizes of the groups
consCRE = finalDF[finalDF$consSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique
alldivs = Reduce(union, list(humdivCRE, hcdivCRE, hcgodivCRE, gapedivCRE, apedivCRE, consCRE))
toplot = data.frame(vars = c('Non-Divergent', 'Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Conserved'),
			vals = c(nrow(finalDF) - length(alldivs), length(humdivCRE), length(hcdivCRE), length(hcgodivCRE), length(gapedivCRE), length(apedivCRE), length(consCRE)))

pdf('GroupsNumbers_Adult.pdf')
ggbarplot(toplot, x = 'vars', y = 'vals', fill = 'vars', color = 'vars') +
ylab('Total number of CREs') + xlab('') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(45) +
NoLegend()
dev.off()


# Plot normalized substitution amounts per group
library(ggrastr)
groups = paste0(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons'), 'Sign')
nodes = paste0(c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape'), '_norm')
fcs = c('hFC', 'hcFC', 'hcgoFC', 'gapeFC', 'apeFC')

for(i in 1:length(groups)){

	toplot = finalDF[finalDF[groups[i]] == 'Sign',]
	toplot = melt(toplot[,c('Human_norm', 'HC_norm', 'HCGo_norm', 'Great_Ape_norm', 'Ape_norm')])

	pdf(paste0(groups[i], '_Adult_SubstitutionsGrouped.pdf'), width = 4)
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

creL = list(humdivCRE, hcdivCRE, hcgodivCRE, gapedivCRE, apedivCRE, consCRE)
names(creL) = c('HumanSubs', 'HCSubs', 'HCGoSubs', 'GreatApeSubs', 'ApeSubs', 'consSubs')

bcg = length(unique(adultSub$CRE))

# Enrichment
library(GeneOverlap)
resgom = newGOM(creL, creL, genome.size=bcg)
pvalmat = getMatrix(resgom, name="pval")
pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
pvalmelt = melt(pvalmat, value.name = 'pval')
pvalmelt$log10_FDR = -log10(pvalmelt$pval)

oddsmat = getMatrix(resgom, name="odds.ratio")
oddsmelt = melt(oddsmat, value.name = 'odds.ratio')

jacmat = getMatrix(resgom, name="Jaccard")
jacmelt = melt(jacmat, value.name = 'jaccard')

resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio, jac = jacmelt$jaccard)
resdf$log10_round = round(resdf$log10_FDR, digits = 2)
resdf$OR_round = round(resdf$OR, digits = 2)
resdf$fdr_Sc = formatC(resdf$pval, format = "e", digits = 2)
resdf$is_sign = ifelse(resdf$pval < 0.05 & resdf$OR > 1.3,  '*', '')

pdf('SubstitutionGroups_Similarity_Adult.pdf', height = 7, width = 9)
ggscatter(resdf, x = 'Var1', y = 'Var2', color = 'OR', size = 'log10_round') +
  geom_text(aes(label = is_sign), vjust = 0.8, colour = "darkgreen", size = 12 ) +
  labs(x="", y="") + scale_size_continuous(range = c(2,15)) +
  scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(45)
dev.off()


resdf2 = resdf[!(resdf$Var1 == resdf$Var2),]
resdf2$jacR = round(resdf2$jac, digits = 2)
pdf('SubstitutionGroups_Similarity_JAC_Adult.pdf', height = 7, width = 9)
ggscatter(resdf2, x = 'Var1', y = 'Var2', color = 'jac', size = 12) +
  geom_label(data = resdf2, aes(label=jacR), color="white",
		label.size = NA, fill = alpha(c("white"),0), fontface = 'bold') + 
  labs(x="", y="") +
  scale_color_gradient2(midpoint = 0.2, low = 'blue', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(45)
dev.off()


# Upset plot
library(UpSetR)

interList = list(c("ApeSubs"), c("GreatApeSubs"), c("HCGoSubs"), c("HCSubs"), c("HumanSubs"),
					c("HumanSubs","HCSubs"), c("HumanSubs","HCGoSubs"), c("HumanSubs","GreatApeSubs"), c("HumanSubs","ApeSubs"),
					c("HCSubs","HCGoSubs"), c("HCSubs","GreatApeSubs"), c("HCSubs","ApeSubs"),
					c("HCGoSubs","GreatApeSubs"), c("HCGoSubs","ApeSubs"),
					c("GreatApeSubs","ApeSubs"))

pdf('Divergence_Nodes_Overlap_Upset_ADULT.pdf', width = 10)
UpSetR::upset(fromList(creL), names(creL), order.by = "degree", keep.order = T,  intersections = interList, nsets = 5,
	number.angles = 30, point.size = 3.5, line.size = 2, 
	mainbar.y.label = "Number of CREs", sets.x.label = "Total CREs per group",
	text.scale = c(2, 1.5, 1.5, 1.5, 2, 1.5),
	decreasing = c(F,F))
dev.off()




