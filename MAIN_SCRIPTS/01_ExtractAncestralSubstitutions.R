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
## Extract the peaks from the alignment
####

# Load the tree
treetop = read.table("UCSC_MULTIZ30_2017/treeforR.txt")
treetop = treetop$V1
tree = read.tree(text = treetop)

# Keep only given species
apes = c('hg38', 'panTro5', 'gorGor5', 'ponAbe2', 'nomLeu3')
owm = c('rheMac8', 'macFas5', 'macNem1', 'cerAty1','papAnu3','chlSab2','manLeu1', 'nasLar1','colAng1','rhiRox1','rhiBie1')
nwm = c('calJac3','saiBol1','cebCap1','aotNan1')
keeps = Reduce(union, list(apes, owm, nwm))

tree = keep.tip(tree, keeps)
tree = unroot(tree)

# Load peaks
pksAll = read.table('peaks_sorted.bed')

# Loop through chromosomes and find substitutions in all lineages
dir.create('Ancestral_Substitutions')
chrs = c(paste0('chr', 1:22), 'chrX', 'chrY')
for(i in 1:length(chrs)){

	# Filter peaks for the given chromosome
	chr = chrs[i]
	pks = pksAll[pksAll[,1] == chr,]
	pks$id = paste0(pks[,1], ':', pks[,2], '-', pks[,3])

	# Load extracted alignments
	subalns = readRDS(paste0('Extracted_CRE_Alignments/', chr, '.RDS'))
	subalns = unlist(subalns, recursive = F)

	print(paste0('# of CREs: ', length(subalns)))

	# Find substitutions
	resL = mclapply(1:length(subalns), mc.cores = 30, function(x){

		# Peak id
		pkid = as.character(pks[x, 'id'])
		
		# Extract peak regions from the alignment
		pkAln = subalns[[x]]

		# Find ancestral sequences
		alStr = ancestralize(pkAln)

		# Find node specific substitutions
		subs = findSubs(alStr = alStr, pks = pks, pos = x)

		if(x%%100==0){print(x)}
		subs
	})


	ancSubsDF = do.call(rbind, resL)

	# Remove the substitutions that took place more than once
	qwL = list()
	for(k in 1:nrow(ancSubsDF)){


		tmp = make.unique(as.character(ancSubsDF[k,5:10]))
		tmp2 = tmp[which(!(grepl('\\..', tmp)))[2]:length(tmp)]

		tmp3 = tmp2[grepl('\\..', tmp2)]
		tmp4 = gsub('\\..', '', tmp3)

		qwL[[k]] = length(unique(tmp4))
		if(k%%1000==0){print(k)}
	}

	qw = unlist(qwL)
	ancSubsDF2 = ancSubsDF[qw != 2, ]

	saveRDS(ancSubsDF2, paste0('Ancestral_Substitutions/', chr, '_AncestralSubstitutions.RDS'))

}





