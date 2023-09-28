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
library(rphast)
library(ape)
library(dplyr)
library(parallel)
library(Biostrings)
library(ggpubr)
library(seqinr)
library(phangorn)
library(msa)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(parallel)
source('Functions.R')

# Set chromosome
chr = 'chr22'

# Enumerate the nodes
mapnames = setNames(c(1:6), c("Human", "HC", "HCGo", "Great_Ape", "Ape", "Anthropoids"))

# Load the motif data bases
opts = list()
opts[["species"]] = "9606"
opts[['collection']] = 'CORE'
PFMatrixList = getMatrixSet(JASPAR2020, opts)
PWMatrixList = toPWM(PFMatrixList, pseudocounts=0.8)


# Load the alignment #

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

# Load the alignment
align = strip.gaps.msa(read.msa(paste0(UCSC_MULTIZ30_2017/", chr, ".maf"), pointer.only=TRUE))

# Keep only given species
align2 = align[keeps,]

####
## FIND DIFFERENTIAL MOTIF OCCURRENCES
####

# Load peaks
pks = read.table('peaks_sorted.bed')
pks = pks[pks[,1] == chr,]
pks$id = paste0(pks[,1], ':', pks[,2], '-', pks[,3])

# Load extracted alignments
subalns = readRDS(paste0('Extracted_CRE_Alignments/', chr, '.RDS'))
subalns = unlist(subalns, recursive = F)
print(length(subalns))

# Loop through extracted alignments to find all motif occurrences
tmpL = mclapply(1:length(subalns), mc.cores = 23, function(x){

	# Extract the sequences and find the ancestral sequences
	pkAln = subalns[[x]]
	alStr = ancestralize(pkAln)
	#pkid = pks[x, 'id']

	# Find motif occurrences in each node per CRE
	motif_ix = matchMotifs(subject = alStr, pwms = PWMatrixList, bg = 'even')
	matchOut = motifMatches(motif_ix)
	if(x%%50 == 0){print(x)}
	matchOut
})

saveRDS(tmpL, paste0('Motif_Evolution/', chr, '_motifMatrix.RDS'))

# Find motif gains and losses
motresL = list()
for(j in 1:length(tmpL)){

	matchOut = tmpL[[j]]
	pkid = pks[j, 'id']

	# Keep only the positions with altered TFBS site
	motsMATCHR = names(which(colSums(matchOut) > 0 & colSums(matchOut) < 6))
	difmat = matchOut[, motsMATCHR, drop=FALSE]

	if(length(difmat) == 0){motresL[[j]] = NULL;next}

	# Remove the positions with more than one change
	dfL1 = list()
	for(i in 1:ncol(difmat)){

		difs = difmat[,i]

		# Subtract the next element from the previous
		# There are 6 elements, 5 potential changes.
		# Four of these should be no change (i.e 0).
		tmp = table(diff(difs))
		if(sum(names(tmp) == '0') == 0){next} # Remove the ones without no change
		if(tmp[names(tmp) == '0'] == 4){dfL1[[i]] = i}
	}

	# If all positions changed more than once, skip it
	if(length(dfL1) == 0){
		motresL[[j]] = NULL;next
	} else{
		difmat2 = difmat[,unlist(dfL1), drop=FALSE]
	}

	# Create a data frame of gained / lost TF sites
	dfL = list()
	for(i in 1:ncol(difmat2)){

		difs = which(difmat2[,i] == T)

		# Gains start with the first node
		if(difs[1] == 1){
			node = difs[length(difs)]
			dfL[[i]] = data.frame(CRE = pkid, motif = colnames(difmat2)[i], type = 'Gain',
							node = names(mapnames[node]))

		# Losses end at the node when the motif occurrence was lost
		} else{
			node = difs[1] - 1
			dfL[[i]] = data.frame(CRE = pkid, motif = colnames(difmat2)[i], type = 'Loss',
							node = names(mapnames[node])) 
		}
	}

	motresL[[j]] = do.call(rbind, dfL)
	if(j%%100 == 0){print(j)}
}

motresDF = do.call(rbind, motresL)

# Add TF names
mots = sapply(PFMatrixList, function(x){x@name})
motresDF$tf = mots[motresDF$motif]

saveRDS(motresDF, paste0('Motif_Evolution/', chr, '_motifEvo.RDS'))






