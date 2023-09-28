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
library(parallel)
source('Functions.R')

# Set chromosome
chr = 'chr22'

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

# Load the alignment
align = strip.gaps.msa(read.msa(paste0("UCSC_MULTIZ30_2017/", chr, ".maf"), pointer.only=TRUE))

# Keep only given species
align2 = align[keeps,]

# Load peaks
pks = read.table('peaks_sorted.bed')
pks = pks[pks[,1] == chr,]
pks$id = paste0(pks[,1], ':', pks[,2], '-', pks[,3])

# Offset of the alignment
ofs = offset.msa(align)

# Extract from the alignment
tmp1 = seq(1,nrow(pks),23)
tmp2 = c(tmp1, nrow(pks))
tmp2 = tmp2[-1]
tmp2[1:length(tmp2)-1] = tmp2[1:length(tmp2)-1] - 1

subalnsL = list()
for(i in 1:length(tmp1)){

subalnsL[[i]] = mclapply(tmp1[i]:tmp2[i], mc.cores = 23, function(x){

	# Extract peak regions from the alignment
	pkAln = sub.msa(align2, start.col = pks[x,2] - ofs, end.col = pks[x,3] - ofs, pointer.only = F)
	pkAln})
	
	print(i)
}

dir.create('Extracted_CRE_Alignments')
saveRDS(subalnsL, paste0('Extracted_CRE_Alignments/', chr, '.RDS'))












