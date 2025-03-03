require(plyr)
require(dplyr)
require(tidyverse)
require(tidyr)
require(ggplot2)
require(ggpubr)
require(reshape2)
require(data.table)
require(dplyr)
require(tidyverse)
require(rio)
require(GenomicRanges)
require(data.table)
require(rphast)
require(ape)
require(dplyr)
require(parallel)
require(Biostrings)
require(ggpubr)
require(seqinr)
require(phangorn)
require(msa)
require(TFBSTools)
require(JASPAR2020)
require(motifmatchr)
require(parallel)
source('SCRIPTS/Functions.R')

args = commandArgs(trailingOnly = TRUE)

for (arg in args) {
	split_arg <- strsplit(arg, "=")[[1]]
	var_name <- split_arg[1]
	var_value <- split_arg[2]
	print(var_name)
	print(var_value)
	assign(var_name, var_value)
}

nthreads = as.integer(nthreads)
dir.create(outdir, , showWarnings = F)

# Enumerate the nodes
mapnames = setNames(c(1:6), c("Human", "HC", "HCGo", "Great_Ape", "Ape", "Anthropoids"))

# Load the motif data bases
opts = list()
opts[["species"]] <- "9606"
opts[['collection']] = 'CORE'
PFMatrixList <- getMatrixSet(JASPAR2020, opts)
PWMatrixList <- toPWM(PFMatrixList, pseudocounts=0.8)

####
## Prepare tree
####

# Load the tree
treetop = read.table(treetop_fn)
treetop = treetop$V1

# Change the names from long version to short version to match with multiz alignment
name_df = read.table("IN_DATA/short_long_species_names.txt", header = T)
for(i in 1:nrow(name_df)){
	treetop = gsub(name_df[i, 'long_name'], name_df[i, 'short_name'], treetop)
}

tree = read.tree(text = treetop)
tree = unroot(tree)

####
## Prepare alignment
####

# Load the alignment
align = strip.gaps.msa(read.msa(paste0("IN_DATA/", chr, ".maf"), pointer.only=TRUE))

# Keep only given species
align2 = align[tree$tip,]


####
## FIND MOTIF OCCURRENCES PER LINEAGE
####

# Load peaks
pks = read.table(pks_fn)
pks = pks[pks[,1] == chr,]
pks$id = paste0(pks[,1], ':', pks[,2], '-', pks[,3])


# Read extracted alignments
subalns = readRDS(extracted_align)
subalns = unlist(subalns, recursive = F)
print(paste0('Number of peaks: ', length(subalns)))

tmpL = mclapply(1:length(subalns), mc.cores = nthreads, function(x){

	# Extract the sequences and find the ancestral sequences
	pkAln = subalns[[x]]
	alStr = ancestralize(pkAln, submodel = submodel, phyltree = tree,
				ancestral_nodes = c('35', '34', '33', '32', '21'))

	# Find motif occurrences in each node per CRE
	motif_ix = matchMotifs(subject = alStr, pwms = PWMatrixList, bg = 'even')
	matchOut = motifMatches(motif_ix)
	if(x%%50 == 0){print(x)}
	matchOut
})

saveRDS(tmpL, paste0(outdir, '/', chr, '_motifMatrix.RDS'))


####
## FIND TFBS EVOLUTION PER LINEAGE
####

motresL = list()
for(j in 1:length(tmpL)){

	matchOut = tmpL[[j]]
	if(grepl('Error', matchOut)){	
		# Extract the sequences and find the ancestral sequences
		pkAln = subalns[[x]]
		alStr = ancestralize(pkAln, submodel = submodel, phyltree = tree,
					ancestral_nodes = c('35', '34', '33', '32', '21'))

		# Find motif occurrences in each node per CRE
		motif_ix = matchMotifs(subject = alStr, pwms = PWMatrixList, bg = 'even')
		matchOut = motifMatches(motif_ix)
	}

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

		# Assign gains: Gains start with the first node
		if(difs[1] == 1){
			node = difs[length(difs)]
			dfL[[i]] = data.frame(CRE = pkid, motif = colnames(difmat2)[i], type = 'Gain',
							node = names(mapnames[node]))

		# Assign losses: Losses end at the node when the motif occurrence was lost
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

saveRDS(motresDF, paste0(outdir, '/', chr, '_motifEvo.RDS'))



