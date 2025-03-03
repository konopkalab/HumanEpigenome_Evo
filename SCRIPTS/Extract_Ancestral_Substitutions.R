require(rphast)
require(ape)
require(dplyr)
require(parallel)
require(Biostrings)
require(ggpubr)
require(seqinr)
require(phangorn)
require(msa)
require(readr)
require(VennDiagram)
source('SCRIPTS/Functions.R')


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


nthreads = as.numeric(nthreads)
pksAll = read.table(pks_fn)
nodes = c('35', '34', '33', '32', '21')
dir.create(outdir, showWarnings = F)

####
## Prepare tree
####

# Load the tree
treetop = read.table(treetop_fn)
treetop = treetop$V1

# Change the names from long version to short version to match with multiz alignment
name_df = read.table("DATA/short_long_species_names.txt", header = T)
for(i in 1:nrow(name_df)){
	treetop = gsub(name_df[i, 'long_name'], name_df[i, 'short_name'], treetop)
}

tree = read.tree(text = treetop)
tree = unroot(tree)

####
## Find substitutions between ancestral states
####

fls = list.files('DATA/Extracted_CRE_Alignments', pattern = 'chr.*RDS', full.names = T)
chrs = gsub('.*/|.RDS', '', fls)
for(i in 1:length(chrs)){

	# Filter peaks
	chr = chrs[i]
	pks = pksAll[pksAll[,1] == chr,]
	pks$id = paste0(pks[,1], ':', pks[,2], '-', pks[,3])

	# Read extracted alignments
	subalns = readRDS(fls[i])
	subalns = unlist(subalns, recursive = F)

	print(paste0('# of CREs: ', length(subalns)))

	# Find substitutions
	resL = mclapply(1:length(subalns), mc.cores = nthreads, function(x){

		# Peak id
		pkid = as.character(pks[x, 'id'])

		# Extract peak regions from the alignment
		pkAln = subalns[[x]]

		alStr = ancestralize(pkAln, submodel = submodel, phyltree = tree,
					ancestral_nodes = nodes)

		# Extract node specific changes
		subs = findSubs(alStr = alStr, pks = pks, pos = x)

		if(x%%100==0){print(x)}
		subs
	})


	ancSubsDF = do.call(rbind, resL)

	# Remove the entries with error. These come from the peaks that do not have multi-species alignment data or even human genome sequence data for hg38 (likely and error during peak calling)
	ancSubsDF = ancSubsDF[!(grepl('Error', ancSubsDF$CRE)),]

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

	saveRDS(ancSubsDF2, paste0(outdir, '/', chr, '_AncestralSubstitutions.RDS'))

}

# Read and combine all substitutions
allcomb = list.files(path = outdir, pattern = 'AncestralSubstitutions.RDS', full.names = T) %>% lapply(., function(x){readRDS(x)}) %>% do.call(rbind, .)

# Save for easier loading
write_rds(allcomb, paste0(outdir, '/AncestralSubstitutions_ALL.RDS'))

