require(rphast)
require(ape)
require(dplyr)
require(parallel)
require(Biostrings)
require(ggpubr)
require(seqinr)
require(phangorn)
require(msa)
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

nthreads = as.numeric(nthreads)


####
## Extract the peaks from the alignment
####

# Load the alignment
align = strip.gaps.msa(read.msa(msa_fn, pointer.only=TRUE))

# Keep only given species
keepspecies = read.table(keepspecies_fn)
keepspecies = keepspecies$V1
align2 = align[keepspecies,]

# Load peaks
pks = read.table(peak_fn)
pks = pks[pks[,1] == chr,]
pks$id = paste0(pks[,1], ':', pks[,2], '-', pks[,3])

# Extract from the alignment
ofs = offset.msa(align)
tmp1 = seq(1,nrow(pks),nthreads)
tmp2 = c(tmp1, nrow(pks))
tmp2 = tmp2[-1]
tmp2[1:length(tmp2)-1] = tmp2[1:length(tmp2)-1] - 1

subalnsL = list()
for(i in 1:length(tmp1)){

	subalnsL[[i]] = mclapply(tmp1[i]:tmp2[i], mc.cores = nthreads, function(x){

		# Extract peak regions from the alignment
		pkAln = sub.msa(align2, start.col = pks[x,2] - ofs, end.col = pks[x,3] - ofs, pointer.only = F)
		pkAln
	})
	
	print(i)
}


dir.create(outdir, showWarnings = F)
saveRDS(subalnsL, paste0(outdir, '/', chr, '.RDS'))

