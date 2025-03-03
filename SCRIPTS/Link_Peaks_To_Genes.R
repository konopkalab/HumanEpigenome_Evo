require(Signac)
require(Seurat)
require(Matrix)
require(EnsDb.Hsapiens.v86)
require(BSgenome.Hsapiens.UCSC.hg38)
require(dplyr)
require(readr)

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

index = as.numeric(index)
dir.create(outdir, showWarnings = F)

print(index)

# Read data
brainSeurat = read_rds(fn_seurat)
gnsList = read_rds(fn_gns)

print(gnsList[[index]])

# Peaks correlated with genes
tmpseurat = LinkPeaks(object = brainSeurat, peak.assay = "ATAC",
		expression.assay = "SCT", genes.use = gnsList[[index]], score_cutoff = 0.01)
randgns_links = tmpseurat@assays$ATAC@links %>% as.data.frame
write_rds(randgns_links, paste0(outdir, '/Multiome_randgns_links_', index, '.RDS'))
