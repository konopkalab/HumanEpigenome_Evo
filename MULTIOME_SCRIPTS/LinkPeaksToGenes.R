rm(list = ls())
library(Signac)
library(Seurat)
library(Matrix)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(readr)

args = commandArgs(trailingOnly = TRUE)
fn_seurat = args[1]
fn_gns = args[2]
outdir = args[3]
index = args[4] %>% as.numeric

print(index)

# Read data
brainSeurat = read_rds(fn_seurat)
gnsList = read_rds(fn_gns)

print(gnsList[[index]])

# Peaks correlated with genes
tmpseurat = LinkPeaks(object = brainSeurat, peak.assay = "ATAC",
		expression.assay = "SCT", genes.use = gnsList[[index]], score_cutoff = 0.01)
randgns_links = tmpseurat@assays$ATAC@links %>% as.data.frame
write_rds(randgns_links, paste0(outdir, 'Trevino_multiome_randgns_links_', index, '.RDS'))

