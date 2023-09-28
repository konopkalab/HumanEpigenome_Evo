rm(list = ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(harmony)
library(Signac)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
require(gprofiler2)

####
## PEAK CONVERSION (Multiome to scATAC)
## This is needed for compatibility with the scATAC-seq results.
## Re-calling the peaks would be better but raw data is not publicly available for Trevino et al. 2021
####

# Multiome peaks
authPeaks = fread('MULTIOME_FETAL/GSE162170_multiome_atac_consensus_peaks.txt.gz')
authPeaks = as.data.frame(authPeaks) %>% makeGRangesFromDataFrame

# scATAC-seq peaks
scatacPeaks = read.table('ATACSEQ_FETAL_CORTEX/fetalBed_Access.bed')
colnames(scatacPeaks) = c('seqnames', 'start', 'end')
scatacPeaks = as.data.frame(scatacPeaks) %>% makeGRangesFromDataFrame

# Determine the minimum overlap to require based on a unique overlap of all peaks to the scatac peaks
ovs = seq(100, 500, 10)
divL = list()
for(i in 1:length(overlaps)){
	tmp = findOverlaps(query=authPeaks, subject=scatacPeaks, minoverlap = ovs[i], select = 'all')
	divL[[i]] = length(unique(tmp@to)) / length(tmp@to)
}
minov = ovs[min(which(unlist(divL) == 1))]

# Get overlapping indices between the two peak sets
ovinds = findOverlaps(query=authPeaks, subject=scatacPeaks, minoverlap = minov, select = 'all')

# Get matching peaks from both
authPeaks_ovs = authPeaks[ovinds@from]
scatacPeaks_ovs = scatacPeaks[ovinds@to]

####
## ATAC-SEQ MATRIX AND METADATA
####

# Load the peak-cell matrix
atacMat = fread('MULTIOME_FETAL/GSE162170_multiome_atac_counts.tsv.gz')

# Convert to matrix format
atacMat = as.data.frame(atacMat) %>% as.matrix() %>% as(., "sparseMatrix")
gc()

# Subset ATAC-seq to the peaks that overlap with scATAC-seq peak list
atacMat= atacMat[ovinds@from,]

# Add scATAC-seq peak peak names
rownames(atacMat) = GRangesToString(scatacPeaks_ovs, sep = c('-', '-'))

# Meta data
atacMeta = fread('MULTIOME_FETAL/GSE162170_multiome_cell_metadata.txt.gz') %>% as.data.frame
atacCB = atacMeta$Cell.ID
rownames(atacMeta) = atacMeta$Cell.ID
colnames(atacMat) = atacCB

####
## RNA-SEQ MATRIX AND SEURAT OBJECT
####

# Load the peak-cell matrix
rnaMat = fread('MULTIOME_FETAL/GSE162170_multiome_rna_counts.tsv.gz')
rnaMat = as.data.frame(rnaMat)
rownames(rnaMat) = rnaMat[,1]
rnaMat[,1] = NULL

# Convert to matrix format
rnaMat = as.matrix(rnaMat) %>% as(., "sparseMatrix")
gc()

# Ensembl id to gene symbols
symbols = ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(rnaMat), keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# Convert to symbol
rownames(rnaMat) = symbols[match(rownames(rnaMat), symbols$GENEID), 'SYMBOL']
rnaMat = rnaMat[!(is.na(rownames(rnaMat))),]

# Create seurat object on RNA-seq
seurObj = CreateSeuratObject(counts = rnaMat, meta.data = atacMeta, min.cells = 0, min.features = 0)

# get gene annotations for hg38
annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) = paste0('chr', seqlevels(annotation))
write_rds(annotation, 'signac_ensdb_hsapiens_v86_annotation.RDS')

# Add ATAC-seq assay
seurObj[["ATAC"]] = CreateChromatinAssay(counts = atacMat, sep = c("-", "-"), annotation = annotation)
write_rds(seurObj, 'MULTIOME_FETAL/seurObj.RDS')

# Compute the GC content for each peak
DefaultAssay(seurObj) = "ATAC"
seurObj = RegionStats(seurObj, genome = BSgenome.Hsapiens.UCSC.hg38)

# Normalize gene expression
DefaultAssay(seurObj) = "RNA"
seurObj = SCTransform(seurObj, ncells = 1000)
write_rds(seurObj, 'MULTIOME_FETAL/seurObj.RDS')

# Get protein coding genes
annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) = paste0('chr', seqlevels(annotation))
annot_pr = annotation[annotation$gene_biotype == 'protein_coding',]
pr_gns = annot_pr$gene_name %>% unique

expgns = rnaMat[rowSums(rnaMat) > 0,] %>% rownames
pr_gns_exp = expgns[expgns %in% pr_gns]
gnsList = split(pr_gns_exp, ceiling(seq_along(pr_gns_exp)/100))

write_rds(gnsList, 'MULTIOME_FETAL/pr_gns.RDS')

####
## PEAKS TO GENES
####

# Run links to peaks script (easier to parallelize)
# MULTIOME_FETAL/01_PeakToGene.sh

# Merge peak-gene links
resdf = list.files('MULTIOME_FETAL/PKS_TO_GNS', pattern = "Tre.*RDS", full.names = T) %>% lapply(., function(x){readRDS(x)}) %>% do.call(rbind, .)

write_rds(resdf, 'MULTIOME_FETAL/PKS_TO_GNS/Trevino_LinkPeaksToGenes_FINAL.RDS')






