rm(list = ls())
library(Signac)
library(Seurat)
library(Matrix)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(readr)
library(tgutil)
library(tidyverse)
library(ggpubr)
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")



# Read adult dataset. Only fetal brain links will be displayed (no fragments)
combinedSeurat = read_rds('MULTIOME_ADULT/PKS_TO_GNS/seurat_object_mypeaks_with_fragments.RDS')
DefaultAssay(combinedSeurat) = 'ATAC'

# Read links
all_links = readRDS('MULTIOME_FETAL/PKS_TO_GNS/Trevino_LinkPeaksToGenes_FINAL.RDS')
all_div_links = readRDS('MULTIOME_FETAL/Linked_Gene_Final_Results_Sign_FDR.RDS')

# Load divergent CREs
hvarsAll = read_rds('02_GroupCREs/FETAL/FetalCREs_SignTested_SD_1_FC_ConsAdded.RDS')
human_div = hvarsAll[hvarsAll$HumanSign == 'Sign', 'CRE'] %>% sub(':', '-', .)
hc_div = hvarsAll[hvarsAll$HCSign == 'Sign', 'CRE'] %>% sub(':', '-', .)

# Plot a human divergent gene
all_links$is_div = ifelse(all_links$peak %in% human_div, 'Human_Divergent', 'NS')
all_links$is_div = factor(all_links$is_div)
all_links$color = 'blue'
top_peaks = all_links[all_links$gene == 'TFEB' & all_links$pvalue < 0.05, ] %>% head(30) %>% tail(20)
highlight = top_peaks[top_peaks$is_div == 'Human_Divergent', 'peak'] %>% StringToGRanges
highlight$color = "darkblue"

links_gr = makeGRangesFromDataFrame(top_peaks, keep.extra.columns=T)
Links(combinedSeurat) = links_gr

interv = "chr6-41306203-41800000"

pdf('TFEB_FetalHumanDivergent.pdf')
CoveragePlot(
  object = combinedSeurat,
  region = interv,
  region.highlight = highlight,
  window = 2000,
  extend.downstream = 0,
  extend.upstream = 10000
)
dev.off()

# Plot substitutions per CRE
allcountsLong = readRDS('02_GroupCREs/FETAL/FetalSubstitutions_Counted_BranchNormalized_LongFormat.RDS')
colors = c('blue', 'lightblue', 'purple', 'pink', 'red')

# Load divergent CREs
hvarsAll = read_rds('02_GroupCREs/FETAL/FetalCREs_SignTested_SD_1_FC_ConsAdded.RDS')
human_div = hvarsAll[hvarsAll$HumanSign == 'Sign', 'CRE'] %>% sub(':', '-', .)

all_links = readRDS('MULTIOME_FETAL/PKS_TO_GNS/Trevino_LinkPeaksToGenes_FINAL.RDS')
all_links$is_div = ifelse(all_links$peak %in% human_div, 'Human_Divergent', 'NS')
all_links$is_div = factor(all_links$is_div)

top_peaks = all_links[all_links$gene == 'TFEB' & all_links$pvalue < 0.05, ] %>%  .[,'peak']
#plotpeaks = top_peaks[top_peaks$is_div == 'Human_Divergent', 'peak'] %>% sub('-', ':', .)
plotpeaks = top_peaks %>% sub('-', ':', .)

for(i in 1:length(plotpeaks)){

	toplot = reshape2::melt(allcountsLong[plotpeaks[i],])

	pdf(paste0('SUBSTITUTIONS_HUMANDIV_TFEB', '_', 'Peak', i, '.pdf'))
	print( ggbarplot(toplot, x = 'variable', y = 'value', color = 'variable', fill = 'variable', palette = colors, ylim = c(0,3)) +
	ylab('Number of substitutions per MY per KB') + xlab('') +
	theme(text = element_text(size=20)) +
	theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
	rotate_x_text(45) +
	NoLegend() )
	dev.off()
}












