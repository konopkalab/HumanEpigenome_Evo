rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(reshape2)
library(data.table)
options(scipen = 999)

####
## EXPAND PEAKS AND WRITE AS BED FILE
####

# Load grouped CREs
cres = readRDS('AdultCREs_SignTested_SD_1_FC_ConsAdded.RDS')

# Top N CREs
n = 20000
type = 'TOP20K_EXPAND_50kb'
dir.create(type)

# For conserved CREs
cres$consFC = cres[, c('hFC', 'hcFC', 'hcgoFC', 'gapeFC', 'apeFC')] %>% rowMeans
cres$consFC = 1/cres$consFC

# Loop through the nodes, create bed file and write bed file for LDSC
nodes = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'cons')
fcs = c('hFC', 'hcFC', 'hcgoFC', 'gapeFC', 'apeFC', 'consFC')
nodes2 = paste0(nodes, 'Sign')
for(i in 1:length(nodes2)){

	# Get CREs
	node_cres = cres[cres[, nodes2[i]] == 'Sign', ]
	node_cres = node_cres[order(node_cres[[fcs[i]]], decreasing = T)[1:n], 'CRE']

	# Make bed
	tmp_bed = node_cres %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
			do.call(rbind, .) %>% as.data.frame() %>%
			mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3)))

	tmp_bed$V2 = tmp_bed$V2 - 25000
	tmp_bed$V3 = tmp_bed$V3 + 25000
	tmp_bed = tmp_bed[tmp_bed$V2 > 0,]

	# Write bed
	dir.create(paste0(type, '/ADULT_', nodes[i]))
	write.table(tmp_bed, col.names = F, row.names = F, sep = '\t', quote = F, paste0(type, '/ADULT_', nodes[i], '/peakset.bed'))
}





