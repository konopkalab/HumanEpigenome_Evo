require(rphast)
require(ape)
require(dplyr)
require(parallel)
require(Biostrings)
require(ggpubr)
require(Seurat)
require(reshape2)
require(pheatmap)
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

bcg_size = as.numeric(bcg_size)
plt_pref = paste0('PLOTS/', plt_pref)

# Load evolutionary groups
evoGroups = readRDS(evo_group_fn)
human_cres = evoGroups[evoGroups$HumanSign == 'Sign', 'CRE']
hc_cres = evoGroups[evoGroups$HCSign == 'Sign', 'CRE']
hcgo_cres = evoGroups[evoGroups$HCGoSign == 'Sign', 'CRE']
gape_cres = evoGroups[evoGroups$Great_ApeSign == 'Sign', 'CRE']
ape_cres = evoGroups[evoGroups$ApeSign == 'Sign', 'CRE']
con_cres = evoGroups[evoGroups$consSign == 'Sign', 'CRE']

evoL = list(human_cres, hc_cres, hcgo_cres, gape_cres, ape_cres, con_cres)
names(evoL) = c('Human', 'HC', 'HCGo', 'Great_Ape', 'Ape', 'Cons')

# Load cell type markers
ctmarks = readRDS(ctmarks_fn)
ctmarks = ctmarks[ctmarks$logFC < -1,]
ctmarksL = split(ctmarks, ctmarks$CellType)
ctmarksL = lapply(ctmarksL, function(x){x$Gene})

## Run enrichment ##

# Fisher's exact test enrichment and depletion
toplot = geneOvEnrDep(gnL1 = evoL, gnL2 = ctmarksL, bcg_size, plot = F, hg = 10, wd = 14)
toplot$sign_label = formatC(toplot$FDR, format = "e", digits = 0)

pdf(paste0(plt_pref, '_EvoGroup_CellTypeEnrichment.pdf'), width = 10, height = 6)
print( ggscatter(toplot, x = 'var1', y = 'var2', color = 'OddsRatio', size = 'log10FDR') +
geom_label(data = toplot, aes(label = sign_label), color="black",
	label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
labs(x="", y="") +
scale_size_continuous(range = c(5,25)) +
scale_color_gradient2(midpoint = 1, low = 'blue', mid = 'white', high = 'red') +
theme_classic() +
theme(text = element_text(size=20, face = 'bold')) +
geom_text(aes(label = sign_label_2), vjust = 1.2, colour = "darkgreen", fontface = 'bold', size = 20 ) +
rotate_x_text(45) )
dev.off()


