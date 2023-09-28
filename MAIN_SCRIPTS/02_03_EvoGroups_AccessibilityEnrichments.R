rm(list = ls())
library(rphast)
library(dplyr)
library(parallel)
library(ggpubr)
library(Seurat)
library(reshape2)
library(pheatmap)
library(GeneOverlap)
source('Functions.R')

####
## ENRICHMENT -- GLM, LENGTH COVARIATE
####

# Load data
adultSub = readRDS('AdultCREs_SignTested_SD_1_FC_ConsAdded.RDS')

darsAll = readRDS('DATASETS/Adult_SpeciesSpecific_DARs_snATACseq.RDS')
darsH = darsAll[darsAll$Evolution == 'Human_Specific', 'Gene'] %>% gsub(':|-', '_', .) %>% unique
darsC = darsAll[darsAll$Evolution == 'Chimp_Specific', 'Gene'] %>% gsub(':|-', '_', .) %>% unique
darsM = darsAll[darsAll$Evolution == 'MvsHC', 'Gene'] %>% gsub(':|-', '_', .) %>% unique
tmp = unique(darsAll$Gene) %>% gsub(':|-', '_', .)
darsCons = setdiff(tmp, unique(c(darsH, darsC, darsM))) %>% gsub(':|-', '_', .) %>% unique

allcres = gsub(':|-', '_', darsAll$Gene) %>% unique

humdivCRE = adultSub[adultSub$HumanSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique %>% intersect(., allcres)
hcdivCRE = adultSub[adultSub$HCSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique %>% intersect(., allcres)
hcgodivCRE = adultSub[adultSub$HCGoSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique %>% intersect(., allcres)
gapedivCRE = adultSub[adultSub$Great_ApeSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique %>% intersect(., allcres)
apedivCRE = adultSub[adultSub$ApeSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique %>% intersect(., allcres)
consCRE = adultSub[adultSub$consSign == 'Sign', 'CRE'] %>% gsub(':|-', '_', .) %>% unique %>% intersect(., allcres)

# Changed specifically in one species
darsH = setdiff(darsH, union(darsC, darsM))

# Create the data frame
allcres = gsub(':|-', '_', darsAll$Gene) %>% unique %>% data.frame
allcres$hAcc = ifelse(allcres[,1] %in% darsH, 1, 0)
allcres$consAcc = ifelse(allcres[,1] %in% darsCons, 1, 0)

allcres$hSub = ifelse(allcres[,1] %in% humdivCRE, 1, 0)
allcres$hcSub = ifelse(allcres[,1] %in% hcdivCRE, 1, 0)
allcres$hcgoSub = ifelse(allcres[,1] %in% hcgodivCRE, 1, 0)
allcres$gapeSub = ifelse(allcres[,1] %in% gapedivCRE, 1, 0)
allcres$apeSub = ifelse(allcres[,1] %in% apedivCRE, 1, 0)
allcres$consSub = ifelse(allcres[,1] %in% consCRE, 1, 0)

# Calculate the length
bed = allcres[,1] %>% strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3)))
allcres$length = bed[,3]-bed[,2]


# Run Enrichment
accs = colnames(allcres)[grepl('Acc', colnames(allcres))]
subs = colnames(allcres)[grepl('Sub', colnames(allcres))]
dfL2 = list()
for(j in 1:length(accs)){
	dfL = list()
	for(i in 1:length(subs)){

		tmp = glm(allcres[,accs[j]]~length + allcres[,subs[i]], allcres, family = 'binomial')
		fc = coef(summary(tmp))[3,1]
		pval = coef(summary(tmp))[3,4]
		dfL[[i]] = data.frame(Accessibility = accs[j], Substitution = subs[i], fc = fc, pval = pval)
	}
	dfL2[[j]] = do.call(rbind, dfL)
	print(j)
}


resdf = do.call(rbind, dfL2)
resdf$FDR = p.adjust(resdf$pval, method = 'fdr')
resdf$log10_round = round(-log10(resdf$FDR), digits = 2)
resdf$fc_log = round(resdf$fc, digits = 2)
resdf$fdr_Sc = formatC(resdf$FDR, format = "e", digits = 2)

resdf$Accessibility = factor(resdf$Accessibility, levels = c('hAcc', 'cAcc', 'mAcc', 'consAcc'))
resdf$Substitution = factor(resdf$Substitution, levels = c('hSub', 'hcSub', 'hcgoSub', 'gapeSub', 'apeSub', 'consSub'))
resdf$is_sign = ifelse(resdf$FDR < 0.05 & abs(resdf$fc) > 0, '*', '')
resdf$absFC = round(exp(resdf$fc_log), digits = 2)

pdf('Substitution_vs_Accessibility_Enrichments.pdf', height = 6, width = 9)
ggscatter(resdf, x = 'Substitution', y = 'Accessibility', color = 'absFC', size = 'log10_round') +
  geom_text(aes(label = absFC), vjust = 0.5, colour = "black", size = 5 ) +
  labs(x="Substitution", y="Accessibility") + scale_size_continuous(range = c(2,20)) +
  scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(45)
dev.off()






