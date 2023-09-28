
findSubs = function(alStr, pks, pos){

	# Convert to data frame. One position per column
	alStrDF = as.data.frame(alStr)
	alStrDF2 = strsplit(alStrDF[,1], "") %>% do.call(rbind, .) %>% as.data.frame
	rownames(alStrDF2) = rownames(alStrDF)

	# Remove the positions with N (ancestral sequence could not reconstructed with certainty)
	unkAnc = apply(alStrDF2, 2, function(x){any(x == 'N')}) %>% which %>% names
	alStrDF2 = alStrDF2[, !(colnames(alStrDF2) %in% unkAnc)]

	alStrDF3 = apply(alStrDF2, 2, function(x){table(x) %>% length})

	# If all positions are conserved (which is rare), skip it
	if(length(table(alStrDF3)) == 1){return(NULL)}
	else if (sum(names(table(alStrDF3)) == '2') == 0){return(NULL)}
	else{

	# Boolean matrix of substitutions
	alStrDFBool = apply(alStrDF2, 2, function(x){duplicated(x)})
	alStrDFBool[1,] = T
	alStrDFBool = !alStrDFBool

	# Positions of the changes within the CRE
	changePos = apply(alStrDFBool, 2, function(x){which(x)}) %>%
			sapply(., function(x){length(x) == 1}) %>%
			which %>% names %>% gsub('V', '', .) %>% as.numeric


	speciesPos = apply(alStrDFBool, 2, function(x){which(x)}) %>%
			.[ lapply(.,length) == 1] %>%
				unlist %>% as.numeric

	speciesPos = speciesPos - 1
	speciesChange = rownames(alStrDF2)[speciesPos]

	# Create the data frame of changes
	cre = paste0(chr, ':', pks[pos,2], '-', pks[pos,3])
	df = data.frame(CRE = cre, pos = changePos, ucscPos = pks[pos,2] + changePos - 1, nodes = speciesChange)

	# Add the nucleotides
	changePos_char = paste0('V', changePos)
	tmpdf = alStrDF2[, changePos_char] %>% t %>% as.data.frame
	colnames(tmpdf) = rownames(alStrDF2)
	res = cbind(df, tmpdf)

	}

	return(res)
}



ancestralize = function(pkAln, tmpdir = '/home2/s422159/workdir/TMP/'){

	require(stringi)
	pkid = stri_rand_strings(1, 10, '[a-z]')
	write.msa(pkAln, file = paste0(tmpdir, pkid, '.phy'), format = 'PHYLIP')
	write.msa(pkAln, file = paste0(tmpdir, pkid, '.fa'), format = 'FASTA')

	# Make unknowns gaps
	togap(paste0(tmpdir, pkid, '.phy'))
	togap(paste0(tmpdir, pkid, '.fa'))

	# Re-read the alignment
	# Different site patterns: substitution patterns, not total number of substitutions
	primates = read.phyDat(paste0(tmpdir, pkid, '.phy'))

	# Construct ancestral sequences
	fit = pml(tree, primates)
	fit = optim.pml(fit, model="F81", control = pml.control(trace=0))
	ancMl = ancestral.pml(fit, "ml")

	# Get human and ancestral nodes
	human = giveAnc(ancMl, node = 'hg38')
	#chimp = giveAnc(ancMl, node = 'panTro5')
	#gorilla = giveAnc(ancMl, node = 'gorGor5')
	#orangutan = giveAnc(ancMl, node = 'ponAbe2')
	#gibbon = giveAnc(ancMl, node = 'nomLeu3')

	hcAnc = giveAnc(ancMl, node = '25')
	hcgoAnc = giveAnc(ancMl, node = '24')
	gapeAnc = giveAnc(ancMl, node = '23')
	apeAnc = giveAnc(ancMl, node = '22')
	allAnc = giveAnc(ancMl, node = '21')

	# Write FASTA
	write.fasta(human, 'Human', paste0(tmpdir, pkid, '_ancestors.fa'), open = "w", nbchar = 70, as.string = FALSE)
	#write.fasta(chimp, 'Chimp', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)
	#write.fasta(gorilla, 'Gorilla', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)
	#write.fasta(orangutan, 'Orangutan', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)
	#write.fasta(gibbon, 'Gibbon', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)

	write.fasta(hcAnc, 'HC', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)
	write.fasta(hcgoAnc, 'HCGo', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)
	write.fasta(gapeAnc, 'Great_Ape', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)
	write.fasta(apeAnc, 'Ape', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)
	write.fasta(allAnc, 'Anthropoids', paste0(tmpdir, pkid, '_ancestors.fa'), open = "a", nbchar = 70, as.string = FALSE)

	# Read the alignment
	fn = paste0(tmpdir, pkid, '_ancestors.fa')
	alStr = readDNAStringSet(paste0(tmpdir, pkid, '_ancestors.fa'))
	rm(fn)

	return(alStr)
}


togap = function(fn){
	tx = readLines(fn)
	tx2 = gsub(pattern = "\\*", replace = "-", x = tx)
	writeLines(tx2, con=fn)
}



giveAnc = function(ancMlRes, node, prob = 0.75){

	# This order is used
	codes1 = c('A', 'C', 'G', 'T')

	# Probabilities of A, C, G, T per different site pattern
	# If no nucleotide's sequence probability is >0.75, assign N.
	# N will be ignored for substitution counting and will not
	# result in false positives for motif evolution analysis.
	# If all probabilities are equal and the sum of probabilities are 4,
	# the sequence was originally N, so assign it again.

	codesFinal = apply(ancMlRes[[node]], 1, function(x){
					if(all(x < prob) | sum(x) == 4){'N'}
					else{ codes1[which(x > prob)] }})

	ancInd = attributes(ancMlRes)$index
	ancSeq = codesFinal[ancInd]
	ancSeq = paste0(ancSeq, collapse = '')
	ancSeq = DNAString(ancSeq)

	#sum(strsplit(as.character(ancSeq), '')[[1]] == 'N')

	return(ancSeq)
}


wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}


geneOvEnrDep = function(gnL1, gnL2, bcg, plot = T, fn = 'Gene_Enrichment_Plot', wd = 20, hg = 10, xlabel = '', ylabel = ''){

	require(ggpubr)
	require(dplyr)

	# Fisher's exact test enrichment
	resL = list()
	for(i in 1:length(gnL2)){

		resL[[i]] = lapply(1:length(gnL1), function(x){
			A = gnL1[[x]]
			B = gnL2[[i]]
			q1 = bcg - length(union(A,B))
			q2 = length(setdiff(A,B))
			q3 = length(setdiff(B,A))
			q4 = length(intersect(A,B))
			mat = matrix(c(q1,q2,q3,q4), nrow=2)
			res = fisher.test(mat)
			or = res$estimate
			pval = res$p.value
			tmp = data.frame(var1 = names(gnL1)[x], var2 = names(gnL2)[i], OddsRatio = or, pval = pval)
			tmp
		}) %>% do.call(rbind, .)

		resL[[i]]$FDR = p.adjust(resL[[i]]$pval, method = 'fdr')
	}

	toplot = do.call(rbind, resL)
	rownames(toplot) = NULL

	# FDR for scaling
	toplot$log10FDR = -log10(toplot$FDR)

	# FDR for labeling
	toplot$sign_label = formatC(toplot$FDR, format = "e", digits = 1)

	# Asterisk for labeling
	toplot$sign_label_2 = ifelse(toplot$FDR < 0.05, '*', '')

	# Plot or return the data frame
	if(plot == T){

	pdf(paste0(fn, '.pdf'), width = wd, height = hg)
		print( ggscatter(toplot, x = 'var1', y = 'var2', color = 'OddsRatio', size = 'log10FDR') +
		geom_label(data = toplot, aes(label = sign_label), color="black",
			label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
		labs(x=xlabel, y=ylabel) +
		scale_size_continuous(range = c(5,25)) +
		scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red') +
		theme_classic() +
		theme(text = element_text(size=20, face = 'bold')) +
		geom_text(aes(label = sign_label_2), vjust = 1.2, colour = "darkgreen", fontface = 'bold', size = 20 ) +
		rotate_x_text(45) )
	dev.off()

	}

	if(plot == F){
		return(toplot)
	}
}


