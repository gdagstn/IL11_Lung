##############################################################################################################
# Function to get TPM values from a count matrix, given the length of each features 						 #
# Arguments:																								 #	
#	counts: count matrix (only numeric values allowed)														 #
#	length: vector of lengths of features in basepairs. Must be in the same order as the count matrix		 #
# Outputs: returns the TPM matrix																			 #
##############################################################################################################

get_tpm <- function(counts, length){
	l <- length/1000
	x <- counts/l
	scalingf <- colSums(x)
	scalingf <- scalingf/1e6
	for(i in 1:length(scalingf)) x[,i] <- x[,i]/scalingf[i]
	return(x)
}

##############################################################################################################
# Function to get results from a given DESeq object, specifying adjusted p-value and log2fold change         #
# thresholds. Must specify the contrast. Log2FoldChange is shrunk by default using lfcShrink from DESeq2	 #
# Arguments:																								 #	
#	dds: a DESeq object																						 #
#	p: maximum adjusted P values 																			 #
#	lfc: minimum absolute log2FoldChange. If shrink = T, will use the shrunk lfc to subset					 #
#	contrast: the contrast to be found in the DESeq object, e.g. c("treatment", "A", "B")					 #
#	shrink: logical, determines if lfcShrink must be used do transform log2foldchange values. Defaults to T	 #
# Outputs: returns a DESeq Result object																	 #
##############################################################################################################

get_results <- function(dds, p = 0.05, lfc = 1, contrast, shrink=T){

	x <- results(dds, contrast = contrast)
	if (shrink == T) {
		x <- lfcShrink(dds, contrast=contrast)
		x <- x[which(abs(x$log2FoldChange) > lfc),]
	}else{x <- x[which(abs(x$log2FoldChange) > lfc),]}
	x <- x[which(x$padj < p),]
	return(x)

}

##############################################################################################################
# Function to get results from a given DESeq object with interaction terms, specifying adjusted p-value      #
# and log2fold change thresholds. Must specify the interaction name. 										 #
# Arguments:																								 #	
#	dds: a DESeq object																						 #
#	p: maximum adjusted P values 																			 #
#	lfc: minimum absolute log2FoldChange. 																	 #
#	interaction: the interaction name to be found in the DESeq object, e.g. "treatmentA.patient0"			 #
# Outputs: returns a DESeq Result object																	 #
##############################################################################################################

get_interaction <- function(dds, p = 0.05, lfc = 1, interaction){

	x <- results(dds, name = interaction)
	x <- x[which(abs(x$log2FoldChange) > lfc),]
	x <- x[which(x$padj < p),]
	return(x)

}

##############################################################################################################
# Function to plot a PCA showing different features such as different Principal Components, different colors #
# symbols (pch) and labels on the plot.                          											 #
# Arguments:																								 #	
#	data: a normalized count matrix (e.g. the result of a rld transformation)								 #
#	comp1, comp2: the two principal components to be visualized together in the plot						 #
#	coldata: a dataframe containing information regarding the samples. must include colors and symbol		 #
#	colors: a character with the coldata$column containing the color vector used in the plot 				 #
#	symbol: a character with the coldata$column containing the pch vector used in the plot 	 				 #
#	ntop: a numeric value indicating the number of genes with displaying the most variance across samples.   #
#		defaults to 500.																					 #
#	lab1, lab2: a character with the coldata$column containing labels for each point. Labels will be         #
#		printed on the left and right hand of each point in the plot.										 #	
#	...: any other graphical parameter passed to plot()														 #
# Outputs: plots a PCA plot																	 				 #
##############################################################################################################


plot_PCA <- function(data, comp1, comp2, coldata, colors, symbol, ntop = 500, lab1 = NA, lab2 = NA, ...){
  rv <- rowVars(assay(data))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca_all <- (prcomp(t(assay(data)[select,])))
  pca_all_data <- pca_all$x
  pca_all_data <- cbind(as.data.frame(pca_all_data), coldata)
  percentVar <- pca_all$sdev^2 / sum(pca_all$sdev^2)
  plot(pca_all_data[,comp1], pca_all_data[,comp2], col=colors, pch=symbol, xlab=paste("PC", comp1, ": ", round(percentVar[comp1], 2)*100, "% variance", sep=""), ylab=paste("PC", comp2, ": ", round(percentVar[comp2], 2)*100, "% variance", sep=""), ...)
    text(pca_all_data[,comp1], pca_all_data[,comp2], pos=2, cex=0.7, labels=lab1)
    text(pca_all_data[,comp1], pca_all_data[,comp2], pos=4, cex=0.7, labels=lab2)
}

plot_PCA_2 <- function(data, comp1, comp2, coldata, colors, symbol, ntop = 500, lab1 = NA, lab2 = NA, ...){
	rv <- rowVars(data)
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
	pca_all <- (prcomp(t(data[select,])))
	pca_all_data <- pca_all$x
  pca_all_data <- cbind(as.data.frame(pca_all_data), coldata)
  percentVar <- pca_all$sdev^2 / sum(pca_all$sdev^2)
  plot(pca_all_data[,comp1], pca_all_data[,comp2], col=colors, pch=symbol, xlab=paste("PC", comp1, ": ", round(percentVar[comp1], 2)*100, "% variance", sep=""), ylab=paste("PC", comp2, ": ", round(percentVar[comp2], 2)*100, "% variance", sep=""), ...)
    text(pca_all_data[,comp1], pca_all_data[,comp2], pos=2, cex=0.7, labels=lab1)
    text(pca_all_data[,comp1], pca_all_data[,comp2], pos=4, cex=0.7, labels=lab2)
}


plot_PCA_3 <- function(data, comp1, comp2, coldata, colors, ntop = 500, outline = NA, evids = evid, evidsrow = evidrow, ...){
	rv <- rowVars(data)
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
	pca_all <- (prcomp(t(data[select,])))
	pca_all_data <- pca_all$x
	label <- evidsrow
  pca_all_data <- cbind(as.data.frame(pca_all_data), coldata)
  percentVar <- pca_all$sdev^2 / sum(pca_all$sdev^2)
  plot(pca_all_data[,comp1], pca_all_data[,comp2], bg=colors, pch=21, col = outline, lwd = 1, xlab=paste("PC", comp1, ": ", round(percentVar[comp1], 2)*100, "% variance", sep=""), ylab=paste("PC", comp2, ": ", round(percentVar[comp2], 2)*100, "% variance", sep=""), ...)
    text(pca_all_data[label,comp1], pca_all_data[label,comp2], pos = thigmophobe(x = pca_all_data[label,comp1], y= pca_all_data[label,comp2]), cex=0.7, labels=evids)
    
}
##############################################################################################################
# Function to plot Fold Change values in a barplot with error bars.  									     #
# Arguments:																								 #	
#	gene: a character value with a HGNC gene symbol 														 #
#	foldchanges_tpm: a data frame in which fold change values are stored. defaults to foldchanges_tpm 		 #
#	sed_tpm: a data frame in which error values are stored. defaults to sed_tpm								 #
#	...: any other graphical parameter passed to barplot()													 #
# Outputs: plots a barplot of FC values with error bars										 				 #
##############################################################################################################

plot_fc <- function(gene, foldchanges_tpm = foldchanges_tpm, sed_tpm = sed_tpm, ...){
	ensid <- conv[which(conv[,2] == gene),1]
	plot_ylim <- range(c(0,foldchanges_tpm[ensid,]+sed_tpm[ensid,]))
	barplot(as.numeric(foldchanges_tpm[ensid, 1:ncol(foldchanges_tpm)]), main=gene, border=NA, ...)
		bpl <- barplot(as.numeric(foldchanges_tpm[ensid, 1:ncol(foldchanges_tpm)]), plot=F)
	abline(h=0)
	axis(1, at=bpl, labels=colnames(foldchanges_tpm), tick=F)
	arrows(bpl, as.numeric(foldchanges_tpm[ensid,] - sed_tpm[ensid,]), bpl, as.numeric(foldchanges_tpm[ensid,] + sed_tpm[ensid,]), code=3, angle=90, length=0.08)
}

##############################################################################################################
# Function to plot TPM values in a barplot.					  											     #
# Arguments:																								 #	
#	gene: a character value with a HGNC gene symbol 														 #
#	tpm: a data frame in which TPM values are stored. defaults to tpm 		 								 #
#	...: any other graphical parameter passed to barplot()													 #
# Outputs: plots a barplot of TPM values													 				 #
##############################################################################################################

plot_tpm <- function(gene, tpm, ...){
	ensid <- conv[which(conv[,2] == gene),1]
	plot_ylim <- range(c(0,tpm[ensid,]))
	barplot(as.numeric(tpm[ensid,]), ylab='TPM', border=NA, ylim=plot_ylim, ...)
		bpl <- barplot(as.numeric(tpm[ensid, ]), plot=F)
	axis(1, at=bpl, labels=colnames(tpm), tick=F, cex.axis=0.6, las=3)
}

##############################################################################################################
# Function to derive standard deviation of mean ratios via error propagation. Uses the formula			     #
# SD(A/B) = A/B * sqrt((SD_A/A)^2 + (SD_B/B)^2) and returns values that are then used by the plot_fc 		 #
# to plot error bars for fold changes of means.	All input data frames must be in the same order				 #
# Arguments:																								 #	
#	avg: a data frame containing averaged values	 														 #
#	sd: a data frame containing standard deviation values 					 								 #
# 	fc: a data frame containing fold changes																 #		
#	a: a numeric value indicating the column in avg and sd containing values for A							 #
#	b: a numeric value indicating the column in avg and sd containing values for B 							 #
#	c: a numeric value indicating the colum in fc containing values for A/B 								 #
# Outputs: a vector of SD(A/B) values 														 				 #
##############################################################################################################

fc_errorprop <- function(avg, sd, fc, a, b, c){
 sd_ab <- sqrt((sd[,a]/avg[,a])^2 + (sd[,b]/avg[,b])^2)*fc[,c]
 return(sd_ab)
}

##############################################################################################################
# Function to rename colnames from the output of a standard NGS pipeline. Takes the underscores out of       #
# sample names to allow for a increased handling and readability. 		    								 #
# Arguments:																								 #	
#	table: a data frame whose colnames need to be replaced													 #
#	sep: the separator in the names, e.g. "_" or "[.]".       				 								 #
# 	fields: a numeric value indicating the fields from the split names to be used							 #		
# Outputs: a vector of characters to be used in lieu of the original colnames.				 				 #
##############################################################################################################


renametable <- function(table, sep, fields){
  j <- colnames(table)
  for (i in 1:ncol(table))
  {
    a <- as.character(unlist(strsplit(j[i], sep)))
    j[i] <- paste(a[fields], collapse="_")
  }
  return(j)
  }



##############################################################################################################
#       																									 #	
# 																											 #
# 																											 #	
#										C O M I N G   S O O N												 #
#																											 #
#  																											 #		
# 																							 				 #
##############################################################################################################

plot_rnaribo <- function(gene){

rnalfc <- rep(0,5)
names(rnalfc) <- names(reslist_RNA)
ribolfc <- rep(0,5)
names(ribolfc) <- names(reslist_Ribo)
ribonly <- rep(0,5)
names(ribonly) <- names(reslist_Int)

for(i in names(rnalfc)) rnalfc[i] <- reslist_RNA[[i]][gene,7]
for(i in names(ribolfc)) ribolfc[i] <- reslist_Ribo[[i]][gene,7]
for(i in names(ribonly)) ribonly[i] <- reslist_Int[[i]][gene,2]

x <- rbind(rnalfc, ribolfc, ribonly)

plot(rnalfc, col="red", ylim=c(min(x),max(x)), main=conv[gene,2], ylab=NA, xlab=NA)
lines(rnalfc, col='red')
points(ribolfc, col="blue")
lines(ribolfc, col="blue")
points(ribonly, col="black")
lines(ribonly, col="black")
abline(h=0, col='gray', lty=2)
#legend(fill=c('red', 'blue', "black"), legend=c('RNA', 'Ribo', 'RiboOnly'), bty='n', border=NA, "top")
}

##############################################################################################################
# Function to plot a Volcano Plot from a DESeq2 results object.												 #	
# Arguments:																								 #
# 	res: the DESeq2 results object																			 #	
#	alpha: the FDR threshold to draw the horizontal line and colour points									 #
#	lfc: the absolute log2foldchange threshold to draw the vertical lines and colour points					 #
#  	maxdot: the number of top genes to label																 #		
# 	cols: a vector of colors to generate the color ramp for positive and negative fold changes				 #	
#	grid: logical: do you want to draw grid lines?															 #
#	labels: logical: do you want to label (maxdot) top genes?												 #	
#	scalefactor: numeric, how far apart are colours in the ramp. Must be > 1								 # 
#	...: arguments passed to the plot function																 #
#	cextest: size of the Labels	 																			 #
#	logZero: value with which FDR = 0 is substituted to avoid having infinite values. Defaults to 10e-300,   #
#			 but if is NA, values with FDR = 0 are kept out of the plot.	 								 #
# Outputs: a pretty volcano plot 																			 #	
# Requirements: plotrix and viridis                                  										 #
##############################################################################################################


plot_volcano <- function(res, alpha = 0.05, lfc = 1, maxdot = 30, cols = c(viridis(25, option="B")[5:10], viridis(25, option="B")[15:20]), grid = T, labels = T, scalefactor = 1.5, cextest = 1.5, logZero = 10e-300, ...){
	require(viridis)
	require(plotrix)
	
	col_lfc <- which(colnames(res) == "log2FoldChange")
	col_padj <- which(colnames(res) == "padj")

	if (scalefactor < 1)
		{
			print("Scale factor must be above 1! Defaulting to 1")
			scalefactor = 1
		}
	
	if (!is.na(logZero)) 
		{
			res$padj[which(res$padj == 0)]	= logZero 
		}	
		else 
		{
			res <- res[which(is.finite(-log10(res$padj))),]
		}
			
	res <- res[order(res$padj, decreasing=F),]

	topsigdf <- res[which(res$padj < alpha & abs(res$log2FoldChange) > lfc),]

	topsigdf <- topsigdf[order(abs(topsigdf$log2FoldChange), decreasing =T),]

	topsig <- rownames(topsigdf)

	dotcol <- colorRampPalette(cols)(length(topsig) * scalefactor)

	dotcol_neg <- dotcol[1:length(which(res[topsig,col_lfc] < 0))]

	dotcol_pos <- dotcol[(length(dotcol) - length(which(res[topsig,col_lfc] > 0))) : length(dotcol)]

	if (length(topsig) < maxdot)
		{ 
			md <-  length(topsig)
		} 
	else 
		{
			md = maxdot
		}

	if ("gene" %in% colnames(res))
		{	
			genecol <- which(colnames(res) == "gene")
			gene_label <- res[topsig,genecol][1:md]
		}	
	else
		{
			gene_label <- rownames(res[topsig,])[1:md]
		}

	plot(
			x = res$log2FoldChange, 
			y = -log10(res$padj), 
			col = "gray", 
			pch = 16, 
			cex = 0,
			ylab = "-log10(FDR)",
			xlab = "log2(Fold Change)",
			...
			)

	if (grid == T)
		
		{
			abline(
					h = seq(0, round(max(-log10(res$padj), na.rm = T),0), by = 2), 
					lwd = 0.5,  
					col = "gray"
					)
			abline(
					v = seq(range(round(res$log2FoldChange, 0))[1], range(round(res$log2FoldChange, 0))[2], by = 0.5), 
					lwd = 0.5,  
					col = "gray"
					)
		}

	points(
			x = res$log2FoldChange, 
			y = -log10(res$padj), 
			col = "gray", 
			pch = 16, 
			cex = 0.5
			)

	points(
			x = res[topsig,col_lfc][order(res[topsig,col_lfc], decreasing=F)], 
			y = -log10(res[topsig,col_padj][order(res[topsig,col_lfc], decreasing=F)]), 
			col = c(dotcol_neg, dotcol_pos), 
			pch = 16, 
			cex = 0.6
			)

	if (labels == T)
		{
			points(
					x = res[topsig,col_lfc][1:md], 	
					y = -log10(res[topsig,col_padj])[1:md], 
					col = "black", 
					pch = 1, 
					cex = 1.2
					)
			text(
					x = res[topsig,col_lfc][1:md], 
					y = -log10(res[topsig,col_padj])[1:md],
					labels = gene_label,
					pos = plotrix::thigmophobe(x = res[topsig,2][1:md],y = -log10(res[topsig,col_padj])[1:md]) ,
					col = "black",
					cex = cextest
					)
		}

	abline(
			v = c(-lfc, lfc), 
			lwd = 1, 
			lty = 2, 
			col = "black"
			)

	abline(
			h = c(-log10(alpha)), 
			lwd = 1, 
			lty = 2, 
			col = "black"
			)
}

##############################################################################################################
#       																									 #	
# 																											 #
# 																											 #	
#										C O M I N G   S O O N												 #
#																											 #
#  																											 #		
# 																							 				 #
##############################################################################################################

categorize_dtg <- function(rna, ribo, ribonly, alpha = 0.05, ...)
	{	
		sig_rna_pos <- list()
		sig_rna_neg <- list()
		sig_ribo_pos <- list()
		sig_ribo_neg <- list()
		sig_int_pos <- list()
		sig_int_neg <- list()

		for(i in names(rna))
			{
				sig_rna_pos[[i]] <- rownames(rna[[i]])[which(rna[[i]]$padj < alpha & rna[[i]]$log2FoldChange > 0)]
				sig_rna_neg[[i]] <- rownames(rna[[i]])[which(rna[[i]]$padj < alpha & rna[[i]]$log2FoldChange < 0)]
			}

		for(i in names(ribo))
			{
				sig_ribo_pos[[i]] <- rownames(ribo[[i]])[which(ribo[[i]]$padj < alpha & ribo[[i]]$log2FoldChange > 0)]
				sig_ribo_neg[[i]] <- rownames(ribo[[i]])[which(ribo[[i]]$padj < alpha & ribo[[i]]$log2FoldChange < 0)]
			}

		for(i in names(ribonly))
			{
				sig_int_pos[[i]] <- rownames(ribonly[[i]])[which(ribonly[[i]]$padj < alpha & ribonly[[i]]$log2FoldChange > 0)]
				sig_int_neg[[i]] <- rownames(ribonly[[i]])[which(ribonly[[i]]$padj < alpha & ribonly[[i]]$log2FoldChange < 0)]
			}

		sig_rna_pos <- unlist(sig_rna_pos, use.names=F)			
		sig_rna_neg <- unlist(sig_rna_neg, use.names=F)

		sig_ribo_pos <- unlist(sig_ribo_pos, use.names=F)			
		sig_ribo_neg <- unlist(sig_ribo_neg, use.names=F)

		sig_int_pos <- unlist(sig_int_pos, use.names=F)			
		sig_int_neg <- unlist(sig_int_neg, use.names=F)


		amp_genes = Reduce(intersect, list(sig_rna_pos, sig_ribo_pos, sig_int_pos))
		weak_genes = Reduce(intersect, list(sig_rna_neg, sig_ribo_neg, sig_int_neg))
		fw_genes = setdiff(setdiff(intersect(union(sig_rna_pos, sig_rna_neg), union(sig_ribo_pos, sig_ribo_neg)), union(sig_int_pos, sig_int_neg)), union(amp_genes, weak_genes))
		tex_genes = setdiff(setdiff(intersect(union(sig_ribo_pos, sig_ribo_neg), union(sig_int_pos, sig_int_neg)), union(sig_rna_pos, sig_rna_neg)), union(union(fw_genes, amp_genes), weak_genes))
		buff_genes = setdiff(setdiff(intersect(union(sig_int_pos, sig_int_neg), union(sig_rna_pos, sig_rna_neg)), union(sig_ribo_pos, sig_ribo_neg)), union(union(union(amp_genes, weak_genes), fw_genes), tex_genes))

		return(list("Amplified" = amp_genes, "Weakened" = weak_genes, "Forwarded" = fw_genes, "Translation_exclusive" = tex_genes, "Buffered" = buff_genes))
	}

##############################################################################################################
# Function to plot a Volcano Plot from a DESeq2 results object.												 #	
# Arguments:																								 #
# 	merged_table: a dataframe with a log2FC value for every time point/condition							 #	
#	timepoint: the number of the column by which values are compared (e.g. 2 for the second time 			 #
#	point/condition)									 													 #
#	categories: a list of gene names divided into categories (named lists are better)					 	 #
#  	col: a vector of colors for the categories																 #		
#	square: logical: do you want to draw have the same scale on both x and y?								 #
#	zoom: numeric: the maximum lfc on both axes to draw a square plot (overrides the value of square)		 #
#	ribonly: logical: do you want to have ribonly (log2TE) on the y axis? This does not make sense and you	 #
#	will be warned		 																					 #	
#	...: arguments passed to the plot function																 #
# Outputs: a quadrant plot to compare fold changes between RNA and ribo										 #	
##############################################################################################################


plot_quadrants <- function(merged_table, timepoint, categories, col=c("slateblue", "magenta","limegreen", "deepskyblue", "orange"), square = T, ribonly = F, zoom = NA, ...)
	
	{	

		colorcats <- as.list(col) 
		names(colorcats) = names(categories)
		if (timepoint > ncol(merged_table)/3) 
		
			{
				print("Whoa there cowboy, that's not a timepoint")
				return()
			}

		if (ribonly == T)
			{
				multi = 2
				lab.y = "RibOnly LFC (log2 TE)"
				print("This does not actually make sense...WHY DO YOU DO THAT?")
			}
		else
			{
				multi = 1
				lab.y = "RIBO LFC"
			}

		if (square == T) 
			{
				lim.y = range(c(range(merged_table[,timepoint]), range(merged_table[,((ncol(merged_table)/3)*multi)+timepoint])))
				lim.x = lim.y
			}
				else
			{
				lim.y = range(merged_table[,((ncol(merged_table)/3)*multi)+timepoint])
				lim.x = range(merged_table[,timepoint])
			}

		if (is.na(zoom) == F & is.numeric(zoom) == T)
			{
				lim.y = c(-as.numeric(zoom), as.numeric(zoom))
				lim.x = c(-as.numeric(zoom), as.numeric(zoom))
			}
		else if	(is.na(zoom) == F & is.numeric(zoom) == F)
			{
				print("Not a valid zoom number (must be numeric)")
				return()
			}
	
		plot(x = merged_table[,timepoint], y = merged_table[,((ncol(merged_table)/3)*multi)+timepoint], col="gray", pch=16, cex=0.3, ylab = lab.y, xlab = "RNA LFC", xlim = lim.x, ylim = lim.y, ...)
		
		for (i in names(categories))
			{
				points(x = merged_table[categories[[i]],timepoint], y = merged_table[categories[[i]],((ncol(merged_table)/3)*multi)+timepoint], col=colorcats[[i]], pch=16, cex=0.8)
			}
		abline(h = 0)
		abline(v = 0)
		abline(0,1, lwd=0.5, lty=2)
		abline(0,-1, lwd=0.5, lty=2)	
		legend(col=col, pch=16, legend=names(colorcats), bty = "n", "topleft")
	}


##############################################################################################################
# Function to test a color palette																			 #	
# Arguments:																								 #
# 	colors: a vector of colors to test																		 #	
#	cex: numeric value with dot size (default is 2)															 #
# Outputs: a plot with squares coloured with the desired palette. Includes color names if they are less than #
# 		   20, and brightness and luminance plots to identify perceptually uniform schemes					 #
##############################################################################################################

test_cols <- function(colors, cex = 2)

	{
		l <- length(colors)
		
		lums <- vector()

		for (i in 1:l)
		{
			lums[i] <- col2rgb(colors[i])[1]*0.299 + col2rgb(colors[i])[2]*0.587 + col2rgb(colors[i])[3]*0.114
		}

		brightness <- vector()

		for (i in 1:l)
		{
			brightness[i] <- col2rgb(colors[i])[1]*0.2126 + col2rgb(colors[i])[2]*0.7152 + col2rgb(colors[i])[3]*0.0722
		}

		ytop <- max(range(c(-60,range(lums), range(brightness))) + 30)

		plot(
			x = 1:l, 
			y = rep(-15,l),
			ylim = c(-60,ytop), 
			cex = cex, 
			pch = 15, 
			col = colors, 
			xaxt = "n", 
			yaxt = "n",
			ylab = "luminance and brightness", 
			xlab = NA, 
			bty = "n", 
			main = paste("Color palette test, n = ", l, sep = "")
			)
		
		lines(
			x = 1:l,
			y = lums
			)

		

		lines(
			x = 1:l,
			y = brightness,
			col = "gray"
			)

		if(l <= 50)
		
		{
			points(
				x = 1:l,
				y = lums,
				pch = 16,
				col = colors,
				)
			
			points(
				x = 1:l,
				y = brightness,
				pch = 17,
				col = colors
				)
			
			legend(
				"top", 
				bty="n", 
				lwd=1, 
				col=c("black", "gray"), 
				pch=c(16, 17), 
				legend=c("Luminance", "Brightness")
				)
		}
		else if(l > 50)
		{
			legend(
				"top", 
				bty="n", 
				lwd=1, 
				col=c("black", "gray"), 
				legend=c("Luminance", "Brightness")
				)
		}
	

		ticks <- seq(0, max(lums), by = 20)
		
		axis(2, at=ticks, labels=ticks)

		if (l < 20) 
		
			{
				text(labels = colors, srt = 90, x = 1:l, y = -50, cex = 0.7)	
			}
		
	}

##############################################################################################################
# Function to plot log2(FC) values by chromosome position													 #	
# Arguments:																								 #
# 	res: the DESeq2 results object or similar data frame													 #	
#	chromosome: the vector (or column) containing the chromosome names in the results object       			 #
#	cols: vector of colors to be used to distinguish among chromosomes (repeated)							 #																				 #
#	lfc: numeric, lower log2(FC) threshold to draw															 #
#  	evid = character, the name(s) of the gene(s) to be highlighted. uses the ensembl conversion table and/or #
#		   needs a "gene" column in the results object      												 #
#	attr = column name in the results df containing a feature that needs to be displayed using a color key	 #		
# 	cexattr = column name in the results df containing a feature that needs to be displayed using 			 #
#			  proportional dot sizes																		 #
#	... = anything else passed to the plotting function 													 #
# Outputs:																									 #
#	A plot in which each gene is a dot, placed within its chromosomal coordinates (x coordinate) and log2(FC)#
#	(y coordinate), and whose color and/or dot size can be proportional to any other numeric value (e.g. FDR)#
##############################################################################################################

plot_chromolfc <- function(res, chromosome = res$chromosome, cols = c("orange", "gray"), lfc = 1, evid = NA, attr = NA, cexattr = NA, ...)

{		
		if(!is.na(evid))
		{
			evids = sapply(evid, function(x) conv[which(conv$gene == x),1])
			evid.df = as.data.frame(res[evids,])
		}

		chrsize <- vector()
		for (i in c(as.character(1:22), "X", "Y"))
		{
			chrsize[i] <- length(which(chromosome == as.character(i)))
		}

		lengths <- c(0,cumsum(chrsize))
		coords <- rep(lengths, each=2)[2:(length(lengths)*2)]
		gaps <- rep(seq(0,(length(lengths)*300),by=300),each=2)[1:((2*length(lengths))-1)]
		cpg <- coords + gaps
		cpglist <- list()
		
		for (i in seq(1,length(cpg),1)) 
		{
			cpglist[[i]] <- c(cpg[i], cpg[i+1])
		}

		cpglist <- cpglist[seq(1,length(cpg),2)]
		cpglist <- cpglist[1:length(cpglist)-1]
		names(cpglist) <- c(as.character(1:22), "X", "Y")
		colorvec <- rep(cols, length.out = length(chrsize))
		names(colorvec) <- names(cpglist)
		
		plot(
			x = 1:max(cpg), 
			y = rep(max(res$log2FoldChange)/2,max(cpg)), 
			cex=0, 
			ylab = "log2(Fold Change)",
			bty = "n",
			xlab = NA,
			xaxt = "n",
			las = 2,
			ylim = c(min(res$log2FoldChange, na.rm = T), max(res$log2FoldChange, na.rm = T)),
			...
			)

		midpoints <- vector()
		for (i in names(cpglist))
			{
				midpoints[i] <- cpglist[[i]][1] + (cpglist[[i]][2] - cpglist[[i]][1]) / 2
			}

		axis(
			at = midpoints, 
			1, 
			labels = c(as.character(1:22), "X", "Y"), 
			pos = min(res$log2FoldChange - 0.5), 
			cex.axis = 0.7, 
			tick = F
			)

		for(i in names(cpglist)) 
		{
			data <- as.data.frame(cbind(res[which(res$chromosome == i),2], rep(i, length(which(res$chromosome == i))), res$padj[which(res$chromosome == i)]), stringsAsFactors =F)
			colnames(data) <- c("lfc", "chromosome", "padj")
			rownames(data) <- rownames(res)[which(res$chromosome == i)]
			data$lfc <- as.numeric(data$lfc)
			data$padj <- as.numeric(data$padj)
			data$xcoord <-seq(cpglist[[i]][1],cpglist[[i]][2]-1,1)	
			if (!is.na(cexattr)) data$cexattr = res[which(res$chromosome == i),which(colnames(res) == cexattr)]
			if (!is.na(attr)) data$attribute = res[which(res$chromosome == i),which(colnames(res) == attr)]

			points(
				x = data$xcoord, 
				y = data$lfc, 
				pch = 16, 
				cex = 0.3, 
				col = "gray"
				)
			data.sub <- data[which(abs(data$lfc) > lfc),]

		if (nrow(data.sub) > 0)
			{
				if (!is.na(attr)) {colpoints = data.sub$attribute}
				else 
				{colpoints = colorvec[i]}
				if (!is.na(cexattr)) {cexpoints = data.sub$cexattr}
				else 
				{cexpoints = 0.7}
			points(
				x = data.sub$xcoord, 
				y = data.sub$lfc, 
				pch = 16, 
				cex = cexpoints,
				col = colpoints,
				#col = "black"
				)
			if(!is.na(evid))
				{	
					evid.now = evid.df[which(evid.df$chromosome == i),]
					data.sub.evid = data.sub[rownames(evid.now),]
				
				if(length(rownames(data.sub.evid)) > 0)
					{
						points(
							x = data.sub.evid$xcoord,
							y = data.sub.evid$lfc,	
							pch = 1, 
							cex = 1.5,
							col = "black",
							)
						text(
							x = data.sub.evid$xcoord,
							y = data.sub.evid$lfc,	
							cex = 0.6,
							col = "black",
							labels = conv[rownames(evid.now),2],
							pos = 2	
							)
					}
				}	

			} else next

		}

		
		abline(h = 0, lwd = 2, col = "black", lty = 2, lend = 3)
		abline(h = c(-lfc, lfc), lwd = 1, col = "red", lty = 2)

}

##############################################################################################################
# Function to generate proportional dot sizes given a vector of numerics									 #
# Arguments:																								 #
#	values = the vector of numerical values that needs to be represented									 #
#	bins = numeric, the amount of different dots to be used													 #
#	minc = numeric, the minimum dot size (defaults to 0.3)													 #
# 	maxc = numeric, the maximum dot size (defaults to 3)													 #
#	na.value = numeric, how to interpret NA values if present (defaults to 0)								 #
# 	inf.value = numeric, how to interpret Infinite values if present (defaults to 400)						 #
# Outputs:																									 #
#	A vector of dot size values binned within minc and maxc values 											 #
##############################################################################################################


prop_cex <- function(values, bins, minc = 0.3, maxc = 3, na.value = 0, inf.value = 400)
{	
	values[is.na(values)] = na.value
	values[!is.finite(values)] = inf.value
	
	ordered.values <- values[order(values, decreasing = T)]

	cex.values <- seq(minc, maxc, length.out = length(unique(values)))
	cex.frame <- as.data.frame(cbind(cex.values, lfc =rev(unique(ordered.values))))
	cex.binned <- seq(minc, maxc, length.out = bins)
	outdf <- as.data.frame(values)
	outdf$cex.values = sapply(outdf$values, function (x) cex.frame[which(cex.frame$lfc == x),1])
	cex.values.binned = as.data.frame(sapply(outdf$cex.values, function(x) cut(x, breaks = seq(minc, maxc, length.out = bins+1), include.lowest = T, labels = cex.binned)), stringsAsFactors = F)
	outdf$binned = as.numeric(levels(cex.values.binned[,1])[cex.values.binned[,1]])
 	return(outdf$binned)
}

##############################################################################################################
# Function to generate a heatmap-like color key																 #
# Arguments:																								 #
#	values = the vector of numerical values that needs to be represented									 #
#	pal = the vector of colors, ordered from min to max value												 #
# Outputs:																									 #
#	A vector of colors binned and assigned to each value 		 											 #
##############################################################################################################



colorKey <- function(values, pal = viridis(25, option = "B"))
{
	require(pheatmap)
	values_sc <- scale(values)
	bks <- pheatmap:::generate_breaks(values_sc, length(pal), center = F)
	cols <- pheatmap:::scale_colours(values_sc, col=pal, breaks=bks, na_col = 'gray')
	cols <- as.character(cols)
	return(cols)
}


#' Bubble map
#' @param valuedf data frame with numerics that will be mapped to colors
#' @param pvaluedf data frame with numerics that will be mapped to dot sizes. treated as a p-value data frame (values < 0.05 are filtered)
#' @param cex.binned logical: bin dot sizes or map them continously to the pvaluedf? Default is TRUE
#' @param cbins numeric, passed to bins in prop_cex
#' @param color_pal vector of characters indicating colors to be mapped
#' @param maplabel character passed to main in the plot
#' @param ... other arguments passed to plot()
#' @return a pretty bubble map where colors are mapped to valuedf, and dot sizes are mapped to pvaluedf. In the current version the function assumes that pvaluedf contains p values, so that all non-significant values are mapped to the same dot size (0.5)



bubbleMap <- function(valuedf, pvaluedf, cex.binned = T, cbins = 5, color_pal = colorRampPalette(c("slateblue", "gray", "orange"))(25), maplabel = 
	"Correlation", ...)
{

	pmap <- pheatmap(valuedf, silent = T)
	valuedf <- as.data.frame(valuedf[rev(pmap$tree_row$order), pmap$tree_col$order])
	pvaluedf <-  as.data.frame(pvaluedf[rev(pmap$tree_row$order), pmap$tree_col$order])

	#pvaluedf$id = rownames(pvaluedf)
	p2 = as.data.frame(t(pvaluedf))
	p2$func = rownames(p2)
	moltenp = suppressMessages({reshape::melt.data.frame(p2)})

	moltenp$logP = as.numeric(-log10(moltenp$value))
	moltenp$logP[moltenp$logP < 1.3] = 1
	if (cex.binned == T)
	{	
		moltenp$cexes = prop_cex(moltenp$logP, bins = cbins, minc = 1, maxc = 4.4)
	}
	else
	{
		moltenp$cexes = moltenp$logP/max(moltenp$logP)
		moltenp$cexes = scales::rescale(moltenp$cexes, from = c(min(moltenp$cexes[moltenp$logP > 1]), max((moltenp$cexes[moltenp$logP > 1]))), to = c(1, 4.4))
	}

	moltenp$cexes[moltenp$logP == 1] = 0.5

	pcex = as.data.frame(cbind(10^-moltenp$logP, moltenp$cexes))
	pcex = pcex[order(pcex[,1]),]
	colnames(pcex) = c("logP", "cexes")

	edgen = round(max(c(abs(min(as.numeric(unlist(valuedf)))), abs(max(as.numeric(unlist(valuedf)))))), digits = 1)
	
	c2 = as.data.frame(t(valuedf))
	c2$func = rownames(c2)
	moltenc = reshape::melt.data.frame(c2)
	colors = colorKey(values = c(- edgen, edgen, moltenc$value), pal = color_pal)
	colors = colors[3:length(colors)]
	moltenp$colors = colors

	nr = nrow(valuedf)
	nc = ncol(valuedf)
	xsc = (1:nc)/nc
	ysc = (1:nr)/nr

	coordf = expand.grid(1:length(unique(moltenp$func)), 1:length(unique(moltenp$variable)))
	plot(coordf$Var1/max(coordf$Var1), coordf$Var2/max(coordf$Var2), cex = 0, bty = "n", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, ylim = c(0,1.2), xlim = c(-0.1, 1.1), ...)
	#abline(h = coordf$Var2)
	points(coordf$Var1/max(coordf$Var1), coordf$Var2/max(coordf$Var2), cex = moltenp$cexes, bg = moltenp$colors, pch = 21, bty = 'n', lwd = 0.5)
	axis(1, at = unique(coordf$Var1)/max(coordf$Var1), labels = colnames(valuedf) , las = 2, cex.axis = 0.9)
	axis(4, at = unique(coordf$Var2)/max(coordf$Var2), labels = rownames(valuedf), las = 2, cex.axis = 0.7)


	rect_series = seq(0.05, 0.4, length.out = length(color_pal) + 1)
	
	for(q in 1:length(color_pal)) 
		{
			rect(xleft = rect_series[q], ybottom = 1.1, xright = rect_series[q+1], ytop = 1.15, col = color_pal[q], border = NA)
		} 

	rect(xleft = 0.05, xright = 0.4, ybottom = 1.1 ,ytop = 1.15)

	text(x = 0.05, y = 1.18, labels = paste("-", edgen))
	text(x = 0.4, y = 1.18, labels = edgen)
	text(x = 0.225, y = 1.18, pos = 3, labels = maplabel)

	
	cvec = vector()
	for(i in 1:length(unique(moltenp$cexes)))
	{
		cvec[i] <- pcex[which(pcex$cexes == unique(pcex$cexes)[i])[length(which(pcex$cexes == unique(pcex$cexes)[i]))],1]	
	}
	cvec = cvec[1:(length(cvec)-1)]
	cxvec = c("p > 0.05", paste("p <=", rev(formatC(cvec, format = "e", digits = 2))))

	if(cex.binned == T)
	{
		point_series = seq(0.45, 1, length.out = length(unique(moltenp$cex)))
		for(i in 1:length(point_series)) points(y = 1.1, x = point_series[i], cex = unique(moltenp$cex)[order(unique(moltenp$cex))][i], pch = 21, lwd = 0.5, bg = "gray")
		text(x = point_series, y = rep(1.12, length(point_series)), labels = cxvec, cex = 0.7, srt = 45, pos = 4)

	} else {
		if(is.na(cbins)) cbins = length(unique(moltenp$cex))
		point_series = seq(0.45, 1, length.out = cbins)
		ccex = unique(moltenp$cex)[order(unique(moltenp$cex))][seq(1, length(unique(moltenp$cex)), length.out = cbins)]
		ccex[cbins] = max(unique(moltenp$cex))
		ccex[1] = 0.5
		for(i in 1:length(point_series)) points(y = 1.1, x = point_series[i], cex = ccex[i], pch = 21, bg = "gray", lwd = 0.5)
			ccxvec = cxvec[seq(1, length(unique(moltenp$cex)), length.out = cbins)]
			ccxvec[cbins] = paste("p < = ", formatC(min(cvec), format = "e", digits = 2))
		text(x = point_series, y = rep(1.14, length(point_series)), labels = ccxvec, cex = 0.7, srt = 45, pos = 4)
	}

}

#' Lollipop plot
#' @param gseadf GSEA data frame with results from fgsea
#' @param pal character vector of colors to be used. defaults to a diverging palette of 20 colors between red, gray and blue
#' @param show.max numeric, number of top up- and down- enriched pathways. If the df has less than 2 * show.max entries, it will show all results. defaults to 10
#' @param toLower logical: convert all names of pathways to lowercase? Default is FALSE
#' @param ... other arguments passed to plot()
#' @return a lollipop plot in which dot size is proportional to the p value of the enrichment, color and position map to the normalized enrichment score (x axis)



lollipop_GSEA =function(gseadf, pal = colorRampPalette(c("red", "gray", "blue"))(20), show.max = 10, toLower = F, ...)
{
	df = data.frame("NES" = gseadf$NES, "padj" = gseadf$padj, "pathway" = gseadf$pathway)
	df = df[order(df$NES, decreasing = T), ]

	if((show.max*2) > nrow(df)) 
		{
			df.sub = df
		} else {
			df.sub = df[c(1:show.max, (nrow(df)-(show.max-1)):nrow(df)), ]
		}

	df.sub$pathway = sapply(df.sub$pathway, function(x) gsub(x, pattern = "_", replacement = " "))
	if (toLower == T) df.sub$pathway = tolower(df.sub$pathway)
	bpl = rev(barplot(df.sub$NES, horiz = T, plot = F))
	plot(0,0,cex = 0, ylim = range(c(0,bpl))*1.10, xlab = "Normalized Enrichment Score", ylab = NA, yaxt = "n", ...)
	segments(x0 = rep(0, length(bpl)), y0 = bpl, x1 = df.sub$NES, y1 = bpl)
	segments(x0 = 0, x1 = 0, y0 = bpl, y1 = bpl[length(bpl)])
	df.sub$logP = -log10(df.sub$padj)
	if(length(table(df.sub$logP)) < 3) {bins = 2} else {bins = 3}
	df.sub$cex = prop_cex(df.sub$logP, bins, minc = 2, maxc = 3)
	df.sub$cex[df.sub$logP < 1.3] = 1
	df.sub$color = colorKey(df.sub$NES, pal)
	points(df.sub$NES, bpl, cex = df.sub$cex, pch = 21, bg = df.sub$color)
	if((show.max*2) <= nrow(df))
		{
			text(x = 0, y = bpl[1:show.max], labels = df.sub$pathway[1:show.max], pos = 2, cex = 0.6)
			text(x = 0, y = bpl[(show.max+1):length(bpl)], pos = 4, cex = 0.6, labels = df.sub$pathway[(show.max+1):length(bpl)])
			text(x = df.sub$NES[1:show.max]+0.1, y = bpl[1:show.max], pos = 4, adj = 1.5, cex = 0.6, labels = formatC(df.sub$padj[1:show.max], digits = 2, format = "e"))
			text(x = df.sub$NES[(show.max+1):length(bpl)]-0.1, y = bpl[(show.max+1):length(bpl)], pos = 2, cex = 0.6, labels = formatC(df.sub$padj[(show.max+1):length(bpl)], digits = 2, format = "e"))
		}
	else
		{
			pos_labels = df.sub$pathway[which(df.sub$NES > 0)]
			neg_labels = df.sub$pathway[which(df.sub$NES < 0)]
			text(x = 0, y = bpl, labels = c(pos_labels, neg_labels), cex = 0.6, pos = c(rep(2,length(pos_labels)), rep(4, length(neg_labels))))
			text(x = c(df.sub$NES[which(df.sub$NES > 0)] + 0.1, df.sub$NES[which(df.sub$NES < 0)] - 0.1), y = bpl, pos = c(rep(4,length(pos_labels)), rep(2, length(neg_labels))), adj = c(rep(1.5,length(pos_labels)), rep(0.5, length(neg_labels))), cex = 0.6, labels = formatC(df.sub$padj, digits = 2, format = "e"))

		}


	}


