
plot_PCA_2 <- function(data, comp1, comp2, coldata, colors, symbol, ntop = 500, lab1 = NA, lab2 = NA, ...){
	rv <-sapply(data, 1, var)
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
	pca_all <- (prcomp(t(data[select,])))
	pca_all_data <- pca_all$x
  pca_all_data <- cbind(as.data.frame(pca_all_data), coldata)
  percentVar <- pca_all$sdev^2 / sum(pca_all$sdev^2)
  plot(pca_all_data[,comp1], pca_all_data[,comp2], col=colors, pch=symbol, xlab=paste("PC", comp1, ": ", round(percentVar[comp1], 2)*100, "% variance", sep=""), ylab=paste("PC", comp2, ": ", round(percentVar[comp2], 2)*100, "% variance", sep=""), ...)
    text(pca_all_data[,comp1], pca_all_data[,comp2], pos=2, cex=0.7, labels=lab1)
    text(pca_all_data[,comp1], pca_all_data[,comp2], pos=4, cex=0.7, labels=lab2)
}



renametable <- function(table, sep, fields){
  j <- colnames(table)
  for (i in 1:ncol(table))
  {
    a <- as.character(unlist(strsplit(j[i], sep)))
    j[i] <- paste(a[fields], collapse="_")
  }
  return(j)
  }



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



colorKey <- function(values, pal = viridis(25, option = "B"))
{
	require(pheatmap)
	values_sc <- scale(values)
	bks <- pheatmap:::generate_breaks(values_sc, length(pal), center = F)
	cols <- pheatmap:::scale_colours(values_sc, col=pal, breaks=bks, na_col = 'gray')
	cols <- as.character(cols)
	return(cols)
}




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
	} 	else 	{
		moltenp$cexes = moltenp$logP/max(moltenp$logP)
		moltenp$cexes = scales::rescale(moltenp$cexes, from = c(min(moltenp$cexes[moltenp$logP > 1]), max((moltenp$cexes[moltenp$logP > 1]))), to = c(1, 4.4))
	}

	moltenp$cexes[moltenp$logP == 1] = 0.8

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
	moltenp$pch = rep(21, nrow(moltenp))
	moltenp$pch[moltenp$logP == 1] <- 4

	nr = nrow(valuedf)
	nc = ncol(valuedf)
	xsc = (1:nc)/nc
	ysc = (1:nr)/nr

	coordf = expand.grid(1:length(unique(moltenp$func)), 1:length(unique(moltenp$variable)))
	plot(coordf$Var1/max(coordf$Var1), coordf$Var2/max(coordf$Var2), cex = 0, bty = "n", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, ylim = c(0,1.2), xlim = c(-0.1, 1.1), ...)
	#abline(h = coordf$Var2)
	points(coordf$Var1/max(coordf$Var1), coordf$Var2/max(coordf$Var2), cex = moltenp$cexes, bg = moltenp$colors, pch = as.numeric(moltenp$pch), bty = 'n', lwd = 0.5)
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

