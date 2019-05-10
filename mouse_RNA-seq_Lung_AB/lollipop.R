
lollipop_GSEA =function(gseadf, pal = colorRampPalette(c("red", "gray", "blue"))(20), show.max = 10, toLower = F, legend.pos = "topleft", ...)
{
	df = data.frame("NES" = gseadf$NES, "padj" = gseadf$padj, "pathway" = gseadf$pathway)
	df = df[order(df$NES, decreasing = T), ]
	df.sub = df[c(1:show.max, (nrow(df)-(show.max-1)):nrow(df)), ]
	df.sub$pathway = sapply(df.sub$pathway, function(x) gsub(x, pattern = "_", replacement = " "))
	if (toLower == T) df.sub$pathway = tolower(df.sub$pathway)
	bpl = rev(barplot(df.sub$NES, horiz = T, plot = F))
	plot(0,0,cex = 0, ylim = range(c(0,bpl))*1.10, xlim = range(df.sub$NES)*1.50, xlab = "Normalized Enrichment Score", ylab = NA, yaxt = "n", ...)
	segments(x0 = rep(0, length(bpl)), y0 = bpl, x1 = df.sub$NES, y1 = bpl)
	segments(x0 = 0, x1 = 0, y0 = bpl, y1 = bpl[length(bpl)])
	df.sub$logP = -log10(df.sub$padj)
	if(length(table(df.sub$logP)) < 3) {bins = 2} else {bins = 3}
	df.sub$cex = prop_cex(df.sub$logP, bins, minc = 2, maxc = 3)
	df.sub$color = colorKey(df.sub$NES, pal)
	points(df.sub$NES, bpl, cex = df.sub$cex, pch = 21, bg = df.sub$color)
	text(x = 0, y = bpl[1:show.max], labels = df.sub$pathway[1:show.max], pos = 2, cex = 0.6)
	text(x = 0, y = bpl[(show.max+1):length(bpl)], pos = 4, cex = 0.6, labels = df.sub$pathway[(show.max+1):length(bpl)])
	text(x = df.sub$NES[1:show.max]+0.1, y = bpl[1:show.max], pos = 4, adj = 1.5, cex = 0.6, labels = formatC(df.sub$padj[1:show.max], digits = 2, format = "e"))
	text(x = df.sub$NES[(show.max+1):length(bpl)]-0.1, y = bpl[(show.max+1):length(bpl)], pos = 2, cex = 0.6, labels = formatC(df.sub$padj[(show.max+1):length(bpl)], digits = 2, format = "e"))
	pvalbounds = vector()
	cexvec = seq(2,3,length.out = bins)
	for(i in 1:bins)
	{
		tmp = df.sub[which(df.sub$cex == cexvec[i]),]
		pvalbounds[i] = min(tmp$padj)
	}
	legend("top", bty = "n", pch = 16, pt.cex = cexvec, cex = 1, ncol = 3, 
		legend = c(paste("p <=", formatC(pvalbounds[1], digits = 2, format = "e")), paste("p <=", formatC(pvalbounds[2], digits = 2, format = "e")), paste("p <=", formatC(pvalbounds[3], digits = 2, format = "e"))))
}
