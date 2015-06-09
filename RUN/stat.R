RData.name <- list.files(pattern=".*local_GATCcounts.RData")

for (every in RData.name) {
	load(file=every)
	assign(sub("(.*)_(edge|inner).*", "\\2", every, perl=T), reads2GATC$count)
	assign(sub("(.*)_(edge|inner).*", "\\2_name", every, perl=T), sub("(.*)_local_GATCcounts.RData", "\\1", every, perl=T))
}

bmp(filename=paste("Correlations_between_", edge_name, "_and_", inner_name, ".bmp", sep=""), width=800, height=800, units = "px")
par(mai=c(1.5, 1.5, 0.5, 0.5))
par(cex=1.3)
Cor.P <- round(cor(edge, inner, method="pearson", use="pairwise.complete.obs"), digits=2)
Cor.S <- round(cor(edge, inner, method="spearman", use="pairwise.complete.obs"), digits=2)
plot(x=edge, y=inner, cex=0.3, xlab=edge_name, ylab=inner_name, las=1, bty="l", pch=".", text(x=min(edge, na.rm=T), y=max(inner, na.rm=T)-400, adj=0, labels=c(paste("pearson = ", Cor.P, "\n\n", sep=""), paste("spearman = ", Cor.S, sep=""))))
dev.off()