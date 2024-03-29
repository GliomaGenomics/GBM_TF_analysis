library("DESeq2")

patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local.txt', what = character())

datain<-read.table('original_data/PvR_geneCounts_all_LS_23062021.txt.txt', header=TRUE)
dat=data.frame(EnsID=datain$EnsID)
for (p in patients){
dat[paste(p,"_p",sep="")]<-datain[,grep(paste(p,"_P",sep=""), colnames(datain))]
dat[paste(p,"_r",sep="")]<-datain[,grep(paste(p,"_R",sep=""), colnames(datain))]
}
rownames(dat)<-dat$EnsID
dat<-subset(dat, select = -EnsID)

meta<-data.frame(condition=rep(c("primary","recurrent"),times=length(patients)), patient=rep(patients,times=rep(2,times=length(patients))), row.names=colnames(dat))

dds <- DESeqDataSetFromMatrix(countData = dat, colData = meta, design = ~ patient + condition)
dds <- DESeq(dds)
dds <- dds[rowSums(counts(dds)) >= 10,]

save(dds,file="deseq2/dds.Rdata")

res <- results( dds )
res<-res[mcols(dds)$betaConv=='TRUE',]

pdf("deseq2/plots.pdf")
plotMA( res,ylim = c(-3, 3) )
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

write.table(res, file="deseq2/results.txt")
