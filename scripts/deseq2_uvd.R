library("DESeq2")

patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local.txt', what = character())

nes<-read.delim('reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv', col.names=c('patient_id','value'), header=FALSE)
dnes<-subset(nes,value < 0)
unes<-subset(nes,value > 0)
dpatients<-patients[match(dnes$patient_id,as.vector(patients))]
upatients<-patients[match(unes$patient_id,as.vector(patients))]
dpatients<-dpatients[!is.na(dpatients)]
upatients<-upatients[!is.na(upatients)]

datain<-read.table('original_data/PvR_geneCounts_all_LS_23062021.txt.txt', header=TRUE)

ddat=data.frame(EnsID=datain$EnsID)
for (p in dpatients){
ddat[paste(p,"_p",sep="")]<-datain[,grep(paste(p,"_P",sep=""), colnames(datain))]
ddat[paste(p,"_r",sep="")]<-datain[,grep(paste(p,"_R",sep=""), colnames(datain))]
}
rownames(ddat)<-ddat$EnsID
ddat<-subset(ddat, select = -EnsID)

udat=data.frame(EnsID=datain$EnsID)
for (p in upatients){
udat[paste(p,"_p",sep="")]<-datain[,grep(paste(p,"_P",sep=""), colnames(datain))]
udat[paste(p,"_r",sep="")]<-datain[,grep(paste(p,"_R",sep=""), colnames(datain))]
}
rownames(udat)<-udat$EnsID
udat<-subset(udat, select = -EnsID)


dmeta<-data.frame(condition=rep(c("primary","recurrent"),times=length(dpatients)), patient=rep(dpatients,times=rep(2,times=length(dpatients))), row.names=colnames(ddat))
umeta<-data.frame(condition=rep(c("primary","recurrent"),times=length(upatients)), patient=rep(upatients,times=rep(2,times=length(upatients))), row.names=colnames(udat))

ddds <- DESeqDataSetFromMatrix(countData = ddat, colData = dmeta, design = ~ patient + condition)
udds <- DESeqDataSetFromMatrix(countData = udat, colData = umeta, design = ~ patient + condition)
ddds <- DESeq(ddds)
udds <- DESeq(udds)
save(ddds,file="deseq2_uvd/ddds.Rdata")
save(udds,file="deseq2_uvd/udds.Rdata")

dres <- results( ddds )
ures <- results( udds )

dres$pvalue[mcols(ddds)$betaConv %in% FALSE]<-rep(c(NA),times=length(dres$pvalue[mcols(ddds)$betaConv %in% FALSE]))
ures$pvalue[mcols(udds)$betaConv %in% FALSE]<-rep(c(NA),times=length(ures$pvalue[mcols(udds)$betaConv %in% FALSE]))
dres$padj[mcols(ddds)$betaConv %in% FALSE]<-rep(c(NA),times=length(dres$pvalue[mcols(ddds)$betaConv %in% FALSE]))
ures$padj[mcols(udds)$betaConv %in% FALSE]<-rep(c(NA),times=length(ures$pvalue[mcols(udds)$betaConv %in% FALSE]))


pdf("deseq2_uvd/plots.pdf")

plotMA( dres,ylim = c(-3, 3) )
# Make a basic volcano plot
with(dres, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(dres, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dres, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

plotMA( ures,ylim = c(-3, 3) )
# Make a basic volcano plot
with(ures, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(ures, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(ures, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

dev.off()

write.table(dres, file="deseq2_uvd/results_down.txt")
write.table(ures, file="deseq2_uvd/results_up.txt")
