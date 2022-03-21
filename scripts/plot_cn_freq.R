library(copynumber)

patients<-scan('patient_lists/glass_gbm_idhwt_rt_tmz_local+stead_cna.txt', what = character())
meta<-read.delim('reports/jarid2_results/outputs_actual_glass_1000_JARID2_results+stead.tsv', col.names=c('patient_id','value'), header=FALSE, row.names=1)

datain<-read.table('glass_data/variants_titan_seg_filtered.txt', header=FALSE)
colnames(datain)<-c("sample","chrom","start.pos","end.pos","n.probes","mean")
datain$sampleID<-substring(datain$sample,1,15)
datain$arm<-datain$sampleID
datain<-datain[,c(7,2,8,3,4,5,6)]

up_prim_pat<-c()
down_prim_pat<-c()
up_recu_pat<-c()
down_recu_pat<-c()

for (p in patients){
if(meta[gsub('-','.',p),'value']>0){
up_prim_pat<-append(up_prim_pat,paste(p,'-TP',sep=''))
up_recu_pat<-append(up_recu_pat,paste(p,'-R1',sep=''))
}else{
down_prim_pat<-append(down_prim_pat,paste(p,'-TP',sep=''))
down_recu_pat<-append(down_recu_pat,paste(p,'-R1',sep=''))
}
}

up_prim<-datain[datain$sampleID %in% up_prim_pat,]
down_prim<-datain[datain$sampleID %in% down_prim_pat,]
up_recu<-datain[datain$sampleID %in% up_recu_pat,]
down_recu<-datain[datain$sampleID %in% down_recu_pat,]

write.table(up_prim,file="copy_number/up_prim_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')
write.table(up_recu,file="copy_number/up_recu_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')
write.table(down_prim,file="copy_number/down_prim_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')
write.table(down_recu,file="copy_number/down_recu_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')


thresh=0.05
pdf(paste("copy_number/plotFreq",thresh,".pdf",sep=""))
plotFreq(up_prim,thres.gain=thresh,thres.loss=-thresh,xlab='up_prim',chrom=6)
plotFreq(down_prim,thres.gain=thresh,thres.loss=-thresh,xlab='down_prim',chrom=6)
plotFreq(up_recu,thres.gain=thresh,thres.loss=-thresh,xlab='up_recu',chrom=6)
plotFreq(down_recu,thres.gain=thresh,thres.loss=-thresh,xlab='down_recu',chrom=6)

chro=6
plotAberration(up_prim,thres.gain=thresh,thres.loss=-thresh,xlab='up_prim',chrom=chro)
plotAberration(down_prim,thres.gain=thresh,thres.loss=-thresh,xlab='down_prim',chrom=chro)
plotAberration(up_recu,thres.gain=thresh,thres.loss=-thresh,xlab='up_recu',chrom=chro)
plotAberration(down_recu,thres.gain=thresh,thres.loss=-thresh,xlab='down_recu',chrom=chro)
dev.off()
