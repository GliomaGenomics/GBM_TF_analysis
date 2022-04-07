#library(copynumber)

patients<-scan('patient_lists/glass_gbm_idhwt_rt_tmz_local+dis_cna.txt', what = character())
meta<-read.delim('reports/jarid2_results/outputs_actual_glass_1000_JARID2_results+dis.tsv', col.names=c('patient_id','value'), header=FALSE, row.names=1)

datain<-read.table('copy_number/variants_titan_seg_filtered.txt', header=FALSE)
colnames(datain)<-c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
datain$Sample<-substring(datain$Sample,1,15)

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

up_prim<-datain[datain$Sample %in% up_prim_pat,]
down_prim<-datain[datain$Sample %in% down_prim_pat,]
up_recu<-datain[datain$Sample %in% up_recu_pat,]
down_recu<-datain[datain$Sample %in% down_recu_pat,]

write.table(up_prim,file="copy_number/up_prim_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')
write.table(up_recu,file="copy_number/up_recu_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')
write.table(down_prim,file="copy_number/down_prim_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')
write.table(down_recu,file="copy_number/down_recu_cn.txt",row.names = FALSE,quote = FALSE, sep='\t')


#thresh=0.1
#pdf("copy_number/plotFreq.pdf")
#plotFreq(up_prim,thres.gain=thresh,thres.loss=-thresh)
#plotFreq(down_prim,thres.gain=thresh,thres.loss=-thresh)
#plotFreq(up_recu,thres.gain=thresh,thres.loss=-thresh)
#plotFreq(down_recu,thres.gain=thresh,thres.loss=-thresh)
#plotAberration(up_prim,thres.gain=thresh,thres.loss=-thresh)
#plotAberration(down_prim,thres.gain=thresh,thres.loss=-thresh)
#plotAberration(up_recu,thres.gain=thresh,thres.loss=-thresh)
#plotAberration(down_recu,thres.gain=thresh,thres.loss=-thresh)
#dev.off()
