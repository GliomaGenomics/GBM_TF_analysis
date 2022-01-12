library("ggplot2")
library(uwot)

patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local.txt', what = character())
datain<-read.table('original_data/combat_pc_TPM_23062021.txt', header=TRUE)
datain$EnsID<-substring(datain$EnsID,1,15)
dat=data.frame(EnsID=datain$EnsID)
for (p in patients){
dat[paste(p,"_p",sep="")]<-datain[,grep(paste(p,"_P",sep=""), colnames(datain))]
dat[paste(p,"_r",sep="")]<-datain[,grep(paste(p,"_R",sep=""), colnames(datain))]
}
rownames(dat)<-dat$EnsID
dat<-subset(dat, select = -EnsID)
dat<-t(dat)

genes<-scan('gene_lists/JARID2_bound_genes.txt', what = character())
filter<-intersect(genes,colnames(dat))
#dat=dat[ , filter]

meta<-read.delim('reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv', col.names=c('patient_id','value'), header=FALSE)
meta<-meta[match(as.vector(patients), meta$patient_id),]

up<-meta[meta$value > 0 ,]$patient_id
down<-meta[meta$value < 0 ,]$patient_id

sample<-data.frame(matrix(ncol = 2, nrow = 0))
colnames(sample)<-c('patient_id','sample')

for(p in up){
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_p',sep=''),sample="Primary_up"))
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_r',sep=''),sample="Recurrent_up"))
}
for(p in down){
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_p',sep=''),sample="Primary_down"))
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_r',sep=''),sample="Recurrent_down"))
}

pca<- prcomp(dat, center = TRUE, scale= TRUE)
um<-as.data.frame(umap(dat,fast_sgd = TRUE, n_neighbors = 50, learning_rate = 0.5, init = "random", scale=TRUE))

cols<-merge(pca$x,sample,by.x=0, by.y='patient_id')

pdf("analysis/all_samples/pca_combined.pdf")


ggplot(um, aes(x=V1, y=V2, colour=cols$sample)) +
  theme_classic() +
  geom_point(size = 2.88, alpha=5/6)

ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2], colour=cols$sample)) +
  theme_classic() +
  geom_point(size = 2.88, alpha=5/6)
ggplot(as.data.frame(pca$x), aes(x=pca$x[,3], y=pca$x[,4], colour=cols$sample)) +
  theme_classic() +
  geom_point(size = 2.88, alpha=5/6)
ggplot(as.data.frame(pca$x), aes(x=pca$x[,5], y=pca$x[,6], colour=cols$sample)) +
  theme_classic() +
  geom_point(size = 2.88, alpha=5/6)




dev.off()


