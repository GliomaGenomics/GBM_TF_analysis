library(ggplot2)

load('analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE.RData')
new1<-read.table('cell_lines/tables/log2fc_A172.txt', header=TRUE, row.names=1)
new2<-read.table('cell_lines/tables/log2fc_GBM63.txt', header=TRUE, row.names=1)
new3<-read.table('cell_lines/tables/log2fc_GBM63_CUTRUN.txt', header=TRUE, row.names=1)
new4<-merge(new1,new2,by='row.names')
new<-merge(new4,new3,by.x='Row.names',by.y='row.names')
row.names(new)<-new$Row.names
new<-subset(new, select = -Row.names)
new=t(new)

#patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local.txt', what = character())
#meta<-read.delim('reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv', col.names=c('patient_id','value'), header=FALSE)
#meta<-meta[match(as.vector(patients), meta$patient_id),]

newpca<- new%*%pca$rotation
pca$x <-rbind(pca$x,newpca)


pdf(paste("cell_lines/pca/pca_cell_lines.pdf", sep=""))
screeplot(pca, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)

plot(pca$x[,1],pca$x[,2], xlab="PC1", ylab = "PC2", main = "PC1 / PC2 - plot")
plot(pca$x[,3],pca$x[,4], xlab="PC3", ylab = "PC4", main = "PC3 / PC4 - plot")
plot(pca$x[,5],pca$x[,6], xlab="PC5", ylab = "PC6", main = "PC5 / PC6 - plot")

ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,1])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),"cell_line","sample")), size = 2.88, alpha=5/6) +
  geom_text(aes(label=ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),as.character(names(pca$x[,1])),'')),hjust=0,vjust=0)+
  theme(legend.title=element_blank())
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),"cell_line","sample")), size = 2.88, alpha=5/6) +
  geom_text(aes(label=ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),as.character(names(pca$x[,1])),'')),hjust=0,vjust=0)+
  theme(legend.title=element_blank())
ggplot(as.data.frame(pca$x), aes(x=pca$x[,3], y=pca$x[,4])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),"cell_line","sample")), size = 2.88, alpha=5/6) +
  geom_text(aes(label=ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),as.character(names(pca$x[,1])),'')),hjust=0,vjust=0)+
  theme(legend.title=element_blank())
ggplot(as.data.frame(pca$x), aes(x=pca$x[,5], y=pca$x[,6])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),"cell_line","sample")), size = 2.88, alpha=5/6) +
  geom_text(aes(label=ifelse(names(pca$x[,1]) %in% c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','GBM63_CUTRUN_0','GBM63_CUTRUN_1'),as.character(names(pca$x[,1])),'')),hjust=0,vjust=0)+
  theme(legend.title=element_blank())
dev.off()
