library(ggplot2)
library(ggrepel)

load('analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE.RData')
files<-list.files(path ='cell_lines/tables/',pattern='*log2fc_[^m]')
data<-data.frame()
for (f in files){
print(f)      
if (dim(data)[2]!=0){
data1<-read.table(paste('cell_lines/tables/',f,sep=""), header=TRUE, row.names=1)
data2<-merge(data,data1,by.x='row.names',by.y='row.names')
data<-data2
row.names(data)<-data$Row.names
data<-subset(data, select = -Row.names)
}else{
data <- read.table(paste('cell_lines/tables/',f,sep=""), header=TRUE, row.names=1)
}
}

cell_lines<-colnames(data)


new=t(data)

patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local_+_cell_lines.txt', what = character())
meta<-read.delim('reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv', col.names=c('patient_id','value'), header=FALSE)
meta<-meta[match(as.vector(patients), meta$patient_id),]

newpca<- new%*%pca$rotation
pca$x <-rbind(pca$x,newpca)


pdf(paste("cell_lines/pca/pca_cell_lines_log2fc.pdf", sep=""))

ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% cell_lines,"cell_line","sample")), size = 2.88, alpha=5/6) +
  geom_text_repel(aes(label=ifelse(names(pca$x[,1]) %in% cell_lines,as.character(names(pca$x[,1])),'')),hjust=0,vjust=0, max.time=10, max.overlaps = Inf) +
  theme(legend.title=element_blank())+
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow")+
  geom_text_repel(aes(label=ifelse(names(pca$x[,1]) %in% cell_lines,as.character(names(pca$x[,1])),'')),hjust=0,vjust=0, max.time=10, max.overlaps = Inf) +
  labs(color = "NES")+
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% cell_lines,"cell_line","sample")), size = 2.88, alpha=5/6) +
  theme(legend.title=element_blank())+
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow")+
  labs(color = "NES")+
  labs(x = "PC1", y= "PC2")

dev.off()
