library(ggplot2)
#library(ggrepel)

load('analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE.RData')
#load('analysis/pca/pca_gbm_idhwt_rt_tmz_local_JARID2_log2fc_all_FALSE.RData')
#load('analysis/pca/pca_gbm_idhwt_rt_tmz_local_le50_log2fc_all_FALSE.RData')
#load('analysis/pca/pca_gbm_idhwt_rt_tmz_local_le70_log2fc_all_FALSE.RData')

#datat <- read.table("tissue_slices/tables/log2fc_tissue_slices.txt", header=TRUE, row.names=1,check.names=FALSE)
#tissue_slices<-colnames(datat)
#datac <- read.table("cell_lines/tables/log2fc_cell_lines.txt", header=TRUE, row.names=1,check.names=FALSE)
#cell_lines<-colnames(datac)
#data<-merge(datat,datac,by.x='row.names',by.y='row.names',all=TRUE)

data <- read.table("cell_lines/tables/log2fc_cell_lines.txt", header=TRUE, row.names=1,check.names=FALSE, fill=TRUE)
cell_lines<-colnames(data)

data[is.na(data)] <- 0
#row.names(data)<-data$Row.names
#data<-subset(data, select = -Row.names)

mdata<-data.frame(row.names(pca$rotation))
datafilt<-merge(data,mdata,by.x='row.names',by.y=1,all.y=TRUE)
row.names(datafilt)<-datafilt$Row.names
data<-subset(datafilt, select = -Row.names)
data[is.na(data)] <- 0
new=t(data)
newpca<- new%*%pca$rotation
pca$x <-rbind(pca$x,newpca)

#save(pca,file="tissue_slices/pca/pca_tissue_slices_log2fc_data.Rdata")
save(pca,file="cell_lines/pca/pca_cell_lines_log2fc_data.Rdata")
write.table(pca$x,file="cell_lines/pca/pca_cell_lines_log2fc_PC_scores.txt", quote=FALSE, col.names=NA)
write.table(pca$rotation,file="cell_lines/pca/pca_cell_lines_log2fc_PC_loadings.txt", quote=FALSE, col.names=NA)

if (TRUE) {stop("End")}


#patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local_+_cell_lines+tissue_slices.txt', what = character())
#meta<-read.delim('tissue_slices/reports/outputs_actual_1000_JARID2_results_cell_tissue.tsv', col.names=c('patient_id','value'), header=FALSE)
#meta<-meta[match(as.vector(patients), meta$patient_id),]
#cell_tissue<-append(cell_lines,tissue_slices)

meta<-read.delim('cell_lines/reports/combined/cell_lines_actual_nes_JARID2.txt', col.names=c('patient_id','value'), header=FALSE)


pdf(paste("tissue_slices/pca/pca_tissue_slices_log2fc.pdf", sep=""))



ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% cell_lines,"cell_line",ifelse(names(pca$x[,1]) %in% tissue_slices,"tissue_slices","sample"))), size = 2.88, alpha=5/6) +
  geom_text_repel(aes(label=ifelse(names(pca$x[,1]) %in% cell_lines,as.character(names(pca$x[,1])),ifelse(names(pca$x[,1]) %in% tissue_slices,as.character(names(pca$x[,1])),''))),hjust=0,vjust=0, max.time=10, max.overlaps = Inf) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(color = meta[match(row.names(pca$x),meta$patient_id),]$value), size = 2.88, alpha=3/6) +
  scale_colour_gradient(low = "blue", high = "yellow") +
  geom_text_repel(aes(label=ifelse(names(pca$x[,1]) %in% cell_tissue,as.character(names(pca$x[,1])),'')),hjust=0,vjust=0, max.time=10, max.overlaps = Inf) +
  labs(color = "NES") +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% cell_lines,"cell_line",ifelse(names(pca$x[,1]) %in% tissue_slices,"tissue_slices","sample"))), size = 2.88, alpha=5/6) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(shape = ifelse(names(pca$x[,1]) %in% cell_lines,"cell_line",ifelse(names(pca$x[,1]) %in% tissue_slices,"tissue_slices","sample")), color = meta[match(row.names(pca$x),meta$patient_id),]$value), size = 1.5, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow") +
  labs(color = "NES",shape="sample") +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(shape = ifelse(names(pca$x[,1]) %in% cell_lines,"cell_line",ifelse(names(pca$x[,1]) %in% tissue_slices,"tissue_slices","sample")), color = meta[match(row.names(pca$x),meta$patient_id),]$value), size = 1, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow") +
  labs(color = "NES",shape="sample") +
  labs(x = "PC1", y= "PC2")+
  xlim(-40,40)
ggplot(as.data.frame(pca$x), aes(x=pca$x[,3], y=pca$x[,4])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pca$x[,1]) %in% cell_lines,"cell_line",ifelse(names(pca$x[,1]) %in% tissue_slices,"tissue_slices","sample"))), size = 2.88, alpha=5/6) +
  theme(legend.title=element_blank()) +
  labs(x = "PC3", y= "PC4")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,3], y=pca$x[,4])) +
  theme_classic() +
  geom_point(aes(color = meta[match(row.names(pca$x),meta$patient_id),]$value), size = 2.88, alpha=3/6) +
  scale_colour_gradient(low = "blue", high = "yellow") +
  labs(color = "NES") +
  labs(x = "PC3", y= "PC4")

dev.off()
