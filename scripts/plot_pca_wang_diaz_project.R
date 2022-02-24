library(ggplot2)
library(ggrepel)

pc <- read.table("downloaded_data/NIHMS1579377-supplement-Supplemental_Table_S2_ens.txt", header=TRUE, row.names=1,check.names=FALSE)
datat <- read.table("tissue_slices/tables/recurrent_tissue_slices.txt", header=TRUE, row.names=1,check.names=FALSE)
tissue_slices<-colnames(datat)
datac <- read.table("cell_lines/tables/recurrent_cell_lines.txt", header=TRUE, row.names=1,check.names=FALSE)
cell_lines<-colnames(datac)

patients<-scan("patient_lists/gbm_idhwt_rt_tmz_local.txt", what = character())
datap<-read.table("tables/recurrent_all.txt", header=TRUE, row.names=1)
datap=datap[ , as.vector(patients) ]
cn<-row.names(datap)
row.names(datap)<-substr(cn,1,15)

data<-merge(datat,datac,by.x='row.names',by.y='row.names',all=TRUE)
row.names(data)<-data$Row.names
data<-subset(data, select = -Row.names)
data2<-merge(data,datap,by.x='row.names',by.y='row.names',all=TRUE)
data2[is.na(data2)] <- 0
row.names(data2)<-data2$Row.names
data<-subset(data2, select = -Row.names)

mdata<-data.frame(row.names(pc))
datafilt<-merge(data,mdata,by.x='row.names',by.y=1,all.y=TRUE)
row.names(datafilt)<-datafilt$Row.names
data<-subset(datafilt, select = -Row.names)
data[is.na(data)] <- 0
new=t(data)

patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local_+_cell_lines+tissue_slices.txt', what = character())
meta<-read.delim('reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv', col.names=c('patient_id','value'), header=FALSE)
meta<-meta[match(as.vector(patients), meta$patient_id),]
cell_tissue<-append(cell_lines,tissue_slices)

pc<-as.matrix(pc)
pcx<- new%*%pc

#pdf(paste("tissue_slices/pca/pca_tissue_sclice_log2fc_wang_diaz.pdf", sep=""))
pdf(paste("tissue_slices/pca/pca_tissue_sclice_recurrent_wang_diaz.pdf", sep=""))

ggplot(as.data.frame(pcx), aes(x=pcx[,1], y=pcx[,2])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pcx[,1]) %in% cell_lines,"cell_line",ifelse(names(pcx[,1]) %in% tissue_slices,"tissue_slices","sample"))), size = 2.88, alpha=5/6) +
  geom_text_repel(aes(label=ifelse(names(pcx[,1]) %in% cell_lines,as.character(names(pcx[,1])),ifelse(names(pcx[,1]) %in% tissue_slices,as.character(names(pcx[,1])),''))),hjust=0,vjust=0, max.time=10, max.overlaps = Inf) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pcx), aes(x=pcx[,1], y=pcx[,2])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=3/6) +
  scale_colour_gradient(low = "blue", high = "yellow") +
  geom_text_repel(aes(label=ifelse(names(pcx[,1]) %in% cell_tissue,as.character(names(pcx[,1])),'')),hjust=0,vjust=0, max.time=10, max.overlaps = Inf) +
  labs(color = "NES") +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pcx), aes(x=pcx[,1], y=pcx[,2])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pcx[,1]) %in% cell_lines,"cell_line",ifelse(names(pcx[,1]) %in% tissue_slices,"tissue_slices","sample"))), size = 2.88, alpha=5/6) +
  theme(legend.title=element_blank()) +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pcx), aes(x=pcx[,1], y=pcx[,2])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=3/6) +
  scale_colour_gradient(low = "blue", high = "yellow") +
  labs(color = "NES") +
  labs(x = "PC1", y= "PC2")
ggplot(as.data.frame(pcx), aes(x=pcx[,3], y=pcx[,4])) +
  theme_classic() +
  geom_point(aes(colour = ifelse(names(pcx[,1]) %in% cell_lines,"cell_line",ifelse(names(pcx[,1]) %in% tissue_slices,"tissue_slices","sample"))), size = 2.88, alpha=5/6) +
  theme(legend.title=element_blank()) +
  labs(x = "PC3", y= "PC4")
ggplot(as.data.frame(pcx), aes(x=pcx[,3], y=pcx[,4])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=3/6) +
  scale_colour_gradient(low = "blue", high = "yellow") +
  labs(color = "NES") +
  labs(x = "PC3", y= "PC4")

dev.off()
