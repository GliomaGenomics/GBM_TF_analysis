library("ggplot2")

patients<-scan('patient_lists/gbm_idhwt_rt_tmz_local_+_cell_lines.txt', what = character())
primdatain<-read.table('cell_lines/tables/primary_merged.txt', header=TRUE)
recudatain<-read.table('cell_lines/tables/recurrent_merged.txt', header=TRUE)
colnames(primdatain) <- paste(colnames(primdatain),"_P",sep="")
colnames(recudatain) <- paste(colnames(recudatain),"_R",sep="")
datain<-merge(primdatain,recudatain, by.x='Rows_P', by.y='Rows_R')

names(datain)[names(datain)=="Rows_P"] <- "EnsID"
dat=data.frame(EnsID=datain$EnsID)
for (p in patients){
dat[paste(p,"_P",sep="")]<-datain[,grep(paste(p,"_P",sep=""), colnames(datain))]
dat[paste(p,"_R",sep="")]<-datain[,grep(paste(p,"_R",sep=""), colnames(datain))]
}
rownames(dat)<-dat$EnsID
dat<-subset(dat, select = -EnsID)
dat<-t(dat)

meta<-read.delim('reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv', col.names=c('patient_id','value'), header=FALSE)
meta<-meta[match(as.vector(patients), meta$patient_id),]

up<-meta[meta$value > 0 ,]$patient_id
down<-meta[meta$value < 0 ,]$patient_id
cell_line<-c('A172_0','A172_1','A172_2','GBM63_0','GBM63_1','GBM63_2','PDspheroids_1','PDspheroids_3','GBM63_CUTRUN_0','GBM63_CUTRUN_1')


sample<-data.frame(matrix(ncol = 2, nrow = 0))
colnames(sample)<-c('patient_id','sample')

for(p in up){
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_P',sep=''),sample="Primary_up"))
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_R',sep=''),sample="Recurrent_up"))
}
for(p in down){
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_P',sep=''),sample="Primary_down"))
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_R',sep=''),sample="Recurrent_down"))
}
for(p in cell_line){
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_P',sep=''),sample="Untreated_cell_line"))
    sample<-rbind(sample,data.frame(patient_id=paste(p,'_R',sep=''),sample="Treated_cell_line"))
}

pca<- prcomp(dat, center = TRUE, scale= TRUE)

cols<-merge(pca$x,sample,by.x=0, by.y='patient_id')

pdf("cell_lines/pca/pca_combined_cell_lines.pdf")

xpc<-1
for (ypc in c(2,3,4,5,6,7,8,9,10)){
print(ggplot(as.data.frame(pca$x), aes(x=pca$x[,xpc], y=pca$x[,ypc], colour=cols$sample)) +
  theme_classic() +
  labs(x = paste("PC",xpc,sep=""), y= paste("PC",ypc,sep=""))+
  geom_point(size = 2.88, alpha=5/6)+
  labs(color = "sample")+
  scale_color_manual(values = c("Primary_up" = "chartreuse2","Recurrent_up"="green4", "Primary_down"="deepskyblue","Recurrent_down"="mediumblue","Untreated_cell_line"="maroon1","Treated_cell_line"="darkviolet")))
}


dev.off()


