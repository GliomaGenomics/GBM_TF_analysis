library(ggplot2)
library(argparse)
library(ggrepel)

parser <- ArgumentParser()
parser$add_argument("--patients", help="file containing list of patients to include")
parser$add_argument("--genes", help="list of genes to filter on")
parser$add_argument("--table", help="table containing data")
parser$add_argument("--colour", help="values to scale colour points by")
parser$add_argument("--categories", help="categories to colour points by")
parser$add_argument("--name", help="name for analysis outputs")
parser$add_argument("--scale", help="scale pca: TRUE or FALSE")
parser$add_argument("--label", help="points to be labelled")
args <- parser$parse_args()

if (args$scale=="TRUE")
{
args$scale<-TRUE
} else {
args$scale<-FALSE
}


patients<-scan(args$patients, what = character())
datain<-read.table(args$table, header=TRUE, row.names=1)
datain=datain[ , as.vector(patients) ]
datain=t(datain)
cn<-colnames(datain)
colnames(datain)<-substr(cn,1,15)


if( length(args$genes)!=0 )
{
genes<-scan(args$genes, what = character())
filter<-intersect(genes,colnames(datain))
datain=datain[ , filter]
}

if( length(args$label)!=0 )
{
labels<-strsplit(args$label, split=',')[[1]]
}else{
labels<-c()
}

pca<- prcomp(datain, center = args$scale, scale= args$scale)

save(pca,file=(paste("analysis/pca/pca_",args$name,".RData", sep="")))

pdf(paste("analysis/pca/pca_",args$name,".pdf", sep=""))
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


if( length(args$colour)!=0 )
{
meta<-read.delim(args$colour, col.names=c('patient_id','value'), header=FALSE)
meta<-meta[match(as.vector(patients), meta$patient_id),]

xpc<-1
for (ypc in c(2,3,4,5,6))
{
print(ggplot(as.data.frame(pca$x), aes(x=pca$x[,xpc], y=pca$x[,ypc])) +
  theme_classic() +
  labs(x = paste("PC",xpc,sep=""), y= paste("PC",ypc,sep=""))+
  geom_point(aes(color = meta$value), size = 2.88, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow")+
  geom_text_repel(aes(label=ifelse(names(pca$x[,1]) %in% labels,as.character(names(pca$x[,1])),'')),hjust=0,vjust=0, max.time=10, max.overlaps = Inf)+
  labs(color = "NES"))
}
}

if( length(args$categories)!=0 )
{
cat<-read.delim(args$categories, col.names=c('patient_id','value'), header=FALSE)
cat<-cat[match(as.vector(patients), cat$patient_id),]

xpc<-1
for (ypc in c(2,3,4,5,6)){
print(ggplot(as.data.frame(pca$x), aes(x=pca$x[,xpc], y=pca$x[,ypc])) +
  theme_classic() +
  labs(x = paste("PC",xpc,sep=""), y= paste("PC",ypc,sep=""))+
  geom_point(aes(color = cat$value), size = 2.88, alpha=5/6)+
  labs(color = "sample"))

}
}


write.table(pca$rotation[,1],file=paste("analysis/pca/pca_",args$name,"_PC1_loadings.txt", sep=""))
write.table(pca$rotation[,2],file=paste("analysis/pca/pca_",args$name,"_PC2_loadings.txt", sep=""))
write.table(pca$rotation[,3],file=paste("analysis/pca/pca_",args$name,"_PC3_loadings.txt", sep=""))
write.table(pca$rotation[,4],file=paste("analysis/pca/pca_",args$name,"_PC4_loadings.txt", sep=""))
x<-pca$rotation[,1] 
qqnorm(x)
qqline(x, col = "red")   




dev.off()
