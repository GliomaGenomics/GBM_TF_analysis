library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--patients", help="file containing list of patients to include")
parser$add_argument("--genes", help="list of genes to filter on")
parser$add_argument("--table", help="table containing data")
parser$add_argument("--colour", help="values to scale colour points by")
parser$add_argument("--name", help="name for analysis outputs")
parser$add_argument("--scale", help="scale pca: TRUE or FALSE")
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
meta<-read.delim(args$colour, col.names=c('patient_id','value'), header=FALSE)
meta<-meta[match(as.vector(patients), meta$patient_id),]


if( length(args$genes)!=0 )
{
genes<-scan(args$genes, what = character())
filter<-intersect(genes,colnames(datain))
datain=datain[ , filter]
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

ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,1])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,1], y=pca$x[,2])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,3], y=pca$x[,4])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow")
ggplot(as.data.frame(pca$x), aes(x=pca$x[,5], y=pca$x[,6])) +
  theme_classic() +
  geom_point(aes(color = meta$value), size = 2.88, alpha=5/6) +
  scale_colour_gradient(low = "blue", high = "yellow")

write.table(pca$rotation[,1],file=paste("analysis/pca/pca_",args$name,"_PC1_loadings.txt", sep=""))
x<-pca$rotation[,1] 
qqnorm(x)
qqline(x, col = "red")   




dev.off()
