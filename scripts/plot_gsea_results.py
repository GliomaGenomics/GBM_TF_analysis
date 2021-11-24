library(enrichplot)
library(ggplot2)
library("cowplot")
library(argparse)

parser <- ArgumentParser(description='plot gsea results')
parser$add_argument('--gmt', dest='gmt', type='character', help='gmt file')
args <- parser$parse_args()


#res<-read.table(Sys.glob(paste("gsea_uvd_outputs/reports/h*table.txt", sep="")),header=FALSE)
#lines<-readLines(paste("downloaded_data/h.all.v7.4.symbols.gmt",sep=""))


#read in gmt
lines<-readLines(paste("downloaded_data/",args$gmt,".gmt",sep=""))
l<-strsplit(lines,"\t")
geneSets<-list()
for (i in l){
geneSets[i[1]]<-list(i[3:length(i)])
}

##read in gsea results
res<-read.table(Sys.glob(paste("gsea_uvd_outputs/reports/",args$gmt,"*_up_*table.txt", sep="")),header=FALSE)
result=data.frame(ID=as.character(res$V1))
result$Description=as.character(res$V2)
result$setSize=as.character(res$V5)
result$NES=as.numeric(res$V7)
result$pvalue=as.numeric(res$V8)
result$p.adjust=as.numeric(res$V9)
result<-result[result$p.adjust<0.05, ]
result<-result[order(-result$NES,result$p.adjust), ]

gens=list()
for(i in result$ID){
gens[i]<-list(geneSets[[i]])
}

up<-read.table("deseq2_uvd/u_ranks.rnk")
down<-read.table("deseq2_uvd/d_ranks.rnk")
geneListup<-up$V2
geneListdown<-down$V2
names(geneListup)<-up$V1
names(geneListdown)<-down$V1

set.seed(123)

p1 <- cnetplot(gens, foldChange=geneListup,cex_label_category=0.5, cex_gene=0.5, node_label='category', colorEdge=FALSE, showCategory=5)+
scale_color_gradient2(name='associated data', low='blue', mid='white',high='red',limits = c(-3, 3))+
theme(legend.text = element_text(size=5),legend.title = element_text(size=5))
p2 <- cnetplot(gens, foldChange=geneListdown,cex_label_category=0.5, cex_gene=0.5, node_label='category', colorEdge=FALSE, showCategory=5)+
scale_color_gradient2(name='associated data', low='blue', mid='white',high='red',limits = c(-3, 3))+
theme(legend.text = element_text(size=5),legend.title = element_text(size=5))

##set all nodes the same size (double vector)
#p1[["layers"]][[2]][["data"]][["size"]]<-rep(c(100),each=length(p1[["layers"]][[2]][["data"]][["size"]])
#p2[["layers"]][[2]][["data"]][["size"]]<-rep(c(100),each=length(p2[["layers"]][[2]][["data"]][["size"]])
##set nodes to be coloured according to NES (character vector)
#p1[["layers"]][[2]][["data"]][["color"]]

pdf("test.pdf", width=14, height=10)
cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])

dev.off()



