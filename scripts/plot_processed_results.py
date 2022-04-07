library(enrichplot)
library(ggplot2)
library("cowplot")
library(argparse)

parser <- ArgumentParser(description='plot gsea results')
parser$add_argument('--processed', dest='pro', type='character', help='processed deseq2 results')
parser$add_argument('--sets', dest='sets', type='character', help='file with list of gene sets to plot')
parser$add_argument('--name', dest='name', type='character', help='output name')
parser$add_argument('--seed', dest='seed', type='character', help='seed for plot')
args <- parser$parse_args()


#read in processed
res<-read.table(args$pro,header=TRUE,row.names=1,sep='\t')

case<-row.names(res)
names(case)<-toupper(row.names(res))


deaup1<-read.table('deseq2_uvd/usign_ranks.rnk', header=FALSE)
deaup2<-read.table('deseq2_uvd/usign_ranks_ens.rnk', header=FALSE)
up<-rbind(deaup1,deaup2)
deadown1<-read.table('deseq2_uvd/dsign_ranks.rnk', header=FALSE)
deadown2<-read.table('deseq2_uvd/dsign_ranks_ens.rnk', header=FALSE)
down<-rbind(deadown1,deadown2)

up$V2[is.na(up$V2)] <- 0
down$V2[is.na(down$V2)] <- 0
geneListup<-up$V2
geneListdown<-down$V2
names(geneListup)<-up$V1
names(geneListdown)<-down$V1


gs<-scan(args$sets, what = character())
gens=list()
for(i in gs){
gens[case[i]]<-list(strsplit(res[case[i],"genes"],", ")[[1]])
}

set.seed(args$seed)

p1 <- cnetplot(gens, layout='kk',foldChange=geneListup,cex_label_category=1, cex_gene=1, color_category='#686868', node_label='category', colorEdge=FALSE, showCategory=50)+
scale_color_gradient2(name='-log10(pval)direction', low='blue', mid='gray95',high='red',limits = c(-2, 2),oob = scales::squish)+
theme(legend.position = "none")
p2 <- cnetplot(gens, layout='kk',foldChange=geneListdown,cex_label_category=1, cex_gene=1,color_category='#686868', node_label='category', colorEdge=FALSE, showCategory=50)+
scale_color_gradient2(name='-log10(pval)direction', low='blue', mid='gray95',high='red',limits = c(-2, 2),oob = scales::squish)+
theme(legend.text = element_text(size=15),legend.title = element_text(size=15),legend.position =c(.9, .1))

##set all nodes the same size (double vector)
p1[["layers"]][[2]][["data"]][["size"]]<-rep(c(100),each=length(p1[["layers"]][[2]][["data"]][["size"]]))
p2[["layers"]][[2]][["data"]][["size"]]<-rep(c(100),each=length(p2[["layers"]][[2]][["data"]][["size"]]))
p1[["layers"]][[2]][["show.legend"]]<-FALSE
p2[["layers"]][[2]][["show.legend"]]<-FALSE
##set nodes to be coloured according to NES (character vector)
#p1[["layers"]][[2]][["data"]][["color"]]

pdf(paste("deseq2_uvd/",args$name,".pdf", sep=""), width=20, height=10)
cowplot::plot_grid(p1, p2, ncol=2, labels=c('Up','Down'))

dev.off()



