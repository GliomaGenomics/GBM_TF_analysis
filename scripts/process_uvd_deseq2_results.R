
library(argparse)
library(qvalue)

parser <- ArgumentParser(description='process deseq2 results')
parser$add_argument('--gmt', dest='gmt', type='character', help='gmt file')
args <- parser$parse_args()

genesets<-list()
lines<-readLines(paste("gene_sets/",args$gmt,'.gmt',sep=""))
l<-strsplit(lines,"\t")
for (i in l){
ii<-gsub(" ", "", i)
genesets[ii[1]]<-list(ii[3:length(ii[ii!= ""])])
}


deaup1<-read.table('deseq2_uvd/usign_ranks.rnk', header=FALSE, row.names=1)
deaup2<-read.table('deseq2_uvd/usign_ranks_ens.rnk', header=FALSE, row.names=1)
deaup<-rbind(deaup1,deaup2)
deadown1<-read.table('deseq2_uvd/dsign_ranks.rnk', header=FALSE, row.names=1)
deadown2<-read.table('deseq2_uvd/dsign_ranks_ens.rnk', header=FALSE, row.names=1)
deadown<-rbind(deadown1,deadown2)

up_gsea<-read.table(Sys.glob(file.path("gsea_uvd_outputs/reports/",paste(args$gmt,'_run1_up.*.txt',sep=""))),header=FALSE,row.names=NULL,fill=TRUE,sep="\t")
up_gsea<-up_gsea[!is.na(up_gsea$V1),]
row.names(up_gsea)<-up_gsea$V1
down_gsea<-read.table(Sys.glob(file.path("gsea_uvd_outputs/reports/",paste(args$gmt,'_run1_down.*.txt',sep=""))),header=FALSE,row.names=NULL,fill=TRUE,sep="\t")
down_gsea<-down_gsea[!is.na(down_gsea$V1),]
row.names(down_gsea)<-down_gsea$V1
results<-data.frame(matrix(ncol = 14, nrow = 0))
colnames(results)<-c("id","gene_set","up_increasing","up_decreasing","up_stable","down_increasing","down_decreasing","down_stable","up_NES","up_FDR","down_NES","down_FDR","chi_squared_pval_2x2","chi_squared_pval_2x3")

for (id in names(genesets)){
print(id)
genes<-genesets[[id]]
upincreasing=c()
updecreasing=c()
upstable=c()
downincreasing=c()
downdecreasing=c()
downstable=c()
for (i in genes){
if(is.na(deaup[i,])){
upstable<-append(upstable,deaup[i,]) 
}else if(deaup[i,]>1.3){
upincreasing<-append(upincreasing,deaup[i,])
}else if(deaup[i,]<(-1.3)){
updecreasing<-append(updecreasing,deaup[i,])
}else{
upstable<-append(upstable,deaup[i,])   
}
if(is.na(deadown[i,])){
downstable<-append(downstable,deadown[i,])   
}else if(deadown[i,]>1.3){
downincreasing<-append(downincreasing,deadown[i,])
}else if(deadown[i,]<(-1.3)){
downdecreasing<-append(downdecreasing,deadown[i,])
}else{
downstable<-append(downstable,deadown[i,])   
}
}
if (sum(c(length(upincreasing),length(updecreasing),length(downincreasing),length(downdecreasing)))!=0){
    pval_2x2<-chisq.test(matrix(c(length(upincreasing),length(updecreasing),length(downincreasing),length(downdecreasing)),ncol=2,byrow=TRUE))$p.value
    pval_2x3<-chisq.test(matrix(c(length(upincreasing),length(updecreasing),length(upstable),length(downincreasing),length(downdecreasing),length(downstable)),ncol=3,byrow=TRUE))$p.value
} else{
    pval_2x2<-1
    pval_2x3<-1
}
results[nrow(results) + 1,] = list(id,length(genes),length(upincreasing),length(updecreasing),length(upstable),length(downincreasing),length(downdecreasing),length(downstable),up_gsea[toupper(id),"V6"],up_gsea[toupper(id),"V8"],down_gsea[toupper(id),"V6"],down_gsea[toupper(id),"V8"],pval_2x2,pval_2x3)
}

results["chi_squared_FDR_2x2"]<-qvalue(results["chi_squared_pval_2x2"],pi0=1)$qvalues
results["chi_squared_FDR_2x3"]<-qvalue(results["chi_squared_pval_2x3"],pi0=1)$qvalues

write.table(results, file=paste('deseq2_uvd/processed_gsea_results',args$gmt,'.txt',sep=""), quote=FALSE, sep='\t',row.names=FALSE) 

