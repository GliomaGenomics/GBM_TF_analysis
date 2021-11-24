
genesets<-list()
file_list <- list.files(path="gene_sets/", pattern="*gmt") 
for (file in file_list){
lines<-readLines(file)
l<-strsplit(lines,"\t")
for (i in l){
genesets[i[1]]<-list(i[3:length(i)])
}
}

deaup1<-read.table(deseq2_uvd/usign_ranks.rnk, header=FALSE, row.names=1)
deaup2<-read.table(deseq2_uvd/u_ranks.rnk, header=FALSE, row.names=1)
deaup<-rbind(deaup1,deaup2)
deadown1<-read.table(deseq2_uvd/dsign_ranks.rnk, header=FALSE, row.names=1)
deadown2<-read.table(deseq2_uvd/d_ranks.rnk, header=FALSE, row.names=1)
deadown<-rbind(deadown1,deadown2)


all_resultsup<-data.frame()
file_list<-list.files("gsea_uvd_outputs/reports/*_up.*table.txt")
for (file in file_list){
res<-read.table(file,header=FALSE)
result=data.frame(ID=as.character(res$V1))
result$setSize=as.character(res$V5)
result$NES_up=as.numeric(res$V7)
result$p.adjust_up=as.numeric(res$V9)
result$collection<-rep(substring(file,25,length(file)-10),nrow(result))
genes=list()
valuesup=list()
valuesdown=list()
for (i in res$V1){
genes<-append(genes,genesets[i])
valsup<-list()
valsdown<-list()
for (ii in genesets[i]({
valsup<-append(valsup,deaup[1,ii])
valsdown<-append(valsdown,deadown[1,ii])
}
valuesup<-append(valuesup,valsup)
valuesdown<-append(valuesdown,valsdown)
}
result$genes<-genes
result$valsup<-valsup
result$valsdown<-valsdown
all_resultsup<-rbind(all_results, result)
} 

all_resultsdown<-data.frame()
file_list<-list.files("gsea_uvd_outputs/reports/*_down.*table.txt")
for (file in file_list){
res<-read.table(file,header=FALSE)
result=data.frame(ID=as.character(res$V1))
result$NES_up=as.numeric(res$V7)
result$p.adjust_up=as.numeric(res$V9)
}
all_resultsdown<-rbind(all_resultsdown, result)

all_results<-merge(all_resultsup,all_resultsdown, by="ID")

write.table(all_results[, ! names(all_results) %in% c("valsup", "valsdown","genes")], file='gsea_uvd_outputs/processed_gsea_results.txt') 


