library(argparse)

parser <- ArgumentParser(description='process deseq2 results')
parser$add_argument('--gmt', dest='gmt', type='character', help='gmt file')
args <- parser$parse_args()

genesets<-list()
lines<-readLines(paste("gene_sets/",args$gmt,'.gmt',sep=""))
l<-strsplit(lines,"\t")
for (i in l){
genesets[i[1]]<-list(i[3:length(i[i!= ""])])
}

deaup1<-read.table('deseq2_uvd/usign_ranks.rnk', header=FALSE, row.names=1)
deaup2<-read.table('deseq2_uvd/usign_ranks_ens.rnk', header=FALSE, row.names=1)
deaup<-rbind(deaup1,deaup2)
deadown1<-read.table('deseq2_uvd/dsign_ranks.rnk', header=FALSE, row.names=1)
deadown2<-read.table('deseq2_uvd/dsign_ranks_ens.rnk', header=FALSE, row.names=1)
deadown<-rbind(deadown1,deadown2)


results<-data.frame(matrix(ncol = 6, nrow = 0))
colnames(results)<-c("id","genes","up_values","down_values","length","proportion_different")

for (id in names(genesets)){
print(id)
genes<-genesets[[id]]
upvalues=c()
upinup_upvalues=c()
downinup_upvalues=c()
upinup_downvalues=c()
downinup_downvalues=c()
downvalues=c()
genescheck=c()
upinup_genescheck=c()
downinup_genescheck=c()
dif<-0
upinup_dif<-0
downinup_dif<-0

for (i in genes){
if(!(is.na(deaup[i,]) & is.na(deadown[i,]))){


upvalues<-append(upvalues,deaup[i,])
downvalues<-append(downvalues,deadown[i,])
genescheck<-append(genescheck,i)
#if((deaup[i,]<1.3 & deadown[i,]>1.3 )|(deaup[i,]>1.3 & deadown[i,]<1.3)|(deaup[i,]< -1.3 & deadown[i,]> -1.3 )|(deaup[i,]> -1.3 & deadown[i,]< -1.3)){
#dif<-dif+1
#}

if(deaup[i,]>0){
upinup_upvalues<-append(upinup_upvalues,deaup[i,])
upinup_downvalues<-append(upinup_downvalues,deadown[i,])
upinup_genescheck<-append(upinup_genescheck, i)
#if((deaup[i,]<1.3 & deadown[i,]>1.3 )|(deaup[i,]>1.3 & deadown[i,]<1.3)|(deaup[i,]< -1.3 & deadown[i,]> -1.3 )|(deaup[i,]> -1.3 & deadown[i,]< -1.3)){
#upinup_dif<-upinup_dif+1
#}

}else{
downinup_upvalues<-append(downinup_upvalues,deaup[i,])
downinup_downvalues<-append(downinup_downvalues,deadown[i,])
downinup_genescheck<-append(downinup_genescheck, i)
#if((deaup[i,]<1.3 & deadown[i,]>1.3 )|(deaup[i,]>1.3 & deadown[i,]<1.3)|(deaup[i,]< -1.3 & deadown[i,]> -1.3 )|(deaup[i,]> -1.3 & deadown[i,]< -1.3)){
#downinup_dif<-downinup_dif+1
#}

}
}
}

len=length(genescheck)
upinup_len=length(upinup_genescheck)
downinup_len=length(downinup_genescheck)
#proportion<-dif/len
proportion<-0
upinup_proportion<-upinup_dif/upinup_len
downinup_proportion<-downinup_dif/downinup_len
genescheck<-paste(genescheck, collapse=", ")
upinup_genescheck<-paste(upinup_genescheck, collapse=", ")
downinup_genescheck<-paste(downinup_genescheck, collapse=", ")
upvalues<-paste(upvalues, collapse=", ")
upinup_upvalues<-paste(upinup_upvalues, collapse=", ")
downinup_upvalues<-paste(downinup_upvalues, collapse=", ")
downvalues<-paste(downvalues, collapse=", ")
upinup_downvalues<-paste(upinup_downvalues, collapse=", ")
downinup_downvalues<-paste(downinup_downvalues, collapse=", ")

results[nrow(results) + 1,] = list(id,genescheck,upvalues,downvalues,len,proportion)
results[nrow(results) + 1,] = list(paste(id,'_1',sep=''),upinup_genescheck,upinup_upvalues,upinup_downvalues,upinup_len,upinup_proportion)
results[nrow(results) + 1,] = list(paste(id,'_2',sep=''),downinup_genescheck,downinup_upvalues,downinup_downvalues,downinup_len,downinup_proportion)
}

write.table(results, file=paste('deseq2_uvd/processed_results',args$gmt,'.txt',sep=""), quote=FALSE, sep='\t',row.names=FALSE) 

