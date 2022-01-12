library(argparse)

parser <- ArgumentParser(description='process deseq2 results')
parser$add_argument('--gmt', dest='gmt', type='character', help='gmt file')
args <- parser$parse_args()

genesets<-list()
lines<-readLines(paste("gene_sets/",args$gmt,'.gmt',sep=""))
l<-strsplit(lines,"\t")
for (i in l){
ii<-as.character(lapply(i, function(x) gsub(" ", "", x)))
genesets[ii[1]]<-list(ii[3:length(ii[ii!= ""])])
}

deaup1<-read.table('deseq2_uvd/usign_ranks.rnk', header=FALSE, row.names=1)
deaup2<-read.table('deseq2_uvd/usign_ranks_ens.rnk', header=FALSE, row.names=1)
deaup<-rbind(deaup1,deaup2)
deadown1<-read.table('deseq2_uvd/dsign_ranks.rnk', header=FALSE, row.names=1)
deadown2<-read.table('deseq2_uvd/dsign_ranks_ens.rnk', header=FALSE, row.names=1)
deadown<-rbind(deadown1,deadown2)


results<-data.frame(matrix(ncol = 2, nrow = 0))
colnames(results)<-c("id","genes")

for (id in names(genesets)){
genes<-genesets[[id]]
genescheck=c()
upinup_genescheck=c()
downinup_genescheck=c()

for (i in genes){
if(!(is.na(deaup[i,]) & is.na(deadown[i,]))){

genescheck<-append(genescheck,i)

if(deaup[i,]>0){
upinup_genescheck<-append(upinup_genescheck, i)
}else{
downinup_genescheck<-append(downinup_genescheck, i)
}
}
}

genescheck<-paste(genescheck, collapse=", ")
upinup_genescheck<-paste(upinup_genescheck, collapse=", ")
downinup_genescheck<-paste(downinup_genescheck, collapse=", ")

results[nrow(results) + 1,] = list(id,genescheck)
results[nrow(results) + 1,] = list(paste(id,'_1',sep=''),upinup_genescheck)
results[nrow(results) + 1,] = list(paste(id,'_2',sep=''),downinup_genescheck)
}

write.table(results, file=paste('deseq2_uvd/processed_results',args$gmt,'_fast.txt',sep=""), quote=FALSE, sep='\t',row.names=FALSE) 

