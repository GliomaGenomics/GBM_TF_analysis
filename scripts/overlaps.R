library(GenomicRanges)
library(argparse)

parser <- ArgumentParser(description='overlaps')
parser$add_argument('-m', '--metaclusters', dest='metaclusters', type='character', help='A meta_clusters.interval file')
parser$add_argument('-p', '--promotor', dest='promotor', type='character', help='A promotor file')
parser$add_argument('-o', '--output', dest='output', type='character', help='Output file')
args <- parser$parse_args()

data<-read.delim(args$metaclusters)
peaks<-makeGRangesFromDataFrame(data,keep.extra.columns=TRUE,starts.in.df.are.0based=TRUE,ignore.strand=TRUE,seqnames.field="X.CHROM", start.field="START", end.field="END")
rm(data)
proms<-read.delim(args$promotor)
promoters<-makeGRangesFromDataFrame(proms,keep.extra.columns=TRUE)
rm(proms)
merged<-mergeByOverlaps(promoters,peaks,type="any",ignore.strand=TRUE)
write.table(merged,args$output,sep="\t",quote=F,row.names=F)



#data<-read.delim(args$metaclusters)
#proms<-read.delim(args$promotor)
#peaks<-makeGRangesFromDataFrame(data,keep.extra.columns=TRUE,starts.in.df.are.0based=TRUE,ignore.strand=TRUE,seqnames.field="X.CHROM", start.field="START", end.field="END")
#promoters<-makeGRangesFromDataFrame(proms,keep.extra.columns=TRUE)
#merged<-mergeByOverlaps(promoters,peaks,type="any",ignore.strand=TRUE)
#write.table(merged,args$output,sep="\t",quote=F,row.names=F)

