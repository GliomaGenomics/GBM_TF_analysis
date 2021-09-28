library(ggplot2)
library(argparse)
library("ComplexHeatmap")

parser <- ArgumentParser()
parser$add_argument("--patients", help="file containing list of patients to include")
parser$add_argument("--genes", help="list of genes to filter on")
parser$add_argument("--table", help="table containing data")
parser$add_argument("--nes_colour", help="NES values to colour patients by")
parser$add_argument("--rna_colour", help="rna library type for patients")
parser$add_argument("--name", help="name for analysis outputs")
args <- parser$parse_args()

datain<-read.table(args$table, header=TRUE, row.names=1)
patients<-scan(args$patients, what = character())
datain=datain[ , as.vector(patients) ]

if( length(args$genes)!=0 )
{
genes<-scan(args$genes, what = character())
datain<-datain[row.names(datain) %in% noquote(genes),]
}

nes<-read.delim(args$nes_colour, col.names=c('patient_id','value'), header=FALSE)
nes<-nes[match(as.vector(patients), nes$patient_id),]
nes<-as.vector(nes$value)
rna<-read.delim(args$rna_colour, col.names=c('patient_id','value'), header=FALSE)
rna<-rna[match(as.vector(patients), rna$patient_id),]
rna<-as.vector(rna$value)

dim(as.matrix(datain))
length(nes)
length(rna)
col = list(Library_type = c("Unstranded_mRNA" = "maroon2", "Stranded_mRNA" = "orchid1", "Stranded_Total" = "green"), JARID2_NES = circlize::colorRamp2(c(-2, 2), c("blue", "yellow")) )
column_ha = HeatmapAnnotation(JARID2_NES = nes, Library_type = rna, col=col)
pdf(paste("analysis/heatmap/heatmap_",args$name,".pdf", sep=""))
Heatmap(as.matrix(datain), name = "Log2FC", bottom_annotation = column_ha, col = circlize::colorRamp2(c(-1, 0, 1), c("royalblue3", "white", "red")))
dev.off()
