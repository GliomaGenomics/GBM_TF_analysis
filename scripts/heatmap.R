library(ggplot2)
library(argparse)
library("ComplexHeatmap")

parser <- ArgumentParser()
parser$add_argument("--patients", help="file containing list of patients to include")
parser$add_argument("--genes", help="list of genes to filter on")
parser$add_argument("--table", help="table containing data")
parser$add_argument("--nes_colour", help="NES values to colour patients by")
parser$add_argument("--rna_colour", help="rna library type for patients")
parser$add_argument("--jarid2", help="list of JARID2 genes to colour genes by")
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

upvdown<-c()
for(i in nes) {
if(i>0){
upvdown<-append(upvdown,"Up")
} else {
upvdown<-append(upvdown,"Down")
}
}

if( length(args$rna_colour)!=0 ) {
rna<-read.delim(args$rna_colour, col.names=c('patient_id','value'), header=FALSE)
rna<-rna[match(as.vector(patients), rna$patient_id),]
rna<-as.vector(rna$value)
col = list(Library_type = c("Unstranded_mRNA" = "palegreen3", "Stranded_mRNA" = "cyan4", "Stranded_Total" = "darkgreen"), JARID2_NES = circlize::colorRamp2(c(-2, 2), c("blue", "yellow")) )
column_ha = HeatmapAnnotation(JARID2_NES = nes, Library_type = rna, col=col,show_annotation_name =FALSE)
 } else {
col = list(JARID2_NES = circlize::colorRamp2(c(-2, 2), c("blue", "yellow")) )
column_ha = HeatmapAnnotation(JARID2_NES = nes, col=col,show_annotation_name =FALSE)
}

pdf(paste("analysis/heatmap/heatmap_",args$name,".pdf", sep=""))

if( length(args$jarid2)!=0 ) {
    jarid2<-scan(args$jarid2, what = character())
    jlist<- c()
    for (i in row.names(datain)) {
        if(i %in% jarid2){
            jlist<-append(jlist, "Yes")
        } else {
            jlist<-append(jlist, "No")
        }
    }
print(length(jlist))
print(jlist)
row_ha =rowAnnotation(JARID2_gene = jlist, col=list(JARID2_gene = c("No" = "grey80", "Yes" = "grey40")),show_annotation_name = FALSE)
colnames(datain)<-NULL
rownames(datain)<-NULL
ht=Heatmap(as.matrix(datain),  show_column_dend = FALSE, show_row_dend = FALSE,clustering_method_rows = "complete", clustering_method_columns = "complete", column_split = factor(upvdown, level=c('Up','Down')),cluster_column_slices = FALSE, clustering_distance_rows = "pearson", clustering_distance_columns = "pearson", name = "Log2FC", bottom_annotation = column_ha, right_annotation = row_ha, col = circlize::colorRamp2(c(-1, 0, 1), c("steelblue2", "white", "indianred2")), heatmap_legend_param = list(labels = c("<-1", "", "0", "", ">1")),use_raster = TRUE,width = unit(8, "cm"), height = unit(10, "cm"))

} else {
colnames(datain)<-NULL
rownames(datain)<-NULL
ht=Heatmap(as.matrix(datain), show_column_dend = FALSE, show_row_dend = FALSE, clustering_method_rows = "complete", clustering_method_columns = "complete", column_split = factor(upvdown, level=c('Up','Down')), cluster_column_slices = FALSE,clustering_distance_rows = "pearson", clustering_distance_columns = "pearson", name = "Log2FC", bottom_annotation = column_ha, col = circlize::colorRamp2(c(-1, 0, 1), c("steelblue2", "white", "indianred2")), heatmap_legend_param = list(labels = c("<-1", "", "0","", ">1")),use_raster = TRUE,width = unit(8, "cm"), height = unit(10, "cm"))
}

draw(ht,merge_legend = TRUE) #legend_grouping = "original")

dev.off()
