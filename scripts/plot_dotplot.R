library(ggplot2)

datain<-read.table('deseq2/enrichment_results_wg_result1638348378_deg_0.01_maxsize_1000_top10.txt', header=TRUE, sep='\t')
pdf(paste("deseq2/webgestalt_0.01_1000.pdf", sep=""), width=5, height=4)

ord<-datain$description[order(datain$enrichmentRatio)]
ord<-as.character(ord)
datain$description <- factor(datain$description,levels=ord)

ggplot(datain, aes(enrichmentRatio, description, colour=FDR, size=size))+ 
           geom_point()+ 
           scale_color_gradient2(name='FDR', low='yellow', mid='green4', high='blue',midpoint=0.09)+
           scale_size_area(name='Set size')+
           scale_size(range=c(3,10))+
           theme(legend.key.size = unit(10, 'points'),legend.title = element_text(size=10), legend.text = element_text(size=8),axis.title.y=element_blank())+
           xlab('Enrichment ratio')+
           scale_x_continuous(breaks=c(1,3,5,7,9), limits=c(1,9)) 

dev.off()

datain<-read.table('analysis/leading_edge/enrichment_results_wg_result1644172684_bioprocesses_0.05_max1000_le70_vs_jbs.txt', header=TRUE, sep='\t')
pdf(paste("analysis/leading_edge/enrichment_bioprocesses_0.05_max1000_le70_vs_jbs.pdf", sep=""), width=6, height=4)

ord<-datain$description[order(datain$enrichmentRatio)]
ord<-as.character(ord)
datain$description <- factor(datain$description,levels=ord)

ggplot(datain, aes(enrichmentRatio, description, colour=FDR, size=size))+ 
           geom_point()+ 
           scale_color_gradient2(name='FDR', low='yellow', mid='green4', high='blue',midpoint=0.025)+
           scale_size_area(name='Set size')+
           scale_size(range=c(3,10))+
           theme(legend.key.size = unit(10, 'points'),legend.title = element_text(size=10), legend.text = element_text(size=8),axis.title.y=element_blank())+
           xlab('Enrichment ratio')+
           scale_x_continuous(breaks=c(0,2,4,6,8,10), limits=c(1,10)) 

dev.off()


datain<-read.table('analysis/leading_edge/enrichment_results_wg_result1644172835_molfunction_0.05_max1000_le70_vs_jbs.txt', header=TRUE, sep='\t')
pdf(paste("analysis/leading_edge/enrichment_molfunction_0.05_max1000_le70_vs_jbs.pdf", sep=""), width=9, height=4)

ord<-datain$description[order(datain$enrichmentRatio)]
ord<-as.character(ord)
datain$description <- factor(datain$description,levels=ord)
ggplot(datain, aes(enrichmentRatio, description, colour=FDR, size=size))+ 
           geom_point()+ 
           scale_color_gradient2(name='FDR', low='yellow', mid='green4', high='blue',midpoint=0.09)+
           scale_size_area(name='Set size')+
           scale_size(range=c(2,5))+
           theme(legend.key.size = unit(10, 'points'),legend.title = element_text(size=10), legend.text = element_text(size=8),axis.title.y=element_blank())+
           xlab('Enrichment ratio')
#           scale_x_continuous(breaks=c(1,3,5,7,9), limits=c(1,9)) 

dev.off()
