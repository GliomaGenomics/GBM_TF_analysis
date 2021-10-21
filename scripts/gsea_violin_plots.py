# Created with code provided by Joe Wilkinson

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(tidyr)

# import NES and FDR datasets by reading .csv files :
patients_stead<-scan("patient_lists/gbm_idhwt_rt_tmz_local.txt", what = character())
patients_glass<-scan("patient_lists/glass_gbm_idhwt_rt_tmz_local.txt", what = character())
patients_stead<-append(patients_stead,"Gene")
patients_glass<-append(patients_glass,"Gene")
NES_data_stead <- read.csv("reports/combined/outputs_absolute_1000_nes.txt", stringsAsFactors = T, sep='\t',na.strings=c("---"))[, as.vector(patients_stead), ]
NES_data_glass <- read.csv("reports/combined/outputs_absolute_glass_1000_nes.txt", stringsAsFactors = T, sep='\t',na.strings=c("---"))[, as.vector(patients_glass), ]
FDR_data_stead <- read.csv("reports/combined/outputs_absolute_1000_fdr.txt", stringsAsFactors = T, sep='\t',na.strings=c("---"))[, as.vector(patients_stead), ]
FDR_data_glass <- read.csv("reports/combined/outputs_absolute_glass_1000_fdr.txt", stringsAsFactors = T, sep='\t',na.strings=c("---"))[, as.vector(patients_glass), ]
means<-rowMeans(apply(subset(NES_data_stead, select = -Gene),2, as.numeric))
NES_data_stead['mean']<-means
FDR_data_stead['mean']<-means
gen_ord<-NES_data_stead[order(-NES_data_stead$mean), ][1:5,]['Gene']
gen_ord<-as.character(apply(gen_ord["Gene"], 2,as.character))
NES_data<-merge(NES_data_stead,NES_data_glass,by="Gene")
FDR_data<-merge(FDR_data_stead,FDR_data_glass,by="Gene")
NES_data<-NES_data[order(-NES_data$mean), ][1:5,]
FDR_data<-FDR_data[order(-FDR_data$mean), ][1:5,]
NES_data<-subset(NES_data, select = -mean)
FDR_data<-subset(FDR_data, select = -mean)
NES_data<-pivot_longer(NES_data, cols = -Gene, names_to = "Cohort", values_to = "NES")
FDR_data<-pivot_longer(FDR_data, cols = -Gene, names_to = "Cohort", values_to = "FDR")
NES_data$Cohort <- ifelse(grepl('\\.', NES_data$Cohort), 'GLASS', 'Stead')
FDR_data$Cohort <- ifelse(grepl('\\.', FDR_data$Cohort), 'GLASS', 'Stead')
NES_data<-as.data.frame(NES_data)
FDR_data<-as.data.frame(FDR_data)

#additional code for function geom_split_violin() which creates split violin plots. 

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

warnings()

# End of code for function : geom_split_violin().
NES_data<-rename(NES_data, Molecule = Gene)
NES_data<-NES_data[ , c('NES', 'Molecule', 'Cohort')] 
FDR_data<-rename(FDR_data, Molecule = Gene)
FDR_data<-FDR_data[ , c('FDR', 'Molecule', 'Cohort')] 

print(NES_data)
#NES_data<-NES_data[order(NES_data$Cohort),]
#FDR_data<-FDR_data[order(FDR_data$Cohort),]

NES_data$NES<-as.double(as.character(NES_data$NES))
FDR_data$FDR<-as.double(as.character(FDR_data$FDR))

NES_data$Molecule<-as.factor(NES_data$Molecule)
FDR_data$Molecule<-as.factor(FDR_data$Molecule)
NES_data$Cohort<-as.factor(NES_data$Cohort)
FDR_data$Cohort<-as.factor(FDR_data$Cohort)
NES_data$Cohort<- relevel(NES_data$Cohort, "Stead")
FDR_data$Cohort<- relevel(FDR_data$Cohort, "Stead")
NES_data$Molecule = factor(NES_data$Molecule, levels=gen_ord)
FDR_data$Molecule = factor(FDR_data$Molecule, levels=gen_ord)
#print(typeof(NES_data$Molecule))
#print(typeof(NES_data$NES))
#print(typeof(NES_data$Cohort))
#print(NES_data$Molecule)
#print(NES_data$NES)
#print(NES_data$Cohort)
print(NES_data$NES)


pdf(paste("analysis/gsea_results/gsea.pdf", sep=""))

# NES PLOT STARTS HERE

p1_NES <- ggplot(NES_data, aes(Molecule, NES, fill=Cohort)) + geom_split_violin(scale = "width") 

# Use bw theme, change font size, change legend position to top right

p2_NES <- p1_NES + theme_bw(base_size = 20) + theme(legend.position = c(0.9, 0.85)) + theme(axis.title.x = element_blank())


# FDR PLOT STARTS HERE


p1_FDR <- ggplot(FDR_data, aes(Molecule, FDR, fill=Cohort)) + geom_split_violin(scale = "width", ylim=c(0,2), show.legend = FALSE)


# Use bw theme, change font size, remove x axis title

p2_FDR <- p1_FDR + theme_bw(base_size = 20)  + theme(axis.title.x = element_blank())

# combine plots vertically

grid.arrange(p2_NES, p2_FDR, nrow = 2)

dev.off()

warnings()



