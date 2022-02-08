# GBM_TF_analysis

## Versions

R:v3.6.1
GSEA:v4.1.0

## Data sources

### downoaded_data

```
#Download databases
cd downloaded_data
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.chr_patch_hapl_scaff.annotation.gtf.gz
wget http://gtrd.biouml.org/downloads/19.10/chip-seq/Homo_sapiens_meta_clusters.interval.gz
gunzip gencode.v27.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip Homo_sapiens_meta_clusters.interval.gz
wget http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/msigdb.v7.4.entrez.gmt
cd ..
```
downloaded_data also contains:
* a file containing {ensembl_id}\t{gene_name}\n for ensembl v75 in order to convert the glass data to IDs: ensembl_v75_geneidtoname.txt
* a file containing {gene_id}\t{gene_name}\n for gencode v27 in order to convert the stead data to names: gencode.v27_geneidtoname.txt
* a file containing {gene_id}\t{hgnc_symbol}\n from https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart: gencode.v27_geneidtohgnc.txt
* the "IlmnID, CHR_hg38, Start_hg38, End_hg38" columns from https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip, with the data starting on row 9: infinium-methylationepic-v-1-0-b5-manifest-file_extract.txt


### original_data 
Contains stead cohort expression tables, count data, batch-corrected-protein-coding-only expression data and metadata

### glass_data
This folder contains glass data from Synapse:
* gene_tpm_matrix_all_samples.tsv: expression data from https://www.synapse.org/#!Synapse:syn23548220
* beta.merged.tsv: methylation data
* variants_anno_20201109.csv: annotations for point variants
* variants_passgeno_20201109_filtered.csv: point variant data filtered for primary and first recurrents of patients in glass_gbm_idhwt_rt_tmz_local+stead.txt, using:
```
while read line ; do cat glass_data/variants_passgeno_20201109.csv  | grep "${line}-TP" | grep 't' >> glass_data/variants_passgeno_20201109_filtered_temp.csv ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+stead.txt
while read line ; do cat glass_data/variants_passgeno_20201109.csv  | grep "${line}-R1" | grep 't' >> glass_data/variants_passgeno_20201109_filtered_temp.csv ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+stead.txt
#remove duplicate patients - preference for WXS with highest variant numbers over WGS or WXS with fewer variants
cat glass_data/variants_passgeno_20201109_filtered_temp.csv | grep -v -E "TCGA-06-0125-R1-11D-WGS-MA69JO|TCGA-06-0125-R1-11D-WXS-8Q4RKD|GLSS-HK-0003-R1-01D-WGS-R7P485|TCGA-06-0125-TP-01D-WXS-Z67AVH|TCGA-06-0125-TP-01D-WGS-FZT8H0|GLSS-HK-0003-TP-01D-WGS-WAGBN9" > glass_data/variants_passgeno_20201109_filtered.csv
```
* variants_titan_seg_filtered.txt: copy number data from variants_titan_seg (syn23554313) filtered for primary and first recurrents of patients in glass_gbm_idhwt_rt_tmz_local+stead_cna.txt, using:
```
while read line ; do cat glass_data/variants_titan_seg.txt  | grep "${line}-TP" | cut -f 1,2,3,4,5,7 >> glass_data/variants_titan_seg_filtered_temp.txt ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+stead.txt 
while read line ; do cat glass_data/variants_titan_seg.txt  | grep "${line}-R1" | cut -f 1,2,3,4,5,7 >> glass_data/variants_titan_seg_filtered_temp.txt ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+stead.txt 
#remove duplicate patients - preference for WGS over WXS
cat glass_data/variants_titan_seg_filtered_temp.txt | grep -v -E "TCGA-06-0125-TP-01-NB-01D-WXS|TCGA-06-0125-TP-02-NB-01D-WXS|TCGA-06-0125-R1-11-NB-01D-WXS|TCGA-06-0125-R1-02-NB-01D-WXS|GLSS-HK-0003-R1-01-NB-01D-WXS" > glass_data/variants_titan_seg_filtered.txt
```

### gene_sets
Contains the following gene set lists from http://www.gsea-msigdb.org/gsea/downloads.jsp:
c2.cgp.v7.4.symbols.gmt
c2.cp.v7.4.symbols.gmt
c3.tft.v7.4.symbols.gmt
c3.mir.mirdb.v7.4.symbols.gmt
c5.go.bp.v7.4.symbols.gmt
c5.go.mf.v7.4.symbols.gmt

## patient_lists
WARNING: Many of the analyses rely on GLASS primary and recurrents being labelled as TP and R1. This is the case for patients inlcuded in the following lists but any alterations need checking.
gbm_idhwt_rt_tmz_local.txt: Stead patients who's primary and first recurrent are GBM_IDHwt, local first recurrent, recieved rt+tmz.
glass_gbm_idhwt.txt: GLASS patients who's primary and first recurrent are GBM_IDHwt or GBM_IDHunknown, have RNA data available and not in the stead cohort.
glass_gbm_idhwt_rt_tmz_local.txt: GLASS patients who's primary and first recurrent are GBM_IDHwt or GBM_IDHunknown, local first recurrent, recieved rt+tmz, have RNA data available and not in the stead cohort.
glass_gbm_idhwt_rt_tmz_local+stead.txt: As with glass_gbm_idhwt_rt_tmz_local.txt but with stead patients included.
glass_gbm_idhwt_rt_tmz_local+stead_cna.txt: As with glass_gbm_idhwt_rt_tmz_local+stead.txt but only those with GLASS copy number alteration data available.
glass_gbm_idhwt_rt_tmz_local_methylation+rna.txt: As with glass_gbm_idhwt_rt_tmz_local.txt but with stead patients included and filtered for those with methylation data also available. 
mixed_gbm_idhwt_rt_tmz_local_methylation+rna.txt: As with glass_gbm_idhwt_rt_tmz_local_methylation+rna.txt but with stead IDs for stead patients.

## Run DEA
```
qsubsec scripts/run_deseq2.qsubsec
SIG=0.05
SIG=0.01
awk -F" " '{if($7!="NA") print}' deseq2/results.txt | cut -d'"' -f2 | cut -d"." -f1 | tail -n+2 > deseq2/background_filtered.txt
awk -F" " '{if($7!="NA") print}' deseq2_uvd/results_up.txt | cut -d'"' -f2 | cut -d"." -f1 | tail -n+2 > deseq2_uvd/background_filtered_up.txt
awk -F" " '{if($7!="NA") print}' deseq2_uvd/results_down.txt | cut -d'"' -f2 | cut -d"." -f1 | tail -n+2 > deseq2_uvd/background_filtered_down.txt
awk -v sig=${SIG} -F" " '{if($7<sig) print}' deseq2/results.txt | cut -d'"' -f2 | cut -d"." -f1 | tail -n+2 > deseq2/deg_${SIG}.txt
awk -v sig=${SIG} -F" " '{if($7<sig) print}' deseq2_uvd/results_up.txt | cut -d'"' -f2 | cut -d"." -f1 | tail -n+2 > deseq2_uvd/deg_up_${SIG}.txt
awk -v sig=${SIG} -F" " '{if($7<sig) print}' deseq2_uvd/results_down.txt | cut -d'"' -f2 | cut -d"." -f1 | tail -n+2 > deseq2_uvd/deg_down_${SIG}.txt
```
## GO enrichment analysis
http://www.webgestalt.org/

Method of interest: Over-Represenataion Analysis

Significance Level: FDR

Max set size: 1000

Set deseq2/background_go_list.txt as reference gene list
```
Rscript scripts/plot_dotplot.R
```

## Expression inputs for GSEA 
get_foldchange.py and get_foldchange_tss.py were run on a previous smaller cohort to get filtered lists of genes and TSSs, separately for those samples processed with total RNA or mRNA libraries: filtered_genelist_mrna.txt,filtered_genelist_total.txt,filtered_tsslist_mrna.txt,filtered_tsslist_total.txt.
get_foldchange_newrealease.py and get_foldchange_tss_newrelease.py were later run on the current cohort expression table using the same lists of genes and TSSs.
get_foldchange_glass.py is used to generate the glass data inputs.
cat ranks/glass_absolute_log2fc/TCGA-06-0125.rnk | cut -f 1 > ranks/filtered_genelist_glass.txt 
WARNING: get_foldchange_glass.py only works if primary and recurrent barcodes are labelled as TP and R1 which is the case for the patients included, but is not guaranteed for other patients.
The resulting inputs are in the 'ranks' folder.


## Create files for running GSEA
```
mkdir gsea_files
SIZE=1000 #2000, 5000

#Create promotor files:
scripts/ExtractPromoterCoords.pl -g downloaded_data/gencode.v27.chr_patch_hapl_scaff.annotation.gtf -o gsea_files -s ${SIZE} 

#Get the overlaps between the promoter file and the metaclusters (requires up to 120GB RAM)
Rscript scripts/overlaps.R —metaclusters downloaded_data/Homo_sapiens_meta_clusters.interval —promotor gsea_files/gencode.v27.chr_patch_hapl_scaff.annotation_PromotersTSS_${SIZE}.txt —output gsea_files/OverlapTFpeaksAndPromoters_${SIZE}_GTRDv19_10_gencodev27.txt

#Filter overlap for peaks from 2+ experiments
awk -F"\t" '{if ($39>1) print}' gsea_files/OverlapTFpeaksAndPromoters_${SIZE}_GTRDv19_10_gencodev27.txt > gsea_files/Min2Exp_OverlapTFpeaksAndPromoters_${SIZE}_GTRDv19_10_gencodev27.txt

#Extract the TF to gene/tss pairs
awk -F'\t' '{print($33,"\t",$11,"\t",$12,"\t",$5)}' gsea_files/Min2Exp_OverlapTFpeaksAndPromoters_${SIZE}_GTRDv19_10_gencodev27.txt | sort -u > gsea_files/TF_to_Gene_${SIZE}_GTRDv19_10_gencodev27.txt
awk -F'\t' '{SUM=($2+$3)/2; print($33"\t"$1":"SUM)}' gsea_files/Min2Exp_OverlapTFpeaksAndPromoters_1000_GTRDv19_10_gencodev27.txt | sort -u > gsea_files/TF_to_TSS_1000_GTRDv19_10_gencodev27.txt 

#Remove hash’s from gene names
perl -pi -e 's/#//g' gsea_files/TF_to_Gene_${SIZE}_GTRDv19_10_gencodev27.txt 
perl -pi -e 's/#//g' gsea_files/TF_to_TSS_1000_GTRDv19_10_gencodev27.txt 

#Make the gmt geneset using Ensembl ID
scripts/MakeGeneSet_ENS.pl -i gsea_files/TF_to_Gene_${SIZE}_GTRDv19_10_gencodev27.txt -o gsea_files/TFs_ENS_${SIZE}_GTRDv19_10_gencodev27.gmt 
scripts/MakeGeneSet_ENS.pl -i gsea_files/TF_to_TSS_1000_GTRDv19_10_gencodev27.txt -o gsea_files/TSS_TFs_ENS_1000_GTRDv19_10_gencodev27.gmt 

#Get only JARID2 entry for running gsea specifically for JARID2
grep 'JARID2' gsea_files/TFs_ENS_1000_GTRDv19_10_gencodev27.gmt > gsea_files/TFs_ENS_1000_GTRDv19_10_gencodev27_JARID2_only.gmt
grep 'JARID2' gsea_files/TSS_TFs_ENS_1000_GTRDv19_10_gencodev27.gmt > gsea_files/TSS_TFs_ENS_1000_GTRDv19_10_gencodev27_JARID2_only.gmt

```

## Get intermediate files
```
awk -F'\t' '{SUM=($2+$3)/2; print($7"\t"$1":"SUM)}' gsea_files/gencode.v27.chr_patch_hapl_scaff.annotation_PromotersTSS_1000.txt | sort -u > intermediate_files/transcript_to_tss_position.txt 
awk -F'\t' '{SUM=($2+$3)/2; print($5"\t"$1":"SUM)}' gsea_files/gencode.v27.chr_patch_hapl_scaff.annotation_PromotersTSS_1000.txt | sort -u > intermediate_files/gene_to_tss_position.txt 
```


## Run GSEA

```
for p in ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_absolute/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute -s ; done 
for p in ranks/actual_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_actual/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual -s ; done 
for p in ranks/absolute_log2fc_tss/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_absolute_tss/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute_tss -s ; done 
for p in ranks/actual_log2fc_tss/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_actual_tss/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual_tss -s ; done 
for p in ranks/glass_absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_absolute_glass/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute_glass -s ; done 
for p in ranks/glass_actual_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_actual_glass/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual_glass -s ; done 

```

## Process GSEA results

```
#Remove unnecesary files
mkdir temp 
for dir in gsea_outputs/*/* ; do ls ${dir}/* | grep -E -v 'PCGF2|CBX2|CBX7|CBX8|EZH2|CBX6|HSD17B8|SS18/SSX1|NR5A1|APOBEC3B|YBX1|SRSF9|PBXIP1|SIRT1|SUPT16H|TFAM|ARID3B|LHX2|JARID2|report' | xargs -I {} mv {} temp/; done 
for dir in gsea_outputs/*/* ; do rm -r ${dir}/edb ; done 
rm -r temp 

### Combine results and generate tables for normalised enrichment score, p-value and FDR
SIZES="1000 2000 5000"
SETS="outputs_actual outputs_actual_glass outputs_absolute outputs_absolute_glass outputs_actual_tss outputs_absolute_tss"

### Get JARID2 NES or ES for all patients
for SET in $SETS ; do for SIZE in 1000 ; do grep 'JARID2' gsea_outputs/${SET}/*_${SIZE}_*/gsea_report_for_na_*tsv | sed "s/_GTRD_${SIZE}_/\//" | tr '/' '\t' | cut -f 3,9,10 | awk '{ if ( $3=="---" ) {print $1"\t"$2} else {print $1"\t"$3} }' > reports/jarid2_results/${SET}_${SIZE}_JARID2_results.tsv ; done ; done
#manually create reports/jarid2_results/outputs_actual_glass_1000_JARID2_results+stead.tsv to include stead patient results as GLASS ids.


### Get tables for all patients and genes
cat gsea_files/TFs_ENS_5000_GTRDv19_10_gencodev27.gmt | cut -f 1 > reports/all_tfs.txt
for SET in $SETS ; do for SIZE in $SIZES ; do mkdir reports/${SET}_${SIZE} ; for f in gsea_outputs/${SET}/*_${SIZE}_* ; do fi=$(basename ${f%_GTRD*}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '=='> reports/${SET}_${SIZE}/${fi}_table.txt ; done ; done ; done

Rscript scripts/gsea_violin_plots.py
```

## cell lines
```
python script/process_cell_lines.py
python script/process_cell_lines_control.py
qsubsec scripts/gsea_cell_lines.qsubsec SIZE=1000 MODE=absolute PATIENT=A172_0,A172_1,A172_2,GBM63_0,GBM63_1,GBM63_2,A172_control
qsubsec scripts/gsea_cell_lines.qsubsec SIZE=1000 MODE=absolute PATIENT=PDspheroids_1,PDspheroids_3
qsubsec scripts/gsea_cell_lines.qsubsec SIZE=1000 MODE=absolute PATIENT=ME_48,ME_72


Rscript scripts/plot_pca_cell_lines_log2fc_project.R

TABLE=primary
TABLE=recurrent
TABLE=log2fc


python scripts/merge.py --input tables/${TABLE}_all.txt,cell_lines/tables/${TABLE}_ME.txt,cell_lines/tables/${TABLE}_A172.txt,cell_lines/tables/${TABLE}_GBM63.txt,cell_lines/tables/${TABLE}_GBM63_CUTRUN.txt,cell_lines/tables/${TABLE}_PDspheroids.txt --output cell_lines/tables/${TABLE}_merged.txt

Rscript scripts/pca.R --patients patient_lists/gbm_idhwt_rt_tmz_local_+_cell_lines.txt --table cell_lines/tables/${TABLE}_merged.txt --categories cell_lines/pca/categories.txt --colour reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv --scale=TRUE --name cell_lines_${TABLE} --label A172_0,A172_1,A172_2,GBM63_0,GBM63_1,GBM63_2,PDspheroids_1,PDspheroids_3,GBM63_CUTRUN_0,GBM63_CUTRUN_1,ME_48,ME_72

Rscript scripts/heatmap_cell_lines.R --patients patient_lists/gbm_idhwt_rt_tmz_local_+_cell_lines.txt --genes analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt --table cell_lines/tables/log2fc_merged.txt --nes_colour reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv --rna_colour original_data/library_types.txt --name cell_lines_le70
```


## Analysis

### Get LE50 and LE70 genes
```
LIST="gbm_idhwt_rt_tmz_local"
LIST="glass_gbm_idhwt_rt_tmz_local"


SET="actual"
SET="actual_tss"
SET="actual_glass"

SIZE="1000"

#get le counts
while read line ; do f=${line}; cat ./gsea_outputs/outputs_${SET}/*${f}*${SIZE}*/JARID2.tsv | grep 'Yes' ; done <patient_lists/${LIST}.txt | cut -f 2 | sort | uniq -c > analysis/leading_edge/counts_${SET}_${SIZE}_${LIST}.txt
#NOTE: This resulted in also counting genes in the Walton50,55,59 outputs when intending to look at genes in the Walton5 output. If rerunning, instead use the below code to fix:
#Run for tss:
#while read line ; do f=${line}; cat ./gsea_outputs/outputs_${SET}/*${f}_*${SIZE}*/JARID2.tsv | grep 'Yes' ; done <patient_lists/${LIST}.txt | cut -f 2 | sort | uniq -c > analysis/leading_edge/counts_${SET}_${SIZE}_${LIST}.txt
 


#calculate minimum number of patient leading edges needed for a gene to be classed as le50 of le70
LE=50
LE=70
TOTAL=$(wc -l patient_lists/${LIST}.txt|cut -d' ' -f1)
NUM=$(awk -v total=$TOTAL -v le=0.$LE 'BEGIN {printf "%.0f", total*le}') ;awk -v  num=$NUM -F" " '{if ($1>=num) print}' analysis/leading_edge/counts_${SET}_${SIZE}_${LIST}.txt | awk '{print $(NF)}' > analysis/leading_edge/le${LE}_${SET}_${SIZE}_${LIST}.txt 

#methylation
while read line ; do f=${line}; cat ./gsea_outputs/outputs_actual_glass/*${f}_*1000*/JARID2.tsv ./gsea_outputs/outputs_actual/*${f}_*1000*/JARID2.tsv | grep 'Yes' ; done <patient_lists/gbm_idhwt_rt_tmz_local_methylation+rna.txt | cut -f 2 | sort | uniq -c > analysis/leading_edge/counts_methylation_1000_gbm_idhwt_rt_tmz_local_methylation+rna.txt
LE=50
LE=70
TOTAL=$(wc -l patient_lists/mixed_gbm_idhwt_rt_tmz_local_methylation+rna.txt|cut -d' ' -f1)
NUM=$(awk -v total=$TOTAL -v le=0.$LE 'BEGIN {printf "%.0f", total*le}') ;awk -v  num=$NUM -F" " '{if ($1>=num) print}' analysis/leading_edge/counts_methylation_1000_mixed_gbm_idhwt_rt_tmz_local_methylation+rna.txt | awk '{print $(NF)}' > analysis/leading_edge/le${LE}_methylation_1000_mixed_gbm_idhwt_rt_tmz_local_methylation+rna.txt 
```

### Sankey plots
Data is in analysis/sankey/
https://www.sankeymatic.com/build/ online generator used to create the plots.

### Get filtered expression tables.
Get filtered expression tables for mRNA, total_RNA, and all, containing the list of genes used for GSEA input. 
All contains all patients but with only the filtered mRNA genes.
```
python scripts/get_foldchange_tables.py
#GLASS tables already made with scripts/get_foldchange_glass.py.
```

### Get gene lists
```
grep 'JARID2' gsea_files/TFs_ENS_1000_GTRDv19_10_gencodev27.gmt | sed 's/\t/\n/g' | cut -d ' ' -f 2| tail -n+3 > gene_lists/JARID2_bound_genes.txt 
cut -f 1 original_data/PvR_genefpkm_all_LS_23062021.txt.txt | cut -d'.' -f1 | tail -n+2 > gene_lists/all_genes.txt
```

### Run PCA

```
PATIENTS=gbm_idhwt_rt_tmz_local

TABLE=log2fc_all
TABLE=primary_all
TABLE=recurrent_all

SCALE=TRUE #also centers 
SCALE=FALSE

Rscript scripts/pca.R --patients patient_lists/${PATIENTS}.txt --genes gene_lists/JARID2_bound_genes.txt --table tables/${TABLE}.txt --colour reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv --scale=${SCALE} --name ${PATIENTS}_JARID2_${TABLE}_${SCALE}
Rscript scripts/pca.R --patients patient_lists/${PATIENTS}.txt --table tables/${TABLE}.txt --colour reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv --scale=${SCALE} --name ${PATIENTS}_${TABLE}_${SCALE}

sort analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings.txt -g -k2  | head -n101 | tail -n+2 | cut -f1 -d" " > analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings_head100.txt
sort analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings.txt -g -k2  | tail -n100 | cut -f1 -d" " > analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings_tail100.txt
cat analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings_head100.txt analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings_tail100.txt > analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings_headtail100.txt 
sort analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings.txt -g -k2  | head -n1001 | tail -n+2 | cut -f1 -d" " > analysis/pca/pca_${PATIENTS}_${TABLE}_${SCALE}_PC1_loadings_head1000.txt
```

### Plot heatmaps
```
Rscript scripts/heatmap.R --patients patient_lists/gbm_idhwt_rt_tmz_local.txt --genes analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt --table tables/log2fc_all.txt --nes_colour reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv --rna_colour original_data/library_types.txt --name gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_le70
Rscript scripts/heatmap.R --patients patient_lists/gbm_idhwt_rt_tmz_local.txt --genes analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_PC1_loadings_headtail100.txt --table tables/log2fc_all.txt --nes_colour reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv --rna_colour original_data/library_types.txt --jarid2 gene_lists/JARID2_bound_genes.txt --name gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_PC1_loadings_headtail100
Rscript scripts/heatmap.R --patients patient_lists/gbm_idhwt_rt_tmz_local.txt --genes analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_PC1_loadings_head1000.txt --table tables/log2fc_all.txt --nes_colour reports/jarid2_results/outputs_actual_1000_JARID2_results.tsv --rna_colour original_data/library_types.txt --jarid2 gene_lists/JARID2_bound_genes.txt --name gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_PC1_loadings_head1000
Rscript scripts/heatmap.R --patients patient_lists/glass_gbm_idhwt_rt_tmz_local.txt --genes analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt --table tables/log2fc_glass.txt  --rna_colour original_data/library_types.txt --nes_colour reports/jarid2_results/outputs_actual_glass_1000_JARID2_results.tsv --name glass_gbm_idhwt_rt_tmz_local_le70


```

## Methylation
```
python scripts/get_probe_promotor_overlap.py
cat gene_lists/JARID2_bound_genes.txt | while read line ; do printf "\n${line}\t" >> methylation/all_JARID2_bound_genes_probes.txt ; grep ${line} methylation/probe_promotor_overlap.txt | cut -f4 | tr "\n" "\t" >> methylation/all_JARID2_bound_genes_probes.txt ; done
cat analysis/leading_edge/le50_methylation_1000_gbm_idhwt_rt_tmz_local_methylation+rna.txt | while read line ; do printf "\n${line}\t" >> methylation/le50_JARID2_bound_genes_probes.txt ; grep ${line} methylation/probe_promotor_overlap.txt | cut -f4 | tr "\n" "\t" >> methylation/le50_JARID2_bound_genes_probes.txt ; done
cat analysis/leading_edge/le70_methylation_1000_gbm_idhwt_rt_tmz_local_methylation+rna.txt | while read line ; do printf "\n${line}\t" >> methylation/le70_JARID2_bound_genes_probes.txt ; grep ${line} methylation/probe_promotor_overlap.txt | cut -f4 | tr "\n" "\t" >> methylation/le70_JARID2_bound_genes_probes.txt ; done
cat gene_lists/all_genes.txt | while read line ; do printf "\n${line}\t" >> methylation/all_genes_probes.txt ; grep ${line} methylation/probe_promotor_overlap.txt | cut -f4 | tr "\n" "\t" >> methylation/all_genes_probes.txt ; done
qsubsec scripts/scripts/get_methylation_values.qsubsec
#WARNING: the get_methylation_values.py script called by the above code only works if primary and recurrent barcodes are labelled as TP and R1 which is the case for the patients included, but is not guaranteed for other patients.

```

## Run UvD DEA
```
qsubsec scripts/run_deseq2_uvd.qsubsec
#plot results
python scripts/plot_deseq2_ranks.py
#get ranks based on both -log10(pvalue) and sign*-log10(pvalue) 
#(duplicate gene values are calculated as the mean and therefore have diferent absolute values between these two methods, ie (-2+3)/2 vs. (2+3)/2) 
python scripts/get_deseq2_ranks.py

```

## Run UvD GSEA
```
cp TFs_ENS_1000_GTRDv19_10_gencodev27.gmt gene_sets/TFs_ENS_1000_GTRDv19_10_gencodev27.gmt
SETS="TFs_ENS_1000_GTRDv19_10_gencodev27 custom_gbm_gene_sets" 
for s in $SETS; do qsubsec scripts/gsea_uvd.qsubsec RUN=run1 GMT=$s GENES=IDs; done
SETS="h.all.v7.4.symbols c2.cgp.v7.4.symbols c2.cp.v7.4.symbols c3.mir.mirdb.v7.4.symbols c3.tft.v7.4.symbols c5.go.bp.v7.4.symbols c5.go.mf.v7.4.symbols" 
for s in $SETS; do qsubsec scripts/gsea_uvd.qsubsec RUN=run1 GMT=$s GENES=symbols; done

SETS="TFs_ENS_1000_GTRDv19_10_gencodev27 custom_gbm_gene_sets h.all.v7.4.symbols c2.cgp.v7.4.symbols c2.cp.v7.4.symbols c3.mir.mirdb.v7.4.symbols c3.tft.v7.4.symbols c5.go.bp.v7.4.symbols c5.go.mf.v7.4.symbols" 
for f in gsea_uvd_outputs/run1/* ; do fi=$(basename ${f}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '==> ' | grep '[0123456789]'> gsea_uvd_outputs/reports/${fi}_table.txt ; done


#Process results, either full or fast method.
Rscript scripts/process_uvd_deseq2_results.R --gmt custom_gbm_gene_sets
for s in $SETS; do qsubsec scripts/process_uvd_deseq2_results.qsubsec SET=$s ; done

#Plot results in a network
#conda activate r4
#source activate
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultscustom_gbm_gene_sets_fast.txt --sets gsea_uvd_outputs/reports/custom_gbm_gene_sets_run1_fdr_0.25.txt --name custom_gbm_gene_sets_run1
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsh.all.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/h.all.v7.4.symbols_run1_fdr_top10_up.txt --name h.all.v7.4.symbols_run1_up
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsh.all.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/h.all.v7.4.symbols_run1_fdr_top10_down.txt --name h.all.v7.4.symbols_run1_down
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsc2.cgp.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/c2.cgp.v7.4.symbols_run1_fdr_1_up.txt --name c2.cgp.v7.4.symbols_run1_up
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsc2.cgp.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/c2.cgp.v7.4.symbols_run1_fdr_1_down.txt --name c2.cgp.v7.4.symbols_run1_down
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsc2.cp.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/c2.cp.v7.4.symbols_run1_top10.txt --name c2.cp.v7.4.symbols_run1
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsc3.mir.mirdb.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/c3.mir.mirdb.v7.4.symbols_run1_fdr_1.txt --name c3.mir.mirdb.v7.4.symbols_run1
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsc3.tft.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/c3.tft.v7.4.symbols_run1_fdr_top5s.txt --name c3.tft.v7.4.symbols_run1
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsc5.go.bp.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/c5.go.bp.v7.4.symbols_run1_fdr_1.txt --name c5.go.bp.v7.4.symbols_run1
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsc5.go.mf.v7.4.symbols_fast.txt --sets gsea_uvd_outputs/reports/c5.go.mf.v7.4.symbols_run1_fdr_top10s.txt --name c5.go.mf.v7.4.symbols_run1
Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsTFs_ENS_1000_GTRDv19_10_gencodev27_fast.txt --sets gsea_uvd_outputs/reports/TFs_ENS_1000_GTRDv19_10_gencodev27_run1_fdr_0.25.txt --name TFs_ENS_1000_GTRDv19_10_gencodev_run1

```
## Variants
```
# Glass copy number

Rscript scripts/plot_cn_freq.R

# Glass point variants

```
#filter for deleterious
grep -E 'DE_NOVO_START_IN_FRAME|DE_NOVO_START_OUT_FRAME|FRAME_SHIFT_DEL|FRAME_SHIFT_INS|IN_FRAME_DEL|IN_FRAME_INS|MISSENSE|NONSENSE|NONSTOP|START_CODON_DEL|START_CODON_INS|START_CODON_SNP' glass_data/variants_anno_20201109.csv > glass_data/variants_anno_20201109_deleterious.csv
grep 'SILENT' glass_data/variants_anno_20201109.csv > glass_data/variants_anno_20201109_silent.csv
```
