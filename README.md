# GBM_TF_analysis

## Versions

R:v3.6.1

## Data sources

### downoaded_data

```
#Download databases
mkdir downloaded_data
cd downloaded_data
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.chr_patch_hapl_scaff.annotation.gtf.gz
wget http://gtrd.biouml.org/downloads/19.10/chip-seq/Homo_sapiens_meta_clusters.interval.gz
gunzip gencode.v27.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip Homo_sapiens_meta_clusters.interval.gz
cd ..
```
downloaded_data also contains a file listing "{Gene ID}\t{Gene name}\n" for each gene in ensembl v75. 

### original_data 
Contains stead cohort expression tables and metadata

### glass_data
This folder contains the gene_tpm_matrix_all_samples.tsv file from the GLASS resources on Synapse: https://www.synapse.org/#!Synapse:syn23548220

patients_gbm_not_in_stead.txt contains a list of 54 patient barcodes for those with IDHwt (or IDH unknown) GBM at both primary and first recurrence, with RNA expression data available, and whose data aren't part of the stead cohort.

Metadata is from clinical_surgeries_100521.tsv

## Expression inputs for GSEA 

get_foldchange.py and get_foldchange_tss.py were run on a previous smaller cohort to get filtered lists of genes and TSSs, separately for those samples processed with total RNA or mRNA libraries: filtered_genelist_mrna.txt,filtered_genelist_total.txt,filtered_tsslist_mrna.txt,filtered_tsslist_total.txt.
get_foldchange_newrealease.py and get_foldchange_tss_newrelease.py were later run on the current cohort expression table using the same lists of genes and TSSs.
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

####################
### Combine results and generate tables for normalised enrichment score, p-value and FDR
SETS="1000_classic_absolute 1000_weighted_actual 2000_classic_absolute 2000_weighted_actual 5000_classic_absolute 5000_weighted_actual" 
#SETS="1000_classic_absolute_tss_stringent 1000_weighted_actual_tss_stringent"
#SETS="1000_classic_absolute_glass,1000_weighted_actual_glass"

### Get NES (or ES if NES not available) scores for JARID2
grep 'JARID2' gsea_outputs/outputs_actual/*_1000_*/gsea_report_for_na_*tsv | sed 's/_GTRD_1000_/\//' | tr '/' '\t' | cut -f 3,9,10 | awk '{ if ( $3=="---" ) {print $1"\t"$2} else {print $1"\t"$3} }' > reports/outputs_actual/JARID2_results.tsv
grep 'JARID2' gsea_outputs/outputs_actual_glass/*_1000_*/gsea_report_for_na_*tsv | sed 's/_GTRD_1000_/\//' | tr '/' '\t' | cut -f 3,9,10 | awk '{ if ( $3=="---" ) {print $1"\t"$2} else {print $1"\t"$3} }' > reports/outputs_actual_glass/JARID2_results.tsv
```
#######################

## Analysis

### Get LE50 and LE70 genes
```
LIST="gbm_idhwt_rt_tmz_local"

SET="actual"
SET="actual_tss"

SIZE="1000"

#get le counts
while read line ; do f=${line}; cat ./gsea_outputs/outputs_${SET}/*${f}*${SIZE}*/JARID2.tsv | grep 'Yes' ; done <patient_lists/${LIST}.txt | cut -f 2 | sort | uniq -c > analysis/leading_edge/counts_${SET}_${SIZE}_${LIST}.txt

#calculate minimum number of patient leading edges needed for a gene to be classed as le50 of le70
LE=50
LE=70
TOTAL=$(wc -l patient_lists/${LIST}.txt|cut -d' ' -f1)
NUM=$(awk -v total=$TOTAL -v le=0.$LE 'BEGIN {printf "%.0f", total*le}') ;awk -v  num=$NUM -F" " '{if ($1>=num) print}' analysis/leading_edge/counts_${SET}_${SIZE}_${LIST}.txt | awk '{print $(NF)}' > analysis/leading_edge/le${LE}_${SET}_${SIZE}_${LIST}.txt 

```

### Get filtered expression tables.
Get filtered expression tables for mRNA, total_RNA, and all, containing the list of genes used for GSEA input. 
All contains all patients but with only the filtered mRNA genes.
```
python scripts/get_foldchange_tables.py
```

### Get gene lists
```
grep 'JARID2' gsea_files/TFs_ENS_1000_GTRDv19_10_gencodev27.gmt | sed 's/\t/\n/g' | tail -n+3 > gene_lists/JARID2_bound_genes.txt 
```

### Run PCA

```
PATIENTS=gbm_idhwt_rt_tmz_local

TABLE=log2fc_all
TABLE=log2fc_primary
TABLE=log2fc_recurrent

SCALE=TRUE
SCALE=FALSE

Rscript scripts/pca.R --patients patient_lists/${PATIENTS}.txt --genes gene_lists/JARID2_bound_genes.txt --table tables/${TABLE}.txt --colour reports/outputs_actual/JARID2_results.tsv --scale=${SCALE} --name ${PATIENTS}_JARID2_${TABLE}_${SCALE}
Rscript scripts/pca.R --patients patient_lists/${PATIENTS}.txt --table tables/${TABLE}.txt --colour reports/outputs_actual/JARID2_results.tsv --scale=${SCALE} --name ${PATIENTS}_${TABLE}_${SCALE}
```
