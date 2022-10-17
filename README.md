# GBM_TF_analysis

#conda activate env/r4


## Versions

R:v4.1.0
GSEA:v4.1.0
DESeq2:1.34.0 (originally v1.26.0)

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
* a file containing {gene_id}\t{gene_name}\n for gencode v27 in order to convert the discovery data to names: gencode.v27_geneidtoname.txt
* a file containing {gene_id}\t{hgnc_symbol}\n from https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart: gencode.v27_geneidtohgnc.txt
* the "IlmnID, CHR_hg38, Start_hg38, End_hg38" columns from https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip, with the data starting on row 9: infinium-methylationepic-v-1-0-b5-manifest-file_extract.txt
* hg38.fasta

### original_data 
Contains discovery cohort expression tables, count data, batch-corrected-protein-coding-only expression data and metadata
Also contains JARID2promotersPeaks.bed which is a control where ChIP-seq was performed on a single recurrent sample using an anti-JARID2 antibody and input DNA which had been extracted and sheared but not immunoprecipitated.

### glass_data
This folder contains glass validation data from Synapse:
* gene_tpm_matrix_all_samples.tsv: expression data from https://www.synapse.org/#!Synapse:syn23548220
* beta.merged.tsv: methylation data
* variants_anno_20201109.csv: annotations for point variants
* variants_passgeno_20201109_filtered.csv: point variant data filtered for primary and first recurrents of patients in glass_gbm_idhwt_rt_tmz_local+dis.txt, using:
```
while read line ; do cat glass_data/variants_passgeno_20201109.csv  | grep "${line}-TP" | grep 't' >> glass_data/variants_passgeno_20201109_filtered_temp.csv ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+dis.txt
while read line ; do cat glass_data/variants_passgeno_20201109.csv  | grep "${line}-R1" | grep 't' >> glass_data/variants_passgeno_20201109_filtered_temp.csv ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+dis.txt
#remove duplicate patients - preference for WXS with highest variant numbers over WGS or WXS with fewer variants
cat glass_data/variants_passgeno_20201109_filtered_temp.csv | grep -v -E "TCGA-06-0125-R1-11D-WGS-MA69JO|TCGA-06-0125-R1-11D-WXS-8Q4RKD|GLSS-HK-0003-R1-01D-WGS-R7P485|TCGA-06-0125-TP-01D-WXS-Z67AVH|TCGA-06-0125-TP-01D-WGS-FZT8H0|GLSS-HK-0003-TP-01D-WGS-WAGBN9" > glass_data/variants_passgeno_20201109_filtered.csv
```
* variants_titan_seg_filtered.txt: copy number data from variants_titan_seg (syn23554313) filtered for primary and first recurrents of patients in glass_gbm_idhwt_rt_tmz_local+dis_cna.txt, using:
```
while read line ; do cat glass_data/variants_titan_seg.txt  | grep "${line}-TP" | cut -f 1,2,3,4,5,7 >> glass_data/variants_titan_seg_filtered_temp.txt ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+dis.txt 
while read line ; do cat glass_data/variants_titan_seg.txt  | grep "${line}-R1" | cut -f 1,2,3,4,5,7 >> glass_data/variants_titan_seg_filtered_temp.txt ; done <patient_lists/glass_gbm_idhwt_rt_tmz_local+dis.txt 
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
Also contains created custom gene sets.

## patient_lists
WARNING: Many of the analyses rely on GLASS primary and recurrents being labelled as TP and R1. This is the case for patients inlcuded in the following lists but any alterations need checking.
gbm_idhwt_rt_tmz_local.txt: Discovery patients who's primary and first recurrent are GBM_IDHwt, local first recurrent, recieved rt+tmz.
glass_gbm_idhwt.txt: GLASS patients who's primary and first recurrent are GBM_IDHwt or GBM_IDHunknown, have RNA data available and not in the discovery cohort.
glass_gbm_idhwt_rt_tmz_local.txt: GLASS patients who's primary and first recurrent are GBM_IDHwt or GBM_IDHunknown, local first recurrent, recieved rt+tmz, have RNA data available and not in the discovery cohort.
glass_gbm_idhwt_rt_tmz_local+dis.txt: As with glass_gbm_idhwt_rt_tmz_local.txt but with discovery patients included.
glass_gbm_idhwt_rt_tmz_local+dis_cna.txt: As with glass_gbm_idhwt_rt_tmz_local+dis.txt but only those with GLASS copy number alteration data available.
glass_gbm_idhwt_rt_tmz_local_methylation+rna.txt: As with glass_gbm_idhwt_rt_tmz_local.txt but with discovery patients included and filtered for those with methylation data also available. 
mixed_gbm_idhwt_rt_tmz_local_methylation+rna.txt: As with glass_gbm_idhwt_rt_tmz_local_methylation+rna.txt but with discovery IDs for discovery patients.

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
#TFs
for p in ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_absolute/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute -s ; done 
for p in ranks/actual_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_actual/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual -s ; done 
for p in ranks/absolute_log2fc_tss/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_absolute_tss/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute_tss -s ; done 
for p in ranks/actual_log2fc_tss/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_actual_tss/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual_tss -s ; done 
for p in ranks/glass_absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_absolute_glass/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute_glass -s ; done 
for p in ranks/glass_actual_log2fc/*.rnk ; do pi=$(basename $p .rnk); [ ! -d gsea_outputs/outputs_actual_glass/${pi}*_${SIZE}* ] && qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual_glass -s ; done 
#c3.tft
for f in ranks/absolute_log2fc/* ; do p=$(basename $f .rnk) ;python scripts/convert_ensids_to_names.py -i $f -o ranks/absolute_log2fc_names/${p}.rnk ; done
for f in ranks/actual_log2fc/* ; do p=$(basename $f .rnk) ;python scripts/convert_ensids_to_names.py -i $f -o ranks/actual_log2fc_names/${p}.rnk ; done
for p in ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea_c3.tft.qsubsec PATIENT=${pi} MODE=absolute -s ; done 
for p in ranks/actual_log2fc/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea_c3.tft.qsubsec PATIENT=${pi} MODE=actual -s ; done 

```

## Process GSEA results

```
#Remove unnecesary files
mkdir temp 
for dir in gsea_outputs/*/* ; do ls ${dir}/* | grep -E -v 'PCGF2|CBX2|CBX7|CBX8|EZH2|CBX6|HSD17B8|SS18/SSX1|NR5A1|APOBEC3B|YBX1|SRSF9|PBXIP1|SIRT1|SUPT16H|TFAM|ARID3B|LHX2|JARID2|report' | xargs -I {} mv {} temp/; done 
for dir in gsea_outputs/*/* ; do rm -r ${dir}/edb ; done 
rm -r temp 

#Combine results and generate tables for normalised enrichment score, p-value and FDR
SIZES="1000 2000 5000"
SETS="outputs_actual outputs_actual_glass outputs_absolute outputs_absolute_glass outputs_actual_tss outputs_absolute_tss"

#Get JARID2 NES or ES for all patients
for SET in $SETS ; do for SIZE in 1000 ; do grep 'JARID2' gsea_outputs/${SET}/*_${SIZE}_*/gsea_report_for_na_*tsv | sed "s/_GTRD_${SIZE}_/\//" | tr '/' '\t' | cut -f 3,9,10 | awk '{ if ( $3=="---" ) {print $1"\t"$2} else {print $1"\t"$3} }' > reports/jarid2_results/${SET}_${SIZE}_JARID2_results.tsv ; done ; done
#manually create reports/jarid2_results/outputs_actual_glass_1000_JARID2_results+dis.tsv to include discovery patient results as GLASS ids.


#Get tables for all patients and genes
cat gsea_files/TFs_ENS_5000_GTRDv19_10_gencodev27.gmt | cut -f 1 > reports/all_tfs.txt
for SET in $SETS ; do for SIZE in $SIZES ; do mkdir reports/${SET}_${SIZE} ; for f in gsea_outputs/${SET}/*_${SIZE}_* ; do fi=$(basename ${f%_GTRD*}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '=='> reports/${SET}_${SIZE}/${fi}_table.txt ; done ; done ; done
for SET in ${SETS} ; do for SIZE in ${SIZES} ; do python scripts/combine_patients.py --genes reports/all_tfs.txt --directory reports/${SET}_${SIZE} -output reports/combined/${SET}_${SIZE}  ;done ; done 
Rscript scripts/gsea_violin_plots.py

#c3.tft
cat gene_sets/c3.tft.v7.4.symbols.gmt | cut -f 1 > reports/all_c3.tft.txt
for f in gsea_outputs/outputs_absolute/*c3.tft* ; do fi=$(basename ${f%_c3.tft*}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '=='> reports/outputs_absolute_c3.tft/${fi}_table.txt ; done 
python scripts/combine_patients.py --genes reports/all_c3.tft.txt --directory reports/outputs_absolute_c3.tft -output reports/combined/c3.tft 

```

## cell lines and tissue slices
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179649/suppl/GSE179649_Neuro_Organo_kallisto_gene_tpm.csv.gz -P tissue_slices/downloaded_data/

python script/process_cell_lines.py
python script/process_tissue_slices.py

Rscript scripts/plot_pca_tissue_slices_log2fc_project.R

TABLE=primary
TABLE=recurrent
TABLE=log2fc

for p in tissue_slices/ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk);qsubsec scripts/gsea_tissue_slices.qsubsec DIR=tissue_slices SIZE=1000 MODE=actual PATIENT=${pi} -s ; done
for p in tissue_slices/ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk);qsubsec scripts/gsea_tissue_slices.qsubsec DIR=tissue_slices SIZE=1000 MODE=absolute PATIENT=${pi} -s ; done
for p in cell_lines/ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk);qsubsec scripts/gsea_tissue_slices.qsubsec DIR=cell_lines SIZE=1000 MODE=actual PATIENT=${pi} -s ; done
for p in cell_lines/ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk);qsubsec scripts/gsea_tissue_slices.qsubsec DIR=cell_lines SIZE=1000 MODE=absolute PATIENT=${pi} -s ; done

for f in cell_lines/gsea_outputs/outputs_actual/* ; do fi=$(basename ${f%_GTRD*}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '=='> cell_lines/reports/outputs_actual/${fi}_table.txt ; done
for f in cell_lines/gsea_outputs/outputs_absolute/* ; do fi=$(basename ${f%_GTRD*}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '=='> cell_lines/reports/outputs_absolute/${fi}_table.txt ; done
for f in tissue_slices/gsea_outputs/outputs_actual/* ; do fi=$(basename ${f%_GTRD*}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '=='> tissue_slices/reports/outputs_actual/${fi}_table.txt ; done
for f in tissue_slices/gsea_outputs/outputs_absolute/* ; do fi=$(basename ${f%_GTRD*}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '=='> tissue_slices/reports/outputs_absolute/${fi}_table.txt ; done

python scripts/combine_patients.py --genes reports/all_tfs.txt --directory cell_lines/reports/outputs_actual --output cell_lines/reports/combined/cell_lines_actual
python scripts/combine_patients.py --genes reports/all_tfs.txt --directory cell_lines/reports/outputs_absolute --output cell_lines/reports/combined/cell_lines_absolute 
python scripts/combine_patients.py --genes reports/all_tfs.txt --directory tissue_slices/reports/outputs_actual --output tissue_slices/reports/combined/tissue_slices_actual 
python scripts/combine_patients.py --genes reports/all_tfs.txt --directory tissue_slices/reports/outputs_absolute --output tissue_slices/reports/combined/tissue_slices_absolute



```

# Wang-Diaz comparison
```
python scripts/convert_names_to_ensids.py -i downloaded_data/NIHMS1579377-supplement-Supplemental_Table_S2.txt -o downloaded_data/NIHMS1579377-supplement-Supplemental_Table_S2_ens.txt
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

#inverse - genes not in the leading edge across patients
while read line ; do f=${line}; cat ./gsea_outputs/outputs_actual/*${f}_*1000*/JARID2.tsv | grep 'No' ; done <patient_lists/gbm_idhwt_rt_tmz_local.txt | cut -f 2 | sort | uniq -c > analysis/leading_edge/counts_actual_1000_gbm_idhwt_rt_tmz_local_inverse.txt
TOTAL=$(wc -l patient_lists/gbm_idhwt_rt_tmz_local.txt|cut -d' ' -f1)
NUM=$(awk -v total=$TOTAL -v le=1 'BEGIN {printf "%.0f", total*le}') ;awk -v  num=$NUM -F" " '{if ($1>=num) print}' analysis/leading_edge/counts_actual_1000_gbm_idhwt_rt_tmz_local_inverse.txt | awk '{print $(NF)}' > analysis/leading_edge/le${LE}_actual_1000_gbm_idhwt_rt_tmz_local_inverse.txt 

#cell_lines
cat cell_lines/gsea_outputs/outputs_actual/*/JARID2.tsv | grep 'Yes' | cut -f 2 | sort | uniq -c > cell_lines/leading_edge/counts_cell_lines.txt
LE=50
LE=70
TOTAL=$(ls -l cell_lines/ranks/absolute_log2fc/ |wc -l)
NUM=$(awk -v total=$TOTAL -v le=0.$LE 'BEGIN {printf "%.0f", total*le}') ;awk -v  num=$NUM -F" " '{if ($1>=num) print}' cell_lines/leading_edge/counts_cell_lines.txt | awk '{print $(NF)}' > cell_lines/leading_edge/le${LE}_cell_lines.txt 

#tissue_slices
cat tissue_slices/gsea_outputs/outputs_actual/*/JARID2.tsv | grep 'Yes' | cut -f 2 | sort | uniq -c > tissue_slices/leading_edge/counts_tissue_slices.txt
LE=50
LE=70
TOTAL=$(ls -l tissue_slices/ranks/absolute_log2fc/ |wc -l)
NUM=$(awk -v total=$TOTAL -v le=0.$LE 'BEGIN {printf "%.0f", total*le}') ;awk -v  num=$NUM -F" " '{if ($1>=num) print}' tissue_slices/leading_edge/counts_tissue_slices.txt | awk '{print $(NF)}' > tissue_slices/leading_edge/le${LE}_tissue_slices.txt 
```

### get annotations for each tss
```
python scripts/get_tss_annotations.py
```
Creates tss_anotations.txt which contains a list of every tss with the following annotations: 
    - gene: the gene that the tss belongs to. Some tsss belong to multiple genes and in these cases the tss is listed multiple times, once for each gene:
    - transcripts: the transcripts that the tss belongs to.
    - gene filter: whether this gene was included in totalRNA, mRNA, or both (all) filtered gene lists.
    - gene status: whether the gene is a jbs, le50 or le70 gene
    - jbs_tss: whether the tss is bound by JARID2



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
Rscript scripts/heatmap.R --patients patient_lists/glass_gbm_idhwt_rt_tmz_local.txt --genes analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt --table tables/log2fc_glass.txt --nes_colour reports/jarid2_results/outputs_actual_glass_1000_JARID2_results.tsv --name glass_gbm_idhwt_rt_tmz_local_le70

python scripts/merge.py -i cell_lines/tables/log2fc_cell_lines.txt,tables/log2fc_all.txt -o cell_lines/tables/log2fc_cell_lines_merged_with_patients.txt
Rscript scripts/heatmap_cell_lines.R --patients patient_lists/gbm_idhwt_rt_tmz_local_+_cell_lines.txt --genes analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt --table cell_lines/tables/log2fc_cell_lines_merged_with_patients.txt --nes_colour tissue_slices/reports/outputs_actual_1000_JARID2_results_cell_tissue.tsv --rna_colour original_data/library_types.txt --name cell_lines_le70
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
#Some custom gene set names were modified after running GSEA. 
```
cp TFs_ENS_1000_GTRDv19_10_gencodev27.gmt gene_sets/TFs_ENS_1000_GTRDv19_10_gencodev27.gmt
SETS="TFs_ENS_1000_GTRDv19_10_gencodev27 custom_gbm_gene_sets" 
for s in $SETS; do qsubsec scripts/gsea_uvd.qsubsec RUN=run1 GMT=$s GENES=IDs; done
SETS="h.all.v7.4.symbols c2.cgp.v7.4.symbols c2.cp.v7.4.symbols c3.mir.mirdb.v7.4.symbols c3.tft.v7.4.symbols c5.go.bp.v7.4.symbols c5.go.mf.v7.4.symbols" 
for s in $SETS; do qsubsec scripts/gsea_uvd.qsubsec RUN=run1 GMT=$s GENES=symbols; done

qsubsec scripts/gsea_uvd.qsubsec RUN=run1 GMT=gbm_invasive_unconnected GENES=IDs


for f in gsea_uvd_outputs/run1/* ; do fi=$(basename ${f}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '==> ' | grep '[0123456789]'> gsea_uvd_outputs/reports/${fi}_table.txt ; done

for f in gsea_uvd_outputs/run1/gbm_inv*3776* ; do fi=$(basename ${f}) ; tail ${f}/gsea_report_for_na_pos*tsv ${f}/gsea_report_for_na_neg*tsv -n +2 | grep -v '==> ' | grep '[0123456789]'> gsea_uvd_outputs/reports/${fi}_table.txt ; done



#Process results to combine GSEA NES and FDR with chi squared of numbers of increasing or decreasing significant DEGs
for s in $SETS; do qsubsec scripts/process_uvd_deseq2_results.qsubsec SET=$s ; done

#Plot results in a network
Rscript scripts/process_uvd_deseq2_results_for_network.R --gmt custom_gbm_gene_sets
Rscript scripts/process_uvd_deseq2_results_for_network.R --gmt h.all.v7.4.symbols

Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultscustom_gbm_gene_sets_for_network.txt --sets gsea_uvd_outputs/reports/custom_gbm_gene_sets_run1_plot.txt --name custom_plot

Rscript scripts/plot_processed_results.py --processed deseq2_uvd/processed_resultsh.all.v7.4.symbols_for_network.txt --sets gsea_uvd_outputs/reports/h.all.v7.4.symbols_run1_plot.txt --name hallmark_plot

```
## Variants
```
# Glass copy number

Rscript scripts/plot_cn_freq.R

# Glass point variants

#filter for deleterious
grep -E 'DE_NOVO_START_IN_FRAME|DE_NOVO_START_OUT_FRAME|FRAME_SHIFT_DEL|FRAME_SHIFT_INS|IN_FRAME_DEL|IN_FRAME_INS|MISSENSE|NONSENSE|NONSTOP|START_CODON_DEL|START_CODON_INS|START_CODON_SNP' glass_data/variants_anno_20201109.csv > glass_data/variants_anno_20201109_deleterious.csv
grep 'SILENT' glass_data/variants_anno_20201109.csv > glass_data/variants_anno_20201109_silent.csv
python scripts/process_glass_variants.py 
```

## Motif sequences
### MEME
```
#get sequences
cat gsea_files/Min2Exp_OverlapTFpeaksAndPromoters_1000_GTRDv19_10_gencodev27.txt | cut -f 11,16,17,18,33 | grep 'JARID2' > sequences/Min2Exp_OverlapTFpeaksAndPromoters_1000_GTRDv19_10_gencodev27_JARID2.txt
mv ./sequences/control_peaks.txt temp.txt ; cat analysis/leading_edge/le100_actual_1000_gbm_idhwt_rt_tmz_local_inverse.txt | while read line ; do grep $line sequences/Min2Exp_OverlapTFpeaksAndPromoters_1000_GTRDv19_10_gencodev27_JARID2.txt >> ./sequences/control_peaks.txt ; done
mv ./sequences/le70_peaks.txt temp.txt ; cat analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt | while read line ; do grep $line sequences/Min2Exp_OverlapTFpeaksAndPromoters_1000_GTRDv19_10_gencodev27_JARID2.txt >> ./sequences/le70_peaks.txt ; done
awk '{$3 = $3-=500 ; $4 = $4+=500}1' OFS='\t' ./sequences/le70_peaks.txt | sort -u > ./sequences/le70_sequence_positions.txt
awk '{$3 = $3-=500 ; $4 = $4+=500}1' OFS='\t' ./sequences/control_peaks.txt | sort -u > ./sequences/control_sequence_positions.txt
python sequences/merge_peak_sequence_positions.py
mv ./sequences/le70_sequences.fa temp.txt ; cat ./sequences/le70_sequence_positions_merged.txt | while read line ; do g=$(echo $line | cut -f1 -d" " | cut -f1 -d".") ; c=$(echo $line | cut -f2 -d" ") ; s=$(echo $line | cut -f3 -d" ") ; e=$(echo $line | cut -f4 -d" ") ; r=${c}:${s}-${e} ; samtools faidx downloaded_data/hg38.fasta $r --mark-strand custom,_${g},_ >>  ./sequences/le70_sequences.fa ; done
mv ./sequences/control_sequences.fa temp.txt ; cat ./sequences/control_sequence_positions_merged.txt | while read line ; do g=$(echo $line | cut -f1 -d" " | cut -f1 -d".") ; c=$(echo $line | cut -f2 -d" ") ; s=$(echo $line | cut -f3 -d" ") ; e=$(echo $line | cut -f4 -d" ") ; r=${c}:${s}-${e} ; samtools faidx downloaded_data/hg38.fasta $r --mark-strand custom,_${g},_ >>  ./sequences/control_sequences.fa ; done
awk '{$2 = $2-=500 ; $3 = $3+=500}1' OFS='\t'  original_data/JARID2promotersPeaks.bed | sort -u > ./sequences/JARID2promotersPeaks_sequence_positions.txt
mv ./sequences/JARID2promotersPeaks_sequences.fa temp.txt ; cat ./sequences/JARID2promotersPeaks_sequence_positions.txt | while read line ; do c=$(echo $line | cut -f1 -d" ") ; s=$(echo $line | cut -f2 -d" ") ; e=$(echo $line | cut -f3 -d" ") ; r=${c}:${s}-${e}  ; samtools faidx downloaded_data/hg38.fasta $r  >>  ./sequences/JARID2promotersPeaks_sequences.fa ; done
#merge did not make a difference for JARID2promotersPeaks

### MEME-ChIP

awk '{N=$3+$4 ; $3=int(N/2-250) ; $4=int(N/2+250)}1' OFS='\t' ./sequences/le70_peaks.txt | sort -u > ./sequences/le70_sequence_positions_meme-chip.txt
awk '{N=$3+$4 ; $3=int(N/2-250) ; $4=int(N/2+250)}1' OFS='\t' ./sequences/control_peaks.txt | sort -u > ./sequences/control_sequence_positions_meme-chip.txt
mv ./sequences/le70_sequences_meme-chip.fa temp.txt ; cat ./sequences/le70_sequence_positions_meme-chip.txt | while read line ; do g=$(echo $line | cut -f1 -d" " | cut -f1 -d".") ; c=$(echo $line | cut -f2 -d" ") ; s=$(echo $line | cut -f3 -d" ") ; e=$(echo $line | cut -f4 -d" ") ; r=${c}:${s}-${e} ; samtools faidx downloaded_data/hg38.fasta $r --mark-strand custom,_${g},_ >>  ./sequences/le70_sequences_meme-chip.fa ; done
mv ./sequences/control_sequences_meme-chip.fa temp.txt ; cat ./sequences/control_sequence_positions_meme-chip.txt | while read line ; do g=$(echo $line | cut -f1 -d" " | cut -f1 -d".") ; c=$(echo $line | cut -f2 -d" ") ; s=$(echo $line | cut -f3 -d" ") ; e=$(echo $line | cut -f4 -d" ") ; r=${c}:${s}-${e} ; samtools faidx downloaded_data/hg38.fasta $r --mark-strand custom,_${g},_ >>  ./sequences/control_sequences_meme-chip.fa ; done
```


