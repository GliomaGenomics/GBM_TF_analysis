
# Data sources

## downoaded_data
'''
#Download databases
cd downloaded_data
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.chr_patch_hapl_scaff.annotation.gtf.gz
wget http://gtrd.biouml.org/downloads/19.10/chip-seq/Homo_sapiens_meta_clusters.interval.gz
gunzip gencode.v27.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip Homo_sapiens_meta_clusters.interval.gz
cd ..
'''

## original_data 
Contains expression tables and metadata

# Expression inputs for GSEA - in 'ranks' folder

get_foldchange.py and get_foldchange_tss.py were run on a previous smaller cohort to get filtered lists of genes and TSSs, separately for those samples processed with total RNA or mRNA libraries: filtered_genelist_mrna.txt,filtered_genelist_total.txt,filtered_tsslist_mrna.txt,filtered_tsslist_total.txt.
get_foldchange_newrealease.py and get_foldchange_tss_newrelease.py were later run on the current cohort expression table using the same lists of genes and TSSs.


# Create files for running GSAE
```
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

# Run GSEA

```
for p in ranks/absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute -s ; done 
for p in ranks/actual_log2fc/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual -s ; done 
for p in ranks/absolute_log2fc_tss/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute_tss_stringent -s ; done 
for p in ranks/actual_log2fc_tss/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual_tss_stringent -s ; done 
for p in ranks/glass_absolute_log2fc/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=absolute_glass -s ; done 
for p in ranks/glass_actual_log2fc/*.rnk ; do pi=$(basename $p .rnk); qsubsec scripts/gsea.qsubsec PATIENT=${pi} SIZE=${SIZE} MODE=actual_glass -s ; done 
```

# Process GSEA results

```
#Remove unnecesary files
mkdir temp 
for dir in gsea_outputs/*/* ; do ls ${dir}/* | grep -E -v 'PCGF2|CBX2|CBX7|CBX8|EZH2|CBX6|HSD17B8|SS18/SSX1|NR5A1|APOBEC3B|YBX1|SRSF9|PBXIP1|SIRT1|SUPT16H|TFAM|ARID3B|LHX2|JARID2|report' | xargs -I {} mv {} temp/; done 
for dir in gsea_outputs/*/* ; do rm -r ${dir}/edb ; done 
rm -r temp 

#Combine results and generate tables for normalised enrichment score, p-value and FDR
SETS="1000_classic_absolute 1000_weighted_actual 2000_classic_absolute 2000_weighted_actual 5000_classic_absolute 5000_weighted_actual" 
#SETS="1000_classic_absolute_tss_stringent 1000_weighted_actual_tss_stringent"
#SETS="1000_classic_absolute_glass,1000_weighted_actual_glass"


```

# Get intermediate files
```
awk -F'\t' '{SUM=($2+$3)/2; print($7"\t"$1":"SUM)}' gsea_files/gencode.v27.chr_patch_hapl_scaff.annotation_PromotersTSS_1000.txt | sort -u > intermediate_files/transcript_to_tss_position.txt 
awk -F'\t' '{SUM=($2+$3)/2; print($5"\t"$1":"SUM)}' gsea_files/gencode.v27.chr_patch_hapl_scaff.annotation_PromotersTSS_1000.txt | sort -u > intermediate_files/gene_to_tss_position.txt 
```
