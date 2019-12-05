#!/bin/bash
# @ Author Charlotte Capitanchik
# Provide peak files and get the merged clusters of all peak files
# Then get crosslink coverage over these (here I am doing exactly what iCount clusters does, but it wasn't working for some reason so implemented with bedtools)
# Then format for DeSeq

# INPUT
files="PHO92_Sc_4TUiCLIP_IME4del_SPO_REP1_G.peaks.bed,PHO92_Sc_4TUiCLIP_IME4del_SPO_REP2_G.peaks.bed,PHO92_Sc_4TUiCLIP_IME4del_SPO_REP3_G.peaks.bed,PHO92_Sc_4TUiCLIP_IME4del_SPO_REP4_G.peaks.bed,PHO92_Sc_4TUiCLIP_WT_SPO_REP1_G.peaks.bed,PHO92_Sc_4TUiCLIP_WT_SPO_REP2_G.peaks.bed,PHO92_Sc_4TUiCLIP_WT_SPO_REP3_G.peaks.bed,PHO92_Sc_4TUiCLIP_WT_SPO_REP4_G.peaks.bed,WT_Sc_4TUiCLIP_WT-xl365_SPO_REP1_G.peaks.bed,PHO92_Sc_4TUiCLIP_IME4del-6Jxl_SPO_REP1_G.peaks.bed"
stem="PHO92_iCLIP"
cluster_dist="20"
folder="."


# MAIN
# Create a new directory to store all these results and a file containing what files were used to generate them
#mkdir ${folder}
#chmod -R 755 ${folder}
# Convert to a list we can iterate over
file_list=$(echo ${files} | tr "," "\n")
# Make a read me that contains the files we are using
echo ${file_list} > ${folder}/README.txt
# Convert to a line separated by spaces instead of commas
file_list_space=$(echo ${files} | sed s'/,/ /g')
# Make a merged megafile of peaks
echo ${file_list_space}
cat ${file_list_space} > ${folder}/${stem}_merged_peaks.bed
# Prepare crosslink sites
file_xlsites=$(echo ${file_list_space} | sed s'/\.peaks.bed/\.Aligned\.out\.sorted\.cdnacounts/g' | sed s'/m6A_Sc/\.\.\/\.\.\/results\/genome_xlsites\/PHO92_Sc/g')
# Make a merged megafile of crosslinks
cat ${file_xlsites} > ${folder}/${stem}_merged_xlsites.bed

# Group the files and sort them
bedtools sort -i ${folder}/${stem}_merged_xlsites.bed | bedtools groupby -g 1,2,3,6 -c 5 -o count | awk '{print $1 "\t" $2 "\t" $3 "\t.\t" $5 "\t" $4}' > ${folder}/${stem}_merged_xlsites.sorted.bed
bedtools sort -i ${folder}/${stem}_merged_peaks.bed | bedtools groupby -g 1,2,3,6 -c 5 -o count | awk '{print $1 "\t" $2 "\t" $3 "\t.\t" $5 "\t" $4}' > ${folder}/${stem}_merged_peaks.sorted.bed

# Call clusters on the merged peaks
bedtools merge -s -d ${cluster_dist} -i ${folder}/${stem}_merged_peaks.sorted.bed -c 4,5,6 -o distinct,sum,distinct > ${folder}/${stem}_${cluster_dist}_merged_clusters.bed
# For each site, find closest cluster to which assign the site to, then filter for sites that dont overlap but are within 3nt distance.
bedtools closest -s -d -t first -a ${folder}/${stem}_merged_xlsites.sorted.bed -b ${folder}/${stem}_${cluster_dist}_merged_clusters.bed | awk '$13 != 0 && $13 < 4 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > ${folder}/${stem}_${cluster_dist}_nearestTo.bed

# Merge crosslinks that are within 3nt of peak clusters
cat ${folder}/${stem}_${cluster_dist}_nearestTo.bed ${folder}/${stem}_${cluster_dist}_merged_clusters.bed | bedtools sort | bedtools merge -s -d 3 -c 4,5,6 -o distinct,sum,distinct > ${folder}/${stem}_${cluster_dist}_merged_3_clusters.bed

# For each crosslink file get the coverage over the joint clusters
for file in $file_xlsites; do
	echo $file
	filename=$(echo ${file} | sed s'/\.\.\/\.\.\/results\/genome_xlsites\///g' | sed s'/\.Aligned\.out\.sorted\.cdnacounts//g')
	echo $filename
	bedtools sort -i ${file} > ${filename}.sorted
	bedtools map -a ${folder}/${stem}_${cluster_dist}_merged_3_clusters.bed -b ${filename}.sorted -c 5 -o sum -s > ${folder}/${filename}.mergedClusterCoverage.${cluster_dist}.bed
done

# Merge all the coverage files into a file for use with DeSeq
Rscript --vanilla formatForDeSeq.R ${folder} mergedClusterCoverage.${cluster_dist}.bed
