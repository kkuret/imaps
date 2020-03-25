library(ggplot2)
library(data.table)
library(dplyr)
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)

dir = args[1] #"m6A_WT_miCLIP" #
searchpat = args[2] #"mergedClusterCoverage.50.bed"#

mfiles = list.files(path=dir, pattern=searchpat, full.names=TRUE)
mmaps = lapply(mfiles, fread)

mmaps_names <- gsub(paste0(".", searchpat),"", mfiles)
mmaps_names <- gsub(".*/","", mmaps_names)

names(mmaps) <- mmaps_names

for (i in 1:length(mmaps)){
  mmaps[[i]]$sample <- mmaps_names[[i]]
  mmaps[[i]] <- mmaps[[i]] %>% select(V1,V2,V3,V6,V7,sample) %>% transmute(gene=paste0(V1,":",V2,"-",V3,";",V6), count=V7, sample=sample)
  rownames(mmaps[[i]]) <- NULL
}

summary_df <- do.call(rbind, mmaps)
rownames(summary_df) <- NULL
summary_df <- reshape(summary_df, idvar = "gene", timevar = "sample", direction = "wide")
names(summary_df) <- gsub("count.", "", names(summary_df))

fwrite(summary_df, file=paste0(dir,"/DeSeq_table.tsv"), sep="\t", col.name=TRUE, quote=FALSE)
