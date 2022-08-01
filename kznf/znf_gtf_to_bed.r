library(rtracklayer)
library(svMisc)

args <- commandArgs(trailingOnly = TRUE)

gtf<-import('/Users/vikas/Downloads/gencode.v35.annotation.gtf')

gtf_df<-as.data.frame(gtf)

#Extracted from supplement 2 of Trono et al 2017
krab_genes<-read.csv("/Users/vikas/Documents/UCSC/rotations/Haussler/HS_ZNF_extraction.csv", header = FALSE)

#removing rows that only have a dot for ensg#
krab_genes<-krab_genes[-which(grepl('[.]', krab_genes$V3)),]

krab_rows<-c()
for(i in 1:nrow(krab_genes)){ 
  #progress(i, nrow(krab_genes))
  temp_index<-which(grepl(krab_genes$V3[i], gtf_df$gene_id))  
  krab_rows<-append(temp_index,krab_rows)
  print(length(krab_rows))
}

gtf_krab_filtered<-gtf_df[krab_rows,]

gtf_krab_filtered_granges <- makeGRangesFromDataFrame(gtf_krab_filtered)

export(gtf_krab_filtered_granges, format = "bed", con = "krab_znf_filtered.bed")
export()