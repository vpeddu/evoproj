library(rGREAT)
library(tidyverse)
library(readr)
library(ggupset)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Repitools)

#tmp <- read_csv(file.path( '/media/vikas/thiccy1/data/TEBAG/TE/TE.merged_species_TE_df.TEBAG.csv'))
# laptop
tmp <- read_csv(file.path('/Users/vikas/Documents/UCSC/lab/kim_lab/TEBAG/local_data/TE.merged_species_TE_df.TEBAG.csv'))

treatments<-c('panc','ctrl','covid')
data_cols <- which(grepl(paste(treatments, collapse = '|'), colnames(tmp)))
tmp<-tmp[which(tmp[,data_cols] > 3),] 
tmp<-tmp[complete.cases(tmp$insertion_x),]

#tmp<-tmp[-which(grepl('rich',tmp$insertion_x)),]
#tmp<-tmp[-which(grepl('tRNA',tmp$insertion_x)),]
#tmp<-tmp[-which(grepl('poly',tmp$insertion_x)),]
#tmp<-tmp[-which(grepl('rRNA',tmp$insertion_x)),]

create_species_list <-function(row){
  species_list <- list('Human')
  if (row[35]){ 
    species_list<-append(species_list,'Chimp')
  }  
  if (row[36]){ 
    species_list<-append(species_list,'Gorilla')
  }
  if (row[37]){ 
    species_list<-append(species_list,'Orangutan')
  }
  if (row[38]){ 
    species_list<-append(species_list,'Bonobo')
  }
  return(species_list)
}

tmp$species<- apply(tmp, 1, create_species_list)

tmp %>% 
  filter(Chimp == F & 
           Gorilla == F & 
           Orangutan == F & 
           Bonobo == F) %>%
  separate(insertion_x, c('TE','location'),'=') %>%
  separate(location, c('chr','position'),':') %>%
  separate(position, c('start','stop'),'-') %>%
  mutate(stop= gsub('_strand',"",stop))  %>%
  select(chr,start,stop,TE) %>%
  transform(start = as.numeric(start)) %>%
  transform(stop = as.numeric(stop)) -> human_bed_df

tmp %>% 
  filter(Chimp != F & 
           Gorilla != F & 
           Orangutan != F & 
           Bonobo != F) %>%
  separate(insertion_x, c('TE','location'),'=') %>%
  separate(location, c('chr','position'),':') %>%
  separate(position, c('start','stop'),'-') %>%
  mutate(stop= gsub('_strand',"",stop))  %>%
  select(chr,start,stop,TE) %>%
  transform(start = as.numeric(start)) %>%
  transform(stop = as.numeric(stop))-> non_human_bed_df
  

human_great_job = submitGreatJob(human_bed_df,
                     species = 'hg38', 
                     rule = 'twoClosest',
                     adv_twoDistance = 2)

human_great_job_tb = getEnrichmentTables(human_great_job, ontology = c("GO Biological Process"))
human_great_job_df<- as.data.frame(human_great_job_tb)
#tb = getEnrichmentTables(job)
human_great_job_res = plotRegionGeneAssociationGraphs(human_great_job)
human_enrichment <- annoGR2DF(human_great_job_res)



non_human_great_job = submitGreatJob(non_human_bed_df,
                                 species = 'hg38', 
                                 rule = 'twoClosest',
                                 adv_twoDistance = 2)

non_human_great_job_tb = getEnrichmentTables(non_human_great_job, ontology = c("GO Molecular Function"))
#tb = getEnrichmentTables(job)
non_human_great_job_res = plotRegionGeneAssociationGraphs(non_human_great_job)
non_human_enrichment <- annoGR2DF(non_human_great_job_res)


