library(rGREAT)
library(tidyverse)
library(ggupset)
library(ggplot2)

tmp <- read_csv(file.path( '/media/vikas/thiccy1/data/TEBAG/TE/TE.merged_species_TE_df.TEBAG.csv'))

treatments<-c('panc','ctrl','covid')
data_cols <- which(grepl(paste(treatments, collapse = '|'), colnames(tmp)))
tmp<-tmp[which(tmp[,data_cols] > 5),] 
tmp<-tmp[complete.cases(tmp$insertion_x),]

#tmp<-tmp[-which(grepl('rich',tmp$insertion_x)),]
#tmp<-tmp[-which(grepl('tRNA',tmp$insertion_x)),]
#tmp<-tmp[-which(grepl('poly',tmp$insertion_x)),]
#tmp<-tmp[-which(grepl('rRNA',tmp$insertion_x)),]

create_species_list <-function(row){
  species_list <- c('Human')
  if (row[35]){ 
    species_list<-c(species_list,'Chimp')
  }  
  if (row[36]){ 
    species_list<-c(species_list,'Gorilla')
  }
  if (row[37]){ 
    species_list<-c(species_list,'Orangutan')
  }
  if (row[38]){ 
    species_list<-c(species_list,'Bonobo')
  }
  return(species_list)
}

tmp$species<- apply(tmp, 1, create_species_list)

tmp %>% 
  separate(insertion_x, c('TE','location'),'=') %>%
  separate(location, c('chr','position'),':') %>%
  separate(position, c('start','stop'),'-') %>%
  mutate(stop= gsub('_strand',"",stop)) %>%
  select(chr,start,stop,TE) -> bed_df

bed_df$start<-as.numeric(bed_df$start)
bed_df$stop<-as.numeric(bed_df$stop)

job = submitGreatJob(bed_df,
                     species = 'hg38')

