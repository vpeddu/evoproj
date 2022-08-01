library(ggupset)
library(reshape2)
library(ggplot2)
library(pbapply)      
library(cowplot)
library(tidyverse)
# run these to replace python True/False with R TRUE/FALSE
#sed 's/True/TRUE/g' copy.merged_species_TE_df.TEBAG.csv > inbetween.txt
#sed 's/False/FALSE/g' inbetween.txt > oof.csv


x<-read.csv('/Users/vikas/Documents/UCSC/lab/kim_lab/TEBAG/oof.csv')
x<-x[-which(grepl('rich',x$insertion_x)),]
x<-x[-which(grepl('tRNA',x$insertion_x)),]
x<-x[-which(grepl('poly',x$insertion_x)),]
x<-x[-which(grepl('rRNA',x$insertion_x)),]

treatments<-c('panc','ctrl','covid')
data_cols <- which(grepl(paste(treatments, collapse = '|'), colnames(x)))
x<-x[which(x[,data_cols] > 5),] 
x<-x[complete.cases(x[,data_cols]),]
x$species <-NA

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


x$species<- pbapply(x, 1, create_species_list)


x$insertion_x <- NULL
x$Human<- NULL
x$Orangutan<- NULL
x$Chimp<- NULL
x$Bonobo<-NULL
x$Gorilla<- NULL
x$X<-NULL

y<-melt(x, id.vars = c('insertion_y', 'species'))

z<-head(y)
z$cohort = strsplit(z$variable, '[.]')[[1]][1]

ssplit<-function(r){ 
  return(strsplit(r$variable, '[.]')[[1]][1])
}


#y<-y[y$value > 0, ]

y$cohort = apply(y, 1, ssplit)
y$value<-as.numeric(y$value)



z<-head(y, 10000)
z$cohort = apply(z, 1, ssplit)

z$value<-as.numeric(z$value)

# ggplot(z, aes(x = species, y = value, color = cohort)) + 
#   geom_jitter(height = 0, width = 0.1) + 
#   theme_classic() +
#   scale_y_log10() +
#   scale_x_upset()
# 
# ggplot(z, aes(x = species, y = value, fill = cohort , group = factor(cohort))) + 
#   geom_boxplot() + 
#   #geom_jitter(height = 0, width = 0.1) + 
#   theme_classic() +
#   scale_x_upset() +
#   scale_y_log10() 

y<-y[which(grepl('Alu',y$insertion_y)),]
panc<-y[y$cohort == 'panc',]
covid <- y[y$cohort == 'covid',]
ctrl <-y[y$cohort == 'ctrl',]

view_intersections = list(c('Human','Chimp',"Orangutan",'Gorilla','Bonobo'),
                            c('Human','Chimp',"Orangutan",'Gorilla'),
                            c('Human','Chimp','Gorilla'),
                            c('Human','Chimp'),
                            c('Human'))

panc_plot<-ggplot(panc, aes(x = species, y = value)) + 
  #geom_boxplot() + 
  geom_point() + 
  theme_classic() +
  scale_x_upset(intersections = view_intersections) +
  scale_y_log10(limits = c(1,1e5)) +
  ggtitle('panc')

covid_plot<-ggplot(covid, aes(x = species, y = value)) + 
  #geom_boxplot() + 
  geom_point() + 
  theme_classic() +
  scale_x_upset(intersections = view_intersections) +
  scale_y_log10(limits = c(1,1e5)) +
  ggtitle('covid')

ctrl_plot<-ggplot(ctrl, aes(x = species, value)) + 
  #geom_boxplot() + 
  geom_point() + 
  theme_classic() +
  scale_x_upset(intersections = view_intersections) +
  scale_y_log10(limits = c(1,1e5)) +
  ggtitle('ctrl')

panel<-plot_grid(panc_plot,covid_plot,ctrl_plot)
panel

ggplot(y, aes(x = species , y = value)) + 
  geom_bar(stat = 'identity') + 
  scale_x_upset()

######
######
######
######
file_list = list.files('/Users/vikas/Documents/UCSC/lab/kim_lab/TEBAG',pattern = '*.TEBAG.csv', full.names = TRUE)

for(i in file_list){
  base = strsplit(basename(i), '[.]')[[1]][1]
  print(base)
  tmp_file = read.csv(i)
  tmp_file$sample <- base
  if(exists("sample_df")){ 
    sample_df<-rbind(tmp_file, sample_df)
  }
  else{ 
    sample_df = tmp_file
    }
}


sample_df[sample_df=='True']<-TRUE
sample_df[sample_df=='False']<-FALSE

sample_df$species <-NA
for(i in 1:nrow(sample_df)){
  species_list = c(sample_df['Human')
  if (sample_df$Chimp[i]){ 
    species_list<-c(species_list,'Chimp')
  }  
  if (sample_df$Gorilla[i]){ 
    species_list<-c(species_list,'Gorilla')
  }
  if (sample_df$Orangutan[i]){ 
    species_list<-c(species_list,'Orangutan')
  }
  if (sample_df$Bonobo[i]){ 
    species_list<-c(species_list,'Bonobo')
  }
  # if(length(species_list) == 0){
  #   species_list = c('Only_Human')
  # }
  print(species_list)
  sample_df$species[i] = list(species_list)
}


ggplot(sample_df, aes(x = species, y = TPM, color = sample)) + 
  geom_jitter(height = 0, width = 0.1) + 
  theme_classic() +
  scale_y_log10() +
  scale_x_upset()
