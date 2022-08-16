library(tidyr)
library(ggupset)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

#args = c('/Volumes/metagenomics_drive/liftover_project/test_files/merged_species_TE_df.TEBAG.csv')
tmp <- read_csv(file.path(args[1]))

print('file read in')

tmp %>% 
  mutate(Human = TRUE) %>% 
  mutate(insertion_status = str_replace(insertion_status, pattern = '\\[', replacement = ''),
         insertion_status = str_replace(insertion_status, pattern = '\\]', replacement = ''),
         insertion_status = str_replace_all(insertion_status, pattern = ' ', replacement = ''),
         insertion_status = as.list(as.vector(insertion_status))) %>% 
  unnest(insertion_status) %>% 
  separate(insertion_status, sep = '[,]', into = c('Chimp', 'Gorilla', 'Orangutan', 'Bonobo')) %>% 
  select(insertion = element, Human, Chimp, Gorilla, Orangutan, Bonobo, `panc.1.2.3`:`ctrl.9.1.2`) %>%
  mutate(Chimp = as.logical(Chimp),
         Gorilla = as.logical(Gorilla),
         Orangutan = as.logical(Orangutan),
         Bonobo = as.logical(Bonobo)) %>% 
  pivot_longer(names_to = 'sample', cols = `panc.1.2.3`:`ctrl.9.1.2`) %>% 
  select(insertion:Bonobo) %>% 
  pivot_longer(names_to = 'species', cols = Human:Bonobo) %>% 
  filter(value == T) -> tmp_insert_df

print('tmp_insert_df created')

tmp %>% 
  mutate(Human = TRUE) %>% 
  mutate(insertion_status = str_replace(insertion_status, pattern = '\\[', replacement = ''),
         insertion_status = str_replace(insertion_status, pattern = '\\]', replacement = ''),
         insertion_status = str_replace_all(insertion_status, pattern = ' ', replacement = ''),
         insertion_status = as.list(as.vector(insertion_status))) %>% 
  unnest(insertion_status) %>% 
  separate(insertion_status, sep = '[,]', into = c('Chimp', 'Gorilla', 'Orangutan', 'Bonobo')) %>% 
  select(insertion = element, Human, Chimp, Gorilla, Orangutan, Bonobo, `panc.1.2.3`:`ctrl.9.1.2`) %>% 
  mutate(Chimp = as.logical(Chimp),
         Gorilla = as.logical(Gorilla),
         Orangutan = as.logical(Orangutan),
         Bonobo = as.logical(Bonobo)) %>% 
  pivot_longer(names_to = 'sample', cols = `panc.1.2.3`:`ctrl.9.1.2`) %>%
  filter(value >= 1) %>%
  select(insertion, sample, value) %>% 
  separate(sample, sep = '[.]', into = c('condition', NA, NA, NA)) -> tmp_value_df

print('tmp_value_df created')

merge(tmp_insert_df, 
      tmp_value_df,
      by = 'insertion') -> tmp_merge_df
  
print('merge finished')

tmp_merge_df %>% 
  #merge(reference_meta_data$gencode_tx_to_gene.df, by.x = 'insertion', by.y = 'enst') %>% 
  group_by(insertion, condition) %>%
  summarize(species = unique(list(species)), 
            # mean_exp = mean(log1p(value.y)/n()), sd_exp = sd(log1p(value.y))/n()) %>%
            # mean_exp = sum(value.y/n())) %>% 
            mean_exp = mean(log1p(value.y)), sd_exp = sd(log1p(value.y))) %>%
  ggplot(aes(x = species,y = mean_exp, color = condition)) + 
  #scale_color_manual(values = identity_color_pal) +
  # geom_errorbar(position = position_dodge(width = 0.45)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15)) + 
  ggupset::scale_x_upset(order_by = 'degree') + 
  ylab('mean ( log( insertion exp + 1) / number of insertions )') -> upset_plot

print('plot created')
ggsave(plot = upset_plot, file = 'upset_plot.pdf', height = 10, width = 10)
