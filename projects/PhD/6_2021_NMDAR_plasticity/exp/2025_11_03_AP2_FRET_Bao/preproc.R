# Sasha's data preprocessing, HPCA+AP2 FRET

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_11_03_AP2_FRET_Bao')


df.full <- read.csv('df_FRET.csv') %>%
  select(-'...1') %>%
  rowwise() %>%
  mutate(parts = list(str_split(id, "_"))[[1]],
         len = length(parts),
         ch = parts[[len]],
         app_time = parts[[len-1]],
         id = str_remove(id, "_[^_]+_[^_]+$")) %>%
  ungroup() %>%
  select(-len, -parts) %>%
  mutate(app_time = case_match(app_time, '2' ~ '0.5', '3' ~ '60')) %>%
  mutate_if(is.character, factor)


write.csv(df.full, 'df_processed.csv')
