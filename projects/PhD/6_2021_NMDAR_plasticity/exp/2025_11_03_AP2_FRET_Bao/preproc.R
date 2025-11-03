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
         l = length(parts),
         ch = parts[[l]],
         app_time = parts[[l-1]])
  # mutate_if(is.character, factor)
  # mutate(roi = as.factor(roi),
  #        app_factor = as.factor(app))


bb <- df.full$parts[1]
