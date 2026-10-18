# Sasha's data preprocessing, HPCA+PSD95, ch0 only

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)

setwd('~/bio_note/projects/post/1_2026_hpca_in_spines/exp/2026_09_hpca_psd_20s')

##### PSD #####
df.psd.260903 <- rbind(read.csv('data/HPCA_2026/26_09_03_cell4_2_ch0_psd.csv')) %>%
  rename(dF0_int = "dF.F0_int") %>%
  mutate(id = str_remove(id, '_2_ch0_psd'),
         lab_id = 'psd', ch = 'ch0',
         rel_time = time - 10) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

ggplot(df.psd.260903 %>% filter(base == 'simple'),
       aes(x = rel_time, y = dF0_int, color = roi)) +
  geom_line()

##### OREOL #####
df.oreol.260903 <- rbind(read.csv('data/HPCA_2026/26_09_03_cell4_2_ch0_oreol.csv')) %>%
  rename(dF0_int = "dF.F0_int") %>%
  mutate(id = str_remove(id, '_2_ch0_oreol'),
         lab_id = 'oreol', ch = 'ch0',
         rel_time = time - 10) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

ggplot(df.oreol.260903 %>% filter(base == 'simple'),
       aes(x = rel_time, y = dF0_int, color = roi)) +
  geom_line()

##### SHAFT #####
df.shaft.260903 <- rbind(read.csv('data/HPCA_2026/26_09_03_cell4_2_ch0_shaft.csv')) %>%
  rename(dF0_int = "dF.F0_int") %>%
  mutate(id = str_remove(id, '_2_ch0_shaft'),
         lab_id = 'psd', ch = 'ch0',
         rel_time = time - 10) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

ggplot(df.shaft.260903 %>% filter(base == 'simple'),
       aes(x = rel_time, y = dF0_int, color = roi)) +
  geom_line()

##### SAVING #####


levels(df.hpca$id)
colnames(df.fret)
colnames(df.fret.psd)
colnames(df.hpca)

df.comb <- rbind(df.fret, df.fret.psd, df.hpca)

write.csv(df.comb, 'data/df_combined.csv')
