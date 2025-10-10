# Sasha's data preprocessing, HPCA+PSD95

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(ggsci)
library(introdataviz)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_06_19_Sasha_HPCA-PSD')

df.fret <- rbind(read.csv('data/FRET/07_05_25_cell4_Fc_oreol.csv'),
                read.csv('data/FRET/07_05_25_cell7_Fc_oreol.csv'),
                read.csv('data/FRET/07_05_25_cell9_Fc_oreol.csv'),
                read.csv('data/FRET/07_05_25_cell11_Fc_oreol.csv'),
                read.csv('data/FRET/07_05_25_cell16_Fc_oreol.csv')) %>%
            mutate(id = str_remove(id, '_Fc'),
                   lab_id = 'oreol', ch = 'fc',
                   app_factor = '20', app = 20,
                   rel_time = index - 30) %>%
            select(-X) %>%
            mutate_if(is.character, factor) %>%
            mutate(roi = as.factor(roi)) %>%
            rename(df = "dF.F0_int" )

levels(df.fret$id)


df.fret.psd <- rbind(read.csv('data/FRET/07_05_25_cell4_Fc_psd.csv'),
                     read.csv('data/FRET/07_05_25_cell7_Fc_psd.csv'),
                     read.csv('data/FRET/07_05_25_cell9_Fc_psd.csv'),
                     read.csv('data/FRET/07_05_25_cell11_Fc_psd.csv'),
                     read.csv('data/FRET/07_05_25_cell16_Fc_psd.csv')) %>%
              mutate(id = str_remove(id, '_Fc'),
                     lab_id = 'psd', ch = 'fc',
                     app_factor = '20', app = 20,
                     rel_time = index - 30) %>%
              select(-X) %>%
              mutate_if(is.character, factor) %>%
              mutate(roi = as.factor(roi)) %>%
              rename(df = "dF.F0_int" )

levels(df.fret.psd$id)


df.hpca <- read.csv('data/HPCA/df.csv') %>%
           select(-X) %>%
           mutate(id = str_remove(id, '_ch0'),
                  ch='hpca',
                  roi = as.factor(roi),
                  app_factor = as.factor(app)) %>%
           mutate_if(is.character, factor) %>%
           mutate(rel_time = index - 30,
                  id = recode(id, 'cell15' = '07_05_25_cell15', 'cell7' = '06_06_25_cell7'))

levels(df.hpca$id)


colnames(df.fret)
colnames(df.fret.psd)
colnames(df.hpca)

df.comb <- rbind(df.fret, df.fret.psd, df.hpca)

write.csv(df.comb, 'data/df_combined.csv')
