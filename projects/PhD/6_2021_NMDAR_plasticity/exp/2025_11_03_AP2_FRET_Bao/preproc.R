# Sasha's data prepossessing, HPCA+AP2B1

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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_11_03_AP2_FRET_Bao')

##### CH0 SHAFT #####
df.ch0.shaft <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch0_shaft.csv'),
                      read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch0_shaft.csv'),
                      read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch0_shaft.csv'),
                      read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch0_shaft.csv'),
                      read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch0_shaft.csv'),
                      read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch0_shaft.csv'),
                      read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch0_shaft.csv')) %>%
                rename(dF0_int = "dF.F0_int" ) %>%
                mutate(id = str_remove(id, '_3_ch0'),
                   lab_id = 'shaft', ch = 'ch0',
                   rel_time = time - 40) %>%
                select(-X) %>%
                mutate_if(is.character, factor) %>%
                mutate(roi = as.factor(roi))

levels(df.ch0.shaft$id)

ggplot(data = df.ch0.shaft %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH1 SHAFT #####
df.ch1.shaft <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch1_shaft.csv'),
                      read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch1_shaft.csv'),
                      read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch1_shaft.csv'),
                      read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch1_shaft.csv'),
                      read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch1_shaft.csv'),
                      read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch1_shaft.csv'),
                      read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch1_shaft.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch1'),
         lab_id = 'shaft', ch = 'ch1',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch1.shaft$id)

ggplot(data = df.ch1.shaft %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH3 SHAFT #####
df.ch3.shaft <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch3_shaft.csv'),
                      read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch3_shaft.csv'),
                      read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch3_shaft.csv'),
                      read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch3_shaft.csv'),
                      read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch3_shaft.csv'),
                      read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch3_shaft.csv'),
                      read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch3_shaft.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch3'),
         lab_id = 'shaft', ch = 'ch3',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch3.shaft$id)

ggplot(data = df.ch3.shaft %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH0 HALO #####
df.ch0.halo <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch0_halo.csv'),
                     read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch0_halo.csv'),
                     read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch0_halo.csv'),
                     read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch0_halo.csv'),
                     read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch0_halo.csv'),
                     read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch0_halo.csv'),
                     read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch0_halo.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch0'),
         lab_id = 'halo', ch = 'ch0',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch0.halo$id)

ggplot(data = df.ch0.halo %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH1 HALO #####
df.ch1.halo <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch1_halo.csv'),
                     read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch1_halo.csv'),
                     read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch1_halo.csv'),
                     read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch1_halo.csv'),
                     read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch1_halo.csv'),
                     read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch1_halo.csv'),
                     read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch1_halo.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch1'),
         lab_id = 'halo', ch = 'ch1',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch1.halo$id)

ggplot(data = df.ch1.halo %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH3 HALO #####
df.ch3.halo <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch3_halo.csv'),
                     read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch3_halo.csv'),
                     read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch3_halo.csv'),
                     read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch3_halo.csv'),
                     read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch3_halo.csv'),
                     read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch3_halo.csv'),
                     read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch3_halo.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch3'),
         lab_id = 'halo', ch = 'ch3',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch3.halo$id)

ggplot(data = df.ch3.halo %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH0 PSD #####
df.ch0.psd <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch0_psd.csv'),
                    read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch0_psd.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch0_psd.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch0_psd.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch0_psd.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch0_psd.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch0_psd.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch0'),
         lab_id = 'psd', ch = 'ch0',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch0.psd$id)

ggplot(data = df.ch0.psd %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH1 PSD #####
df.ch1.psd <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch1_psd.csv'),
                    read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch1_psd.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch1_psd.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch1_psd.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch1_psd.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch1_psd.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch1_psd.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch1'),
         lab_id = 'psd', ch = 'ch1',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch1.psd$id)

ggplot(data = df.ch1.psd %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH3 PSD #####
df.ch3.psd <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_ch3_psd.csv'),
                    read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch3_psd.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch3_psd.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch3_psd.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch3_psd.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch3_psd.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch3_psd.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch3'),
         lab_id = 'psd', ch = 'ch3',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch3.psd$id)

ggplot(data = df.ch3.psd %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH0 CCP #####
df.ch0.ccp <- rbind(read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch0_ccp.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch0_ccp.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch0_ccp.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch0_ccp.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch0_ccp.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch0_ccp.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch0'),
         lab_id = 'ccp', ch = 'ch0',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch0.ccp$id)

ggplot(data = df.ch0.ccp %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH1 CCP #####
df.ch1.ccp <- rbind(read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch1_ccp.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch1_ccp.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch1_ccp.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch1_ccp.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch1_ccp.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch1_ccp.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch1'),
         lab_id = 'ccp', ch = 'ch1',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch1.ccp$id)

ggplot(data = df.ch1.ccp %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### CH3 CCP #####
df.ch3.ccp <- rbind(read.csv('df/df_25_09_10_02/25_09_10_cell2_3_ch3_ccp.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_ch3_ccp.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_ch3_ccp.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_ch3_ccp.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_ch3_ccp.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_ch3_ccp.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_ch3'),
         lab_id = 'ccp', ch = 'ch3',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.ch3.ccp$id)

ggplot(data = df.ch3.ccp %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### FC SHAFT #####
df.fc.shaft <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_Fc_shaft.csv'),
                     read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Fc_shaft.csv'),
                     read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Fc_shaft.csv'),
                     read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Fc_shaft.csv'),
                     read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Fc_shaft.csv'),
                     read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Fc_shaft.csv'),
                     read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Fc_shaft.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Fc'),
         lab_id = 'shaft', ch = 'Fc',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.fc.shaft$id)

ggplot(data = df.fc.shaft %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### FC HALO #####
df.fc.halo <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_Fc_halo.csv'),
                    read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Fc_halo.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Fc_halo.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Fc_halo.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Fc_halo.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Fc_halo.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Fc_halo.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Fc'),
         lab_id = 'halo', ch = 'Fc',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.fc.halo$id)

ggplot(data = df.fc.halo %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### FC PSD #####
df.fc.psd <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_Fc_psd.csv'),
                    read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Fc_psd.csv'),
                    read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Fc_psd.csv'),
                    read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Fc_psd.csv'),
                    read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Fc_psd.csv'),
                    read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Fc_psd.csv'),
                    read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Fc_psd.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Fc'),
         lab_id = 'psd', ch = 'Fc',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.fc.psd$id)

ggplot(data = df.fc.psd %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### FC CCP #####
df.fc.ccp <- rbind(read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Fc_ccp.csv'),
                   read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Fc_ccp.csv'),
                   read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Fc_ccp.csv'),
                   read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Fc_ccp.csv'),
                   read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Fc_ccp.csv'),
                   read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Fc_ccp.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Fc'),
         lab_id = 'ccp', ch = 'Fc',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.fc.ccp$id)

ggplot(data = df.fc.ccp %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### EAPP SHAFT #####
df.eapp.shaft <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_Eapp_shaft.csv'),
                       read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Eapp_shaft.csv'),
                       read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Eapp_shaft.csv'),
                       read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Eapp_shaft.csv'),
                       read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Eapp_shaft.csv'),
                       read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Eapp_shaft.csv'),
                       read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Eapp_shaft.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Eapp'),
         lab_id = 'shaft', ch = 'Eapp',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.eapp.shaft$id)

ggplot(data = df.eapp.shaft %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### EAPP HALO #####
df.eapp.halo <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_Eapp_halo.csv'),
                       read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Eapp_halo.csv'),
                       read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Eapp_halo.csv'),
                       read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Eapp_halo.csv'),
                       read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Eapp_halo.csv'),
                       read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Eapp_halo.csv'),
                       read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Eapp_halo.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Eapp'),
         lab_id = 'halo', ch = 'Eapp',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.eapp.halo$id)

ggplot(data = df.eapp.halo %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### EAPP PSD #####
df.eapp.psd <- rbind(read.csv('df/df_25_09_10_01/25_09_10_cell1_3_Eapp_psd.csv'),
                      read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Eapp_psd.csv'),
                      read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Eapp_psd.csv'),
                      read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Eapp_psd.csv'),
                      read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Eapp_psd.csv'),
                      read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Eapp_psd.csv'),
                      read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Eapp_psd.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Eapp'),
         lab_id = 'psd', ch = 'Eapp',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.eapp.psd$id)

ggplot(data = df.eapp.psd %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### EAPP CCP #####
df.eapp.ccp <- rbind(read.csv('df/df_25_09_10_02/25_09_10_cell2_3_Eapp_ccp.csv'),
                     read.csv('df/df_25_09_10_04/25_09_10_cell4_3_Eapp_ccp.csv'),
                     read.csv('df/df_25_09_10_05/25_09_10_cell5_3_Eapp_ccp.csv'),
                     read.csv('df/df_25_09_11_01/25_09_11_cell1_3_Eapp_ccp.csv'),
                     read.csv('df/df_25_10_09_01/25_10_09_cell1_3_Eapp_ccp.csv'),
                     read.csv('df/df_25_10_09_02/25_10_09_cell2_3_Eapp_ccp.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_Eapp'),
         lab_id = 'ccp', ch = 'Eapp',
         rel_time = time - 40) %>%
  select(-X) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

levels(df.eapp.ccp$id)

ggplot(data = df.eapp.ccp %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  theme(legend.position = 'None') +
  facet_wrap(~id, ncol = 3)

##### FIN BIND #####
df.output <- rbind(df.ch0.ccp,
                   df.ch0.halo,
                   df.ch0.psd,
                   df.ch0.shaft,
                   df.ch1.ccp,
                   df.ch1.halo,
                   df.ch1.psd,
                   df.ch1.shaft,
                   df.ch3.ccp,
                   df.ch3.halo,
                   df.ch3.psd,
                   df.ch3.shaft,
                   df.fc.ccp,
                   df.fc.halo,
                   df.fc.psd,
                   df.fc.shaft,
                   df.eapp.ccp,
                   df.eapp.halo,
                   df.eapp.psd,
                   df.eapp.shaft)

levels(df.output$id)
write.csv(df.output, 'df_combined.csv')
