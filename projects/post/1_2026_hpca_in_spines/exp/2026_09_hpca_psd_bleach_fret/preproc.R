# Sasha's data prepossessing, HPCA+PSD95 donor dequenching

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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2026_03_5_HPCA_PSD_FRET_bleach')

##### CH0 #####
df.ch0.data <- rbind(read.csv('data/25_11_27_cell1_3_ch0_25_11_27_cell1_3_ch3_up-labels.csv'),
                     read.csv('data/25_11_27_cell2_3_ch0_25_11_27_cell2_3_ch3_up-labels.csv')) %>%
                rename(dF0_int = "dF.F0_int" ) %>%
                mutate(id = str_remove(id, '_3_algn_ch0'), ch = 'ch0', data = 'target') %>%
                select(-X, -lab_id) %>%
                mutate_if(is.character, factor) %>%
                mutate(roi = as.factor(roi))

df.ch0.back <- rbind(read.csv('data/25_11_27_cell1_3_ch0_25_11_27_cell1_back-labels.csv'),
                     read.csv('data/25_11_27_cell2_3_ch0_25_11_27_cell2_back-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_ch0'), ch = 'ch0', data = 'back') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.ch0 <- rbind(df.ch0.data, df.ch0.back)
remove(df.ch0.data, df.ch0.back)

levels(df.ch0$id)

ggplot(data = df.ch0 %>% filter(base == 'simple'),
       aes(x = time, y = abs_int, color = id, fill = id, linetype = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal()


##### CH1 #####
df.ch1.data <- rbind(read.csv('data/25_11_27_cell1_3_ch1_25_11_27_cell1_3_ch3_up-labels.csv'),
                     read.csv('data/25_11_27_cell2_3_ch1_25_11_27_cell2_3_ch3_up-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_ch1'), ch = 'ch1', data = 'target') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.ch1.back <- rbind(read.csv('data/25_11_27_cell1_3_ch1_25_11_27_cell1_back-labels.csv'),
                     read.csv('data/25_11_27_cell2_3_ch1_25_11_27_cell2_back-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_ch1'), ch = 'ch1', data = 'back') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.ch1 <- rbind(df.ch1.data, df.ch1.back)
remove(df.ch1.data, df.ch1.back)

levels(df.ch1$id)

ggplot(data = df.ch1 %>% filter(base == 'simple'),
       aes(x = time, y = abs_int, color = id, fill = id, linetype = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal()


##### CH3 #####
df.ch3.data <- rbind(read.csv('data/25_11_27_cell1_3_ch3_25_11_27_cell1_3_ch3_up-labels.csv'),
                     read.csv('data/25_11_27_cell2_3_ch3_25_11_27_cell2_3_ch3_up-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_ch3'), ch = 'ch3', data = 'target') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.ch3.back <- rbind(read.csv('data/25_11_27_cell1_3_ch3_25_11_27_cell1_back-labels.csv'),
                     read.csv('data/25_11_27_cell2_3_ch3_25_11_27_cell2_back-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_ch3'), ch = 'ch3', data = 'back') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.ch3 <- rbind(df.ch3.data, df.ch3.back)
remove(df.ch3.data, df.ch3.back)

levels(df.ch3$id)

ggplot(data = df.ch3 %>% filter(base == 'simple'),
       aes(x = time, y = abs_int, color = id, fill = id, linetype = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal()


##### Fc #####
df.fc.data <- rbind(read.csv('data/25_11_27_cell1_3_Fc_25_11_27_cell1_3_ch3_up-labels.csv'),
                    read.csv('data/25_11_27_cell2_3_Fc_25_11_27_cell2_3_ch3_up-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_Fc'), ch = 'Fc', data = 'target') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.fc.back <- rbind(read.csv('data/25_11_27_cell1_3_Fc_25_11_27_cell1_back-labels.csv'),
                     read.csv('data/25_11_27_cell2_3_Fc_25_11_27_cell2_back-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_Fc'), ch = 'Fc', data = 'back') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.fc <- rbind(df.fc.data, df.fc.back)
remove(df.fc.data, df.fc.back)

levels(df.fc$id)

ggplot(data = df.fc %>% filter(base == 'simple'),
       aes(x = time, y = abs_int, color = id, fill = id, linetype = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal()


##### Ea #####
df.ea.data <- rbind(read.csv('data/25_11_27_cell1_3_E_A_25_11_27_cell1_3_ch3_up-labels.csv'),
                    read.csv('data/25_11_27_cell2_3_E_A_25_11_27_cell2_3_ch3_up-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_E_A'), ch = 'Ea', data = 'target') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.ea.back <- rbind(read.csv('data/25_11_27_cell1_3_E_A_25_11_27_cell1_back-labels.csv'),
                    read.csv('data/25_11_27_cell2_3_E_A_25_11_27_cell2_back-labels.csv')) %>%
  rename(dF0_int = "dF.F0_int" ) %>%
  mutate(id = str_remove(id, '_3_algn_E_A'), ch = 'Ea', data = 'back') %>%
  select(-X, -lab_id) %>%
  mutate_if(is.character, factor) %>%
  mutate(roi = as.factor(roi))

df.ea <- rbind(df.ea.data, df.ea.back)
remove(df.ea.data, df.ea.back)

levels(df.ea$id)

ggplot(data = df.ea %>% filter(base == 'simple'),
       aes(x = time, y = abs_int, color = id, fill = id, linetype = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal()


##### FIN BIND #####
df.output <- rbind(df.ch0,
                   df.ch1,
                   df.ch3,
                   df.fc,
                   df.ea)

levels(df.output$id)
levels(df.output$ch)
levels(df.output$data)

ggplot(data = df.output %>% filter(base == 'simple'),
       aes(x = time, y = abs_int, color = id, fill = id,
           linetype = data, shape = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  facet_wrap(~ch, scale = 'free')

write.csv(df.output, 'df_combined.csv')
