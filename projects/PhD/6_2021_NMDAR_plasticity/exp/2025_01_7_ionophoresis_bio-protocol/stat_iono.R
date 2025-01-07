# NMDA ionophoresis, ionophoresis cloud analysis for bio-protocol
# Copyright © 2025 Borys Olifirov

require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol')

##### DATA PREPROCESSING #####
df.abs.raw <- rbind(read.csv('cloud/24_05_30_cell1_ch1_zone_mask_abs.csv'),
                    read.csv('cloud/24_05_30_cell2_ch1_zone_mask_abs.csv'),
                    read.csv('cloud/24_05_30_cell3_ch1_zone_mask_abs.csv'),
                    read.csv('cloud/24_05_30_cell4_ch1_zone_mask_abs.csv'),
                    read.csv('cloud/24_05_30_cell5_ch1_zone_mask_abs.csv'))

baseline.frames <- 5
start.idx.1st <- 5
start.idx.2nd <- 57
start.idx.3d <- 109
start.idx.4th <- 162

df.dF.raw <- rbind(read.csv('cloud/24_05_30_cell1_ch1_zone_mask_ΔF.csv'),
                   read.csv('cloud/24_05_30_cell2_ch1_zone_mask_ΔF.csv'),
                   read.csv('cloud/24_05_30_cell3_ch1_zone_mask_ΔF.csv'),
                   read.csv('cloud/24_05_30_cell4_ch1_zone_mask_ΔF.csv'),
                   read.csv('cloud/24_05_30_cell5_ch1_zone_mask_ΔF.csv')) %>%
             select(-X)
             # filter(index <= 213) %>%
             # mutate(roi = as.factor(roi),
             #        i_app = as.factor(case_when(index %in% seq(start.idx.1st - baseline.frames,start.idx.2nd) ~ '25',
             #                                    index %in% seq(start.idx.2nd - baseline.frames,start.idx.3d) ~ '50',
             #                                    index %in% seq(start.idx.3d - baseline.frames,start.idx.4th) ~ '75',
             #                                    index %in% seq(start.idx.4th - baseline.frames,213) ~ '100',
             #                                    .default = 'out'))) %>%
             # group_by(id,roi,i_app) %>%
             # mutate(index_app = seq(0,57))

df.dF.25 <- df.dF.raw %>%
  filter(index %in% seq(start.idx.1st - baseline.frames,start.idx.2nd)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '25') %>%
  ungroup()
df.dF.50 <- df.dF.raw %>%
  filter(index %in% seq(start.idx.2nd - baseline.frames,start.idx.3d)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '50') %>%
  ungroup()
df.dF.75 <- df.dF.raw %>%
  filter(index %in% seq(start.idx.3d - baseline.frames,start.idx.4th-1)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '75') %>%
  ungroup()
df.dF.100 <- df.dF.raw %>%
  filter(index %in% seq(start.idx.4th - baseline.frames,214)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '100') %>%
  ungroup()

df.dF.sweep <- rbind(df.dF.25, df.dF.50, df.dF.75, df.dF.100) %>%
  mutate(i_app = as.factor(i_app)) %>%
  filter(index_app < 57)
remove(df.dF.25, df.dF.50, df.dF.75, df.dF.100)

ggplot(data = df.dF.sweep %>% filter(id == '24_05_30_cell5_ch1', roi == '1')) +
  geom_point(aes(x = index_app, y = int, color = i_app), size = 1) +
  geom_line(aes(x = index_app, y = int, color = i_app))

plot(df.dF.raw$int[(df.dF.raw$id == '24_05_30_cell3_ch1') & (df.dF.raw$roi == 1)])
length(df.dF.raw$int[(df.dF.raw$id == '24_05_30_cell3_ch1') & (df.dF.raw$roi == 1)])
