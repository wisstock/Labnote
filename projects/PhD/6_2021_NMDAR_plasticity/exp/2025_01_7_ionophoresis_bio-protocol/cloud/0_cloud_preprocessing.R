# NMDA ionophoresis, ionophoresis cloud data preprocessing for bio-protocol
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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/cloud')

##### ΔF DATA PREPROCESSING #####
# df.abs.raw <- rbind(read.csv('cloud/24_05_30_cell1_ch1_zone_mask_abs.csv'),
#                     read.csv('cloud/24_05_30_cell2_ch1_zone_mask_abs.csv'),
#                     read.csv('cloud/24_05_30_cell3_ch1_zone_mask_abs.csv'),
#                     read.csv('cloud/24_05_30_cell4_ch1_zone_mask_abs.csv'),
#                     read.csv('cloud/24_05_30_cell5_ch1_zone_mask_abs.csv'))

baseline.frames <- 5
start.idx.1st <- 5
start.idx.2nd <- 57
start.idx.3d <- 109
start.idx.4th <- 162

df.dF.raw.cell2_3 <- rbind(read.csv('24_05_30_cell2_ch1_zone_mask_abs.csv'),
                   read.csv('24_05_30_cell3_ch1_zone_mask_abs.csv')) %>%
             select(-X)

df.dF.25 <- df.dF.raw.cell2_3 %>%
  filter(index %in% seq(start.idx.1st - baseline.frames,start.idx.2nd)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '25') %>%
  ungroup()
df.dF.50 <- df.dF.raw.cell2_3 %>%
  filter(index %in% seq(start.idx.2nd - baseline.frames,start.idx.3d)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '50') %>%
  ungroup()
df.dF.75 <- df.dF.raw.cell2_3 %>%
  filter(index %in% seq(start.idx.3d - baseline.frames,start.idx.4th-1)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '75') %>%
  ungroup()
df.dF.100 <- df.dF.raw.cell2_3 %>%
  filter(index %in% seq(start.idx.4th - baseline.frames,214)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '100') %>%
  ungroup()

df.dF.sweep.cell2_3 <- rbind(df.dF.25, df.dF.50, df.dF.75, df.dF.100) %>%
  mutate(i_app = as.factor(i_app)) %>%
  filter(index_app < 55)
remove(df.dF.25, df.dF.50, df.dF.75, df.dF.100)

ggplot(data = df.dF.sweep.cell2_3 %>% filter(i_app == '25', roi == '1')) +
  geom_point(aes(x = index_app, y = int, color = id), size = 1) +
  geom_line(aes(x = index_app, y = int, color = id)) + 
  annotate('rect', xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black')


df.dF.raw.cell4 <- read.csv('cloud/24_05_30_cell4_ch1_zone_mask_abs.csv') %>%
                   select(-X)

df.dF.25 <- df.dF.raw.cell4 %>%
  filter(index %in% seq(6 - baseline.frames, 58)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '25') %>%
  ungroup()
df.dF.50 <- df.dF.raw.cell4 %>%
  filter(index %in% seq(66 - baseline.frames, 118)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '50') %>%
  ungroup()
df.dF.75 <- df.dF.raw.cell4 %>%
  filter(index %in% seq(126 - baseline.frames, 178)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '75') %>%
  ungroup()
df.dF.100 <- df.dF.raw.cell4 %>%
  filter(index %in% seq(186 - baseline.frames,238)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '100') %>%
  ungroup()

df.dF.sweep.cell4 <- rbind(df.dF.25, df.dF.50, df.dF.75, df.dF.100) %>%
  mutate(i_app = as.factor(i_app)) %>%
  filter(index_app < 55)
remove(df.dF.25, df.dF.50, df.dF.75, df.dF.100)

ggplot(data = df.dF.sweep.cell4 %>% filter(roi == '1')) +
  geom_point(aes(x = index_app, y = int, color = i_app), size = 1) +
  geom_line(aes(x = index_app, y = int, color = i_app)) +
  annotate('rect', xmin = 4, xmax = 24, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black')



df.dF.raw.cell5 <- read.csv('cloud/24_05_30_cell5_ch1_zone_mask_abs.csv') %>%
  select(-X)

df.dF.25 <- df.dF.raw.cell5 %>%
  filter(index %in% seq(5 - baseline.frames, 57)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '25') %>%
  ungroup()
df.dF.50 <- df.dF.raw.cell5 %>%
  filter(index %in% seq(55 - baseline.frames, 107)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '50') %>%
  ungroup()
df.dF.75 <- df.dF.raw.cell5 %>%
  filter(index %in% seq(105 - baseline.frames, 157)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '75') %>%
  ungroup()
df.dF.100 <- df.dF.raw.cell5 %>%
  filter(index %in% seq(155 - baseline.frames,207)) %>%
  select(-index, -time) %>%
  group_by(id,roi) %>%
  mutate(index_app = seq(0,57), i_app = '100') %>%
  ungroup()

df.dF.sweep.cell5 <- rbind(df.dF.25, df.dF.50, df.dF.75, df.dF.100) %>%
  mutate(i_app = as.factor(i_app)) %>%
  filter(index_app < 55)
remove(df.dF.25, df.dF.50, df.dF.75, df.dF.100)

ggplot(data = df.dF.sweep.cell5 %>% filter(roi == '1')) +
  geom_point(aes(x = index_app, y = int, color = i_app), size = 1) +
  geom_line(aes(x = index_app, y = int, color = i_app)) +
  annotate('rect', xmin = 4, xmax = 20, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black')


df.dF.sweep <- rbind(df.dF.sweep.cell2_3, df.dF.sweep.cell4, df.dF.sweep.cell5)
remove(df.dF.sweep.cell2_3, df.dF.sweep.cell4, df.dF.sweep.cell5)


write.csv(df.dF.sweep, 'cloud_sweeps_abs.csv')

####
ggplot(data = df.dF.sweep %>% filter(roi == '5', i_app == '25')) +
  geom_point(aes(x = index_app, y = int, color = id), size = 1) +
  geom_line(aes(x = index_app, y = int, color = id))

