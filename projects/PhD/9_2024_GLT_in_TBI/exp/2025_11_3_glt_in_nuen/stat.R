# GLT in neyrons 
#
# COLUMNS FROM PLUGIN
# cell_row = [glt_img.name,                        # id
#             group,                               # group
#             a_region.label,                      # cell_num
#             one_cell_area,                       # cell_area
#             one_dot_num,                         # dot_num
#             one_dot_area,                        # dot_area  це сумарна площа на один астроцит, тобто це не можна використати як параметр окремої точки
#             round(one_dot_rel_area, 3),          # dot_rel_area
#             one_dot_sum_int,                     # dot_sum_int
#             int(one_dot_mean_int),               # dot_mean_int
#             int(one_dot_sum_int / one_dot_num),  # dot_men_int_per_dot
#             int(one_dot_mean_int_dens)]          # dot_mean_int_dens


library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggsci)


setwd('/home/wisstock/bio_note/projects/PhD/9_2024_GLT_in_TBI/exp/2025_11_3_glt_in_nuen')

df <- read.csv('df_NeuN.csv') %>%
  mutate_if(is.character, as.factor) %>%
  mutate(group = factor(group, levels = c('Cont', 'TBI_3D', 'TBI_7D', 'TBI_14D')))

df.med <- df %>%
  select(group, treat, relative_area, relative_intensity) %>%
  group_by(group, treat) %>%
  mutate(med_rel_area = median(relative_area),
         med_rel_int = median(relative_intensity)) %>%
  select(-relative_area, -relative_intensity) %>%
  distinct() %>%
  ungroup()
 
##### REL AREA #####
ggplot(data = df %>% filter(group != 'Cont'),
       aes(x = group, y = relative_area,
           color = treat, group = treat)) +
  geom_hline(yintercept = 0.36, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75)

##### REL INT #####
ggplot(data = df %>% filter(group != 'Cont'),
       aes(x = group, y = relative_intensity,
           color = treat, group = treat)) +
  geom_hline(yintercept = 0.4, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75)
