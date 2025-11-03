# Sasha's data analysis, HPCA+AP2 FRET

library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(ggsci)
library(introdataviz)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_11_03_AP2_FRET_Bao')


df.full <- read.csv('df_FRET.csv') %>%
           select(-'...1') %>%
           mutate_if(is.character, factor)
           mutate(roi = as.factor(roi),
                  app_factor = as.factor(app)) %>%
           filter(app_factor %in% c('0.5','2.5', '5', '10', '20', '30', '60')) %>%
           mutate(rel_time = index - 30)


spines_id <- df.full %>%
  filter(lab_id != 'shaft') %>%
  droplevels() %>%
  mutate(roi_id = interaction(roi, id, sep = '_')) %>%
  group_by(roi_id) %>%
  summarise(has_all = n_distinct(lab_id) == 2) %>%
  filter(has_all) %>%
  droplevels() %>%
  pull(roi_id)

df.spines <- df.full %>%
  filter(lab_id != 'shaft') %>%
  mutate(roi_id = interaction(roi, id, sep = '_')) %>%
  filter(roi_id %in% spines_id) %>%
  droplevels() %>%
  group_by(roi_id) %>%
  mutate(rise_group =  ifelse(df[lab_id == 'psd'][3] > df[lab_id == 'oreol'][3], TRUE, FALSE),
         rise_group = as.factor(ifelse(!all(rise_group == FALSE), 'invert', 'direct')),
         psd_dist = if_else(lab_id == 'oreol', dist[lab_id == 'psd'][1], dist),
         psd_um = psd_dist * 0.16) %>%
  droplevels() %>%
  select(-index, -time, -app, -dist) %>%
  ungroup()

# df.spines.wide <- df.spines %>%
#   select(-abs_int, -dF_int) %>%
#   pivot_wider(names_from = lab_id, values_from = df) %>%
#   drop_na() %>%
#   group_by(roi_id) %>%
#   mutate(filter_group = ifelse(((rel_time == 2) & (oreol > psd)), TRUE, FALSE),
#          rise_group = ifelse(!all(filter_group == FALSE), 'up', 'down')) %>%
#   select(-filter_group) %>%
#   droplevels() %>%
#   ungroup()


##### EXPLORATORY ANALYSIS #####
df.spines.summary <- df.spines %>%
  group_by(lab_id, app_factor, psd_um_group) %>%
  summarise(n_cell = n_distinct(id), n_roi = n_distinct(roi),
            max = max(df))
remove(df.spines.summary)


##### ONE FRAME DF vs DIST #####
ggplot(data = df.spines %>% filter(ch == 'fc',
                                   rel_time %in% c(0,1,2,3),
                                   app_factor %in% c('5', '10', '20', '30')),
       aes(x = psd_um, y = df, color = lab_id)) +
  # geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = c(20, 40, 60), linetype = 2) +
  geom_line(aes(group = roi_id), linewidth = 0.01, linetype = 2, color = 'black') +
  geom_point(alpha = .25, size = 0.75) +
  geom_smooth(method = 'loess', span = 0.8) +
  # geom_smooth(method = 'lm', se = FALSE) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(-0.2, 0.3)) +
  scale_color_manual(values = c('darkgreen', 'magenta3')) +
  scale_fill_manual(values = c('darkgreen', 'magenta3')) +
  facet_wrap(~rel_time, ncol = 4)


###### SELECTED DIST vs FRAME #####
dist.1 <- c(0,20)    # 10
dist.2 <- c(20,40)   # 15
dist.3 <- c(40,60)   # 20
dist.4 <- c(60,100)  # 30

time.frames <- c(-1,1,3,5)

df.sel.dist <- df.spines %>%
  mutate(dist_seg = as.factor(case_when((psd_um > dist.1[1]) & (psd_um <= dist.1[2]) ~ '0-20 um',
                                        (psd_um > dist.2[1]) & (psd_um <= dist.2[2]) ~ '20-40 um',
                                        (psd_um > dist.3[1]) & (psd_um <= dist.3[2]) ~ '40-60 um',
                                        (psd_um > dist.4[1]) & (psd_um <= dist.4[2]) ~ '60-100 um',
                                        .default = '0')),
         fac_time = as.factor(rel_time)) %>%
  filter(dist_seg != '0',
         fac_time %in% time.frames) %>%
  droplevels()


df.sel.dist.stat <- df.sel.dist %>%
  filter(ch == 'fc') %>%
  group_by(fac_time, dist_seg) %>%
  wilcox_test(dF_int~lab_id) %>%
  add_significance() %>%
  add_y_position(fun = 'max') %>%
  mutate(lab_id = as.factor('oreol'),
         y.position = y.position + y.position*0.1)

ggplot(data = df.sel.dist,
       aes(x = dist_seg, y = dF_int, fill = lab_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_split_violin(alpha = 0.3, trim = FALSE, scale = "width") +
  geom_boxplot(alpha = 0.75, width = 0.2) +
  stat_summary(aes(group = lab_id, color = lab_id),
               fun = median,
               geom = 'line', size = .75,
               position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = lab_id, color = lab_id),
               fun = median,
               geom = 'point', size = 1,
               position = position_dodge(width = 0.2)) +
  geom_text(data = df.sel.dist.stat,
            aes(x = dist_seg, y = y.position, label = p.signif), size = 5) +
  facet_wrap(~fac_time, ncol = 4) +
  scale_fill_manual(values = c('darkgreen', 'magenta3')) +
  scale_color_manual(values = c('darkgreen', 'magenta3')) +
  theme_classic()


df.sel.dist.roi.stat <- df.sel.dist %>%
  group_by(lab_id, dist_seg) %>%
  pairwise_wilcox_test(df~fac_time, p.adjust.method = 'BH', ref.group = '-1') %>%
  add_significance() %>%
  add_xy_position()

ggplot(data = df.sel.dist,
       aes(x = fac_time, y = df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(aes(fill = fac_time)) +
  facet_grid(lab_id~dist_seg) +
  stat_pvalue_manual(data = df.sel.dist.roi.stat, label = 'p.adj.signif') +
  theme_classic()

regeom_point()remove(df.sel.dist)


##### OLD ROI VS DIST #####
# df.sel.dist <- df.spines %>%
#   mutate(dist_seg = as.factor(case_when(psd_um %in% dist.1 ~ 'd1',
#                                         psd_um %in% dist.2 ~ 'd2',
#                                         psd_um %in% dist.3 ~ 'd3',
#                                         psd_um %in% dist.4 ~ 'd4',
#                                         .default = '0'))) %>%
#   filter(dist_seg != '0')
#   group_by(dist_seg, as.factor(rel_time), lab_id) %>%
#   mutate(dist_df = median(df)) %>%
#   ungroup() %>%
#   select(dist_df, dist_seg, lab_id, rel_time, df)
# 
# ggplot(data = df.sel.dist %>%
#               filter(rel_time %in% seq(-5,5)),
#        aes(x = rel_time, y = df, color = dist_seg, linetype = lab_id)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point() +
#   geom_line()
