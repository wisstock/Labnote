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


df.full <- read.csv('df_processed.csv') %>%
           select(-X) %>%
           mutate_if(is.character, factor) %>%
           mutate(app_time = as.factor(app_time),
                  roi = as.factor(roi),
                  roi_id = interaction(id, roi, sep = '_'),
                  dF_F0_int = dF.F0_int) %>%
           select(-dF.F0_int)


spines_id <- df.full %>%
  filter(lab_id != 'shaft', lab_id != 'psd_dots') %>%
  droplevels() %>%
  group_by(roi_id) %>%
  summarise(has_all = n_distinct(lab_id) == 2) %>%
  filter(has_all) %>%
  droplevels() %>%
  pull(roi_id)

df.spines <- df.full %>%
  filter(lab_id != 'shaft', lab_id != 'psd_dots') %>%
  filter(roi_id %in% spines_id) %>%
  droplevels() %>%
  group_by(roi_id) %>%
  mutate(dist_um = dist * 0.16) %>%
  droplevels() %>%
  select(-index) %>%
  ungroup()

df.spines.05 <- df.spines %>%
  filter(app_time == '0.5') %>%
  select(-app_time) %>%
  mutate(rel_time = time - 1)

df.spines.60 <- df.spines %>%
  filter(app_time == '60') %>%
  select(-app_time) %>%
  mutate(rel_time = time - 40)


##### EXPLORATORY ANALYSIS #####
df.spines.summary <- df.spines %>%
  filter(base == 'simple') %>%
  group_by(lab_id, app_time) %>%
  summarise(n_cell = n_distinct(id), n_roi = n_distinct(roi),
            max = max(abs_int),
            min = min(abs_int))
            # median = median(dF.F0_int),
            # iqr = IQR(dF.F0_int))
remove(df.spines.summary)

# PDFs
ggplot(df.spines.60 %>% filter(base == 'simple',
                               rel_time == 10),
       aes(x = dF_F0_int, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_density(alpha = .5) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free_y')

ggplot(df.spines.05 %>% filter(base == 'simple',
                               rel_time == 0.25),
       aes(x = dF_F0_int, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_density(alpha = .5) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free_y')

# 0.5 S PROF
ggplot(data = df.spines.05 %>% filter(base == 'simple'),
       aes(x = rel_time, y = dF_F0_int, color = id, fill = id)) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free_y')

# 60 S PROF
ggplot(data = df.spines.60 %>% filter(base == 'dietrich'),
       aes(x = rel_time, y = dF_int, color = id, fill = id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free_y')


# summary plot
ggplot(data = df.spines.05 %>% filter(base == 'simple'),
       aes(x = rel_time, y = dF_F0_int, color = lab_id, fill = lab_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scales = 'free_y')




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
