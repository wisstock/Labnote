# NMDA ionophoresis, ionophoresis cloud data analysis for bio-protocol
# Copyright © 2025 Borys Olifirov

require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)
require(minpack.lm)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/hpca')

df <- read.csv('0_df_hpca.csv') %>%
  mutate(id = as.factor(id),
         roi = as.factor(roi),
         cell_id = as.factor(cell_id),
         app_time = as.factor(app_time)) %>%
  select(-X)

font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6


##### MAX STAT #####
df.max <- df %>%
  filter(app_time != '60') %>%
  group_by(cell_id, app_time, roi) %>%
  mutate(int_max = max(int)) %>%
  select(app_time, int_max, roi, cell_id) %>%
  distinct() %>%
  droplevels() %>%
  ungroup()

df.stat <- df.max %>%
  pairwise_wilcox_test(int_max ~ app_time, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = 'max')


ggplot(data = df.max,
       aes(x = app_time, y = int_max)) +
  geom_point(color = 'grey35', size = 1.5) +
  geom_boxplot(fill = 'black', alpha = .3) +
  stat_summary(aes(group = cell_id),
               fun = median,
               geom = 'line', size = 0.75, linetype = 'dashed', color = 'grey25') +
  stat_summary(aes(group = cell_id),
               fun = median,
               geom = 'point', size = 1.5, color = 'grey25') +
  stat_pvalue_manual(data = df.stat, hide.ns = TRUE)



boxplot_rise_tau <- ggplot(data = df.rise.fit, aes(x = roi_type, y = tau)) +
  geom_boxplot(aes(fill = roi_type), alpha = box.alpha) +
  geom_point(color = 'grey35', size = 1.5) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75, linetype = 'dashed', color = 'grey25') +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5, color = 'grey25') +
  geom_text(data = df.rise.tau.stat, aes(label = title)) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  scale_y_continuous(limits = c(0,30), breaks = seq(0,30,5)) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 2/4/35 (cultures/cells/ROIs)',
       x = 'ROI type',
       y = 'Rise \u2CA7, s')

boxplot_rise_tau
save_plot('0_pic_boxplot_rise_tau.png', boxplot_rise_tau, base_width = 3.25, base_height = 4, dpi = 300)
remove(boxplot_rise_tau)
