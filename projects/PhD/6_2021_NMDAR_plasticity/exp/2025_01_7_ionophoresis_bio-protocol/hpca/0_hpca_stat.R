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

color.05 <- 'coral3'
color.25 <- 'cadetblue'
color.50 <- 'purple3'

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
  pairwise_wilcox_test(int_max ~ app_time,
                       p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = 'max') %>%
  mutate(y.position = c(0.67, 0.73, 0))

boxplot_max_amp <- ggplot(data = df.max,
       aes(x = app_time, y = int_max)) +
  geom_boxplot(aes(fill = app_time), alpha = .5) +
  geom_point(color = 'grey35', size = 1.5) +
  stat_summary(aes(group = cell_id),
               fun = median,
               geom = 'line', size = 0.75, linetype = 'dashed', color = 'grey25') +
  stat_summary(aes(group = cell_id),
               fun = median,
               geom = 'point', size = 1.5, color = 'grey25') +
  stat_pvalue_manual(data = df.stat, size = font.size - 10, hide.ns = TRUE, tip.length = 0.01) +
  scale_y_continuous(limits = c(0,0.75), breaks = seq(0,30,0.25)) +
  scale_fill_manual(values = c('0.5' = color.05, '2.5' = color.25, '5' = color.50)) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 2/5/71 (cultures/cells/ROIs)',
       x = 'App. duration, s',
       y = expression(ΔF/F[0]))

boxplot_max_amp
save_plot('0_pic_boxplot_hpca_max_amp.png', boxplot_max_amp, base_width = 3.25, base_height = 4, dpi = 300)
remove(boxplot_max_amp)


##### REPRESENTATIVE PROFILES #####
df.prof <- df %>%
  filter(app_time != '60', cell_id == '03_04_24_cell4') %>%
  mutate(index = index - 5) %>%
  droplevels()

profile_hpca_rep <- ggplot(data = df.prof,
       aes(x = index, y = int, color = app_time, fill = app_time, group = app_time)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(fun = median,
               geom = 'point', size = 2) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  annotate("segment", x = -0.5, xend = 0, y = 0.38, yend = 0.38, size = 5,
           colour = color.05) +
  annotate("segment", x = -0.5, xend = 2, y = 0.41, yend = 0.41, size = 5,
           colour = color.25) +
  annotate("segment", x = -0.5, xend = 4.5, y = 0.44, yend = 0.44, size = 5,
           colour = color.50) +
  scale_fill_manual(name = "App. duration, s",
                    values = c('0.5' = color.05, '2.5' = color.25, '5' = color.50)) +
  scale_color_manual(name = "App. duration, s",
                     values = c('0.5' = color.05, '2.5' = color.25, '5' = color.50)) +
  scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15),
                     limits = c(-5, 15)) +
  theme_classic() +
  theme(legend.position = c(0.82,0.85),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 1/1/14 (cultures/cells/ROIs)',
       x = 'Time, s',
       y = expression(ΔF/F[0]))


profile_hpca_rep
save_plot('0_pic_profile_hpca_rep.png', profile_hpca_rep, base_width = 5.6, base_height = 4, dpi = 300)
remove(profile_hpca_rep)