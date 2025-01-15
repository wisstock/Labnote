# NMDA ionophoresis, FRET between EYFP-mem and HPCA(WT)-ECFP
# Copyright © 2024 Borys Olifirov

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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_08_29_AP2_FRET_60|60|180_stat')

font.size <- 17
font.fam <- 'Arial'

base.indexes <- seq(0,5)
mid.indexes <- seq(10,12) 
end.indexes <- seq(24,29)


##### DATA PREPROCESSING #####
df.mask <- bind_rows(read.csv('./fret_eyfp-mem/24_03_14_07_Eapp_24_03_14_07_Eapp_norm_red-green_up-labels_abs.csv'),  # ΔF
                     read.csv('./fret_eyfp-mem/24_03_14_06_Eapp_24_03_14_06_Eapp_norm_red-green_up-labels_abs.csv'),
                     read.csv('./fret_eyfp-mem/24_03_14_04_Eapp_24_03_14_04_Fc_norm_red-green_up-labels_abs.csv')) %>%
  mutate(id = as.factor(str_remove(id, '_xform_Eapp')),
         roi = as.factor(roi),
         cell_id = id,
         rel_time = time - 60) %>%
  filter(index != 29) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  select(-X)

df.df <- bind_rows(read.csv('./fret_eyfp-mem/24_03_14_07_Eapp_24_03_14_07_Eapp_norm_red-green_up-labels_ΔF.csv'),
                     read.csv('./fret_eyfp-mem/24_03_14_06_Eapp_24_03_14_06_Eapp_norm_red-green_up-labels_ΔF.csv'),
                     read.csv('./fret_eyfp-mem/24_03_14_04_Eapp_24_03_14_04_Fc_norm_red-green_up-labels_ΔF.csv')) %>%
  mutate(id = as.factor(str_remove(id, '_xform_Eapp')),
         roi = as.factor(roi),
         cell_id = id,
         rel_time = time - 60) %>%
  filter(index != 29) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  select(-X)

##### PROFILES #####

plot.abs <- ggplot(data = df.mask,
       aes(x = rel_time, y = int)) +
  annotate('rect', xmin = 0, xmax = 60, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black') +
  stat_summary(aes(group = roi_id, color = cell_id),
               fun = median,
               geom = 'line', size = .3, alpha = .5) +
  stat_summary(aes(group = roi_id, color = cell_id),
               fun = median,
               geom = 'point', size = 0.75, alpha = .5)  +
  stat_summary(aes(group = cell_id),
               color = 'grey15',
               fun = median,
               geom = 'point', size = 1.5) +
  stat_summary(aes(group = cell_id),
               color = 'grey15',
               fun = median,
               geom = 'line', size = 0.8) +
  stat_summary(aes(group = cell_id),
               color = 'grey15',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.6, width = 2.5) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam)) +
  scale_x_continuous(breaks = seq(-60, 230, 30)) +
  labs(x = 'Time, s',
       y = expression(E[app]))

plot.abs
setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/ca_dyn')
save_plot('0_pic_abs_fret.png', plot.abs, base_width = 6, base_height = 2.4, dpi = 300)


plot.df <- ggplot(data = df.df,
       aes(x = rel_time, y = int)) +
  annotate('rect', xmin = 0, xmax = 60, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black') +
  annotate('rect',
           xmin = (base.indexes[1]*10)-60, xmax = (rev(base.indexes)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', size = 0) +
  annotate('rect',
           xmin = (mid.indexes[1]*10)-60, xmax = (rev(mid.indexes)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', size = 0) +
  annotate('rect',
           xmin = (end.indexes[1]*10)-60, xmax = (rev(end.indexes)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', size = 0) +
  annotate("text", x=-35, y=0.6, label= "I",
           color = 'black', size = font.size-12) +
  annotate("text", x=50, y=0.6, label= "II",
           color = 'black', size = font.size-12) +
  annotate("text", x=205, y=0.6, label= "III",
           color = 'black', size = font.size-12) +
  geom_hline(yintercept = 0, lty = 2) +
  stat_summary(aes(group = roi_id, color = cell_id),
               fun = median,
               geom = 'line', size = .3, alpha = .5) +
  stat_summary(aes(group = roi_id, color = cell_id),
               fun = median,
               geom = 'point', size = 0.75, alpha = .5)  +
  stat_summary(color = 'grey15',
               fun = median,
               geom = 'point', size = 1.5) +
  stat_summary(color = 'grey15',
               fun = median,
               geom = 'line', size = 0.8) +
  stat_summary(color = 'grey15',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.6, width = 2.5) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4)) +
  scale_x_continuous(breaks = seq(-60, 230, 30)) +
  labs(caption = 'n = 1/3/115 (cultures/cells/ROIs)',
       x = 'Time, s',
       y = expression(ΔE[app]/E[app[0]]))

plot.df
setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/ca_dyn')
save_plot('0_pic_df_fret.png', plot.df, base_width = 6, base_height = 2.65, dpi = 300)

##### BOXPLOT #####
df.for.box <- df.df %>%
  select(-dist, -time, -rel_time) %>%
  group_by(roi_id) %>%
  mutate(time_interval = case_when(index %in% base.indexes ~ 'I',
                                   index %in% mid.indexes ~ 'II',
                                   index %in% end.indexes ~ 'III',
                                   .default = 'out')) %>%
  filter(time_interval != 'out') %>%
  droplevels() %>%
  ungroup() %>%
  mutate(time_interval = factor(time_interval, c('I', 'II', 'III'), ordered = TRUE)) %>%
  group_by(roi_id, time_interval) %>%
  mutate(int_interval = median(int)) %>%
  select(-index, -int) %>%
  ungroup() %>%
  distinct()

df.box.stat <- df.for.box %>%
  pairwise_wilcox_test(int_interval ~ time_interval, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = 'max') %>%
  mutate(y.position = c(0.6, 0, 0.65))


plot.box <- ggplot(data = df.for.box,
       aes(x = time_interval, y = int_interval)) +
  geom_point(color = 'grey35', size = 1.5) +
  geom_boxplot(fill = 'black', alpha = .3) +
  geom_hline(yintercept = 0, lty = 2) +
  stat_summary(aes(group = cell_id),
               fun = median,
               geom = 'line', size = 0.75, linetype = 'dashed', color = 'grey25') +
  stat_summary(aes(group = cell_id),
               fun = median,
               geom = 'point', size = 1.5, color = 'grey25') +
  stat_pvalue_manual(df.box.stat, label = 'p.adj.signif',
                     size = font.size - 10, hide.ns = TRUE, tip.length = 0.01) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam)) +
  labs(x = 'Time interval',
       y = expression(ΔE[app]/E[app[0]]))

plot.box
setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/ca_dyn')
save_plot('0_pic_box_fret.png', plot.box, base_width = 2.75, base_height = 4, dpi = 300)

##### FIN PLOTS #####

# prof
draw.abs <- ggdraw(plot.abs) +
  draw_plot_label(c("F"),
                  c(0.075),
                  c(1),
                  size = font.size + 3)

draw.df <- ggdraw(plot.df) +
  draw_plot_label(c("G"),
                  c(0.075),
                  c(1),
                  size = font.size + 3)

# box
draw.box <- ggdraw(plot.box) +
  draw_plot_label(c("H"),
                  c(0),
                  c(1),
                  size = font.size + 3)

# prod
left <- plot_grid(draw.abs, draw.df, ncol = 1)

draw.fin <- plot_grid(left, draw.box, rel_widths = c(1,0.3))
draw.fin

save_plot('0_plot_eyfp-mem.png', draw.fin, base_width = 14, base_height = 4)
