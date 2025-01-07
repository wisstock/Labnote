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

font.size <- 20
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
           alpha = 0.15, fill = 'red') +
  stat_summary(aes(group = roi_id, color = cell_id),
               fun = median,
               geom = 'line', size = .5, alpha = .5) +
  stat_summary(aes(group = roi_id, color = cell_id),
               fun = median,
               geom = 'point', size = 1, alpha = .5)  +
  stat_summary(aes(group = cell_id),
               color = 'black',
               fun = median,
               geom = 'point', size = 2) +
  stat_summary(aes(group = cell_id),
               color = 'black',
               fun = median,
               geom = 'line', size = 1) +
  stat_summary(aes(group = cell_id),
               color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.75, width = 3) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  labs(x = 'Time, s',
       y = expression(E[app]))


plot.df <- ggplot(data = df.df,
       aes(x = rel_time, y = int)) +
  annotate('rect', xmin = 0, xmax = 60, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = 'red') +
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
               geom = 'line', size = .5, alpha = .5) +
  stat_summary(aes(group = roi_id, color = cell_id),
               fun = median,
               geom = 'point', size = 1, alpha = .5)  +
  stat_summary(color = 'black',
               fun = median,
               geom = 'point', size = 2) +
  stat_summary(color = 'black',
               fun = median,
               geom = 'line', size = 1) +
  stat_summary(color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.75, width = 3) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  labs(caption = 'n = 1/3/115 (cultures/cells/ROIs)',
       x = 'Time, s',
       y = expression(ΔE[app]/E[app[0]]))


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
  add_xy_position(fun = "max") 


plot.box <- ggplot(data = df.for.box,
       aes(x = time_interval, y = int_interval)) +
  geom_boxplot(fill = 'black', alpha = .3) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(size=2, shape=21) +
  geom_line(aes(group = roi_id),
            size = .25, lty = 2, alpha = .25) +
  stat_pvalue_manual(df.box.stat, label = 'p.adj.signif',
                     hide.ns = TRUE, size = font.size - 12) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  labs(x = 'Time interval',
       y = expression(ΔE[app]/E[app[0]]))



##### FIN PLOTS #####

# prof
draw.abs <- ggdraw(plot.abs) +
  draw_plot_label(c("C"),
                  c(0.075),
                  c(1),
                  size = font.size + 3)

draw.df <- ggdraw(plot.df) +
  draw_plot_label(c("D"),
                  c(0.075),
                  c(1),
                  size = font.size + 3)

# box
draw.box <- ggdraw(plot.box) +
  draw_plot_label(c("E"),
                  c(0),
                  c(1),
                  size = font.size + 3)

# prod
left <- plot_grid(draw.abs, draw.df, ncol = 1)

draw.fin <- plot_grid(left, draw.box, rel_widths = c(1,0.3))
save_plot('0_plot_eyfp-mem.png', draw.fin, base_width = 14, base_height = 4)
