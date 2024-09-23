# NMDA ionophoresis, FRET between AP2B1-EYFP and HPCA(WT)-ECFP
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

df.mask <- read.csv('df_full.csv')

##### SHAFT PROFILES #####

selected.mask <- 'up'
base.indexes <- seq(0,5)
mid.indexes.f <- seq(10,12)
mid.indexes.h <- seq(10,12)  # seq(7,9)
end.indexes <- seq(24,28)
# base.indexes <- c(0,1,2,3,4,5,6)
# mid.indexes <- c(12,13,14,15)
# end.indexes <- c(26,27,28)

# selected frames
df.fret.plot <- df.mask %>%
  filter(channel == 'Eapp',
         int_val == 'abs',
         mask == selected.mask) %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  select(roi_id, int, index, time, mask, cell_id)

ggplot(data = df.fret.plot,
       aes(x = index, y = int, color = roi_id, group = roi_id)) +
  stat_summary(fun = median,
               geom = 'line', size = .5) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'point', size = 0.75) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.15, width = 0.75) +
  annotate('rect', xmin = 6, xmax = 12, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  annotate('rect',
           xmin = base.indexes[1], xmax = rev(base.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = mid.indexes.f[1], xmax = rev(mid.indexes.f)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = end.indexes[1], xmax = rev(end.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  labs(title = 'FRET profiles in individual ROIs, shaft',
       caption = 'Red rect - NMDA app., black rect - time intervals for analysis',
       y = expression(E[app]),
       x = 'Frame idx') +
  scale_fill_manual(values = rainbow(length(df.fret.plot$roi_id))) +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(~cell_id)

df.hpca.plot <- df.mask %>%
  filter(channel == 'ch0',
         int_val == 'df',
         mask == selected.mask) %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  select(roi_id, int, index, time, mask, cell_id)

ggplot(data = df.hpca.plot,
       aes(x = index, y = int, color = roi_id, group = roi_id)) +
  stat_summary(fun = median,
               geom = 'line', size = .5) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.hpca.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'point', size = 0.75) +
  stat_summary(data = df.hpca.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(data = df.hpca.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.15, width = 0.75) +
  annotate('rect', xmin = 6, xmax = 12, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  annotate('rect',
           xmin = base.indexes[1], xmax = rev(base.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = mid.indexes.h[1], xmax = rev(mid.indexes.h)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = end.indexes[1], xmax = rev(end.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  labs(title = 'HPCA insertions profiles in individual ROIs, shaft',
       caption = 'Red rect - NMDA app., black rect - time intervals for analysis, shaft',
       y = expression(ΔF/F[0]),
       x = 'Frame idx') +
  scale_fill_manual(values = rainbow(length(df.fret.plot$roi_id))) +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(~cell_id)
