# NMDA ionophoresis, FRET between AP2B1-EYFP and HPCA(WT)-ECFP, ROIs profiles for representative spine
# Copyright © 2024 Borys Olifirov, Olexandra Fedchenko

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

base.indexes <- seq(-60,-10)
mid.indexes.f <- seq(40,60)
end.indexes <- seq(190,230)


df.ch0 <- read.csv('./fret_sites/24_05_22_07_ch0_24_05_22_07_representative_ΔF.csv') %>%
  mutate(roi = as.factor(roi),
         rel_time = time - 60)
df.ch3 <- read.csv('./fret_sites/24_05_22_07_ch3_24_05_22_07_representative_ΔF.csv') %>%
  mutate(roi = as.factor(roi),
         rel_time = time - 60)
df.eapp <- read.csv('./fret_sites/24_05_22_07_Eapp_24_05_22_07_representative_abs.csv') %>%
  mutate(roi = as.factor(roi),
         rel_time = time - 60)

plot.eapp <- ggplot(data = df.eapp,
       aes(x = rel_time, y = int, color = roi)) +
  geom_hline(yintercept = 0, lty = 2) +
  annotate('rect',
           xmin = 0, xmax = 60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red', size = 0) +
  annotate('text', label = 'Baseline', x = -35, y = 0.028, color = 'black', size = 5) +
  annotate('rect',
           xmin = base.indexes[1], xmax = rev(base.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  scale_color_manual(values = c('1' = 'green2',
                                '2' = 'green4',
                                '3' = 'magenta2',
                                '4' = 'red',
                                '5' = 'green3',
                                '6' = 'magenta3',
                                '7' = 'violet')) +
  labs(x = 'Time, s',
       y = expression(E[app])) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size=font.size, family=font.fam, face="bold"))


plot.ch0 <- ggplot(data = df.ch0,
       aes(x = rel_time, y = int, color = roi)) +
  geom_hline(yintercept = 0, lty = 2) +
  annotate('rect',
           xmin = 0, xmax = 60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red', size = 0) +
  annotate('text', label = 'Baseline', x = -35, y = 0.5, color = 'black', size = 5) +
  annotate('rect',
           xmin = base.indexes[1], xmax = rev(base.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  # annotate('text', label = 'II', x = 50, y = 0.53, color = 'black', size = 9) +
  # annotate('rect',
  #          xmin = mid.indexes.f[1], xmax = rev(mid.indexes.f)[1],
  #          ymin = -Inf, ymax = Inf,
  #          alpha = 0.075, color = 'black', size = 0) +
  # annotate('text', label = 'III', x = 210, y = 0.53, color = 'black', size = 9) +
  # annotate('rect',
  #          xmin = end.indexes[1], xmax = rev(end.indexes)[1],
  #          ymin = -Inf, ymax = Inf,
  #          alpha = 0.075, color = 'black', size = 0) +
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  scale_color_manual(values = c('1' = 'green2',
                                '2' = 'green4',
                                '3' = 'magenta2',
                                '4' = 'red',
                                '5' = 'green3',
                                '6' = 'magenta3',
                                '7' = 'violet')) +
  labs(x = 'Time, s',
       y = expression(ΔF/F[0])) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size=font.size, family=font.fam, face="bold"))

plot.ch3 <- ggplot(data = df.ch3,
       aes(x = rel_time, y = int, color = roi)) +
  geom_hline(yintercept = 0, lty = 2) +
  annotate('rect',
           xmin = 0, xmax = 60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red', size = 0) +
  annotate('text', label = 'Baseline', x = -35, y = 0.5, color = 'black', size = 5) +
  annotate('rect',
           xmin = base.indexes[1], xmax = rev(base.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  scale_color_manual(values = c('1' = 'green2',
                                '2' = 'green4',
                                '3' = 'magenta2',
                                '4' = 'red',
                                '5' = 'green3',
                                '6' = 'magenta3',
                                '7' = 'violet')) +
  labs(x = 'Time, s',
       y = expression(ΔF/F[0])) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size=font.size, family=font.fam, face="bold"))

draw.ch0 <- ggdraw(plot.ch0) +
  draw_plot_label(c("Da"),
                  c(0.08),
                  c(1),
                  size = font.size + 3)
draw.ch3 <- ggdraw(plot.ch3) +
  draw_plot_label(c("Db"),
                  c(0.08),
                  c(1),
                  size = font.size + 3)
draw.eapp <- ggdraw(plot.eapp) +
  draw_plot_label(c("Dc"),
                  c(0.08),
                  c(1),
                  size = font.size + 3)

draw.fret <- plot_grid(draw.ch0, draw.ch3, draw.eapp, ncol = 1)
draw.fret

save_plot('0_plot_fret_sites.png', draw.fret, base_width = 12, base_height = 8)
