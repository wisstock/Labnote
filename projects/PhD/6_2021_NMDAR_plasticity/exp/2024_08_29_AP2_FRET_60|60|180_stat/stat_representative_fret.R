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

font.size <- 25
font.fam <- 'Arial'

df.ch0 <- read.csv('./fret_spine/24_05_16_09_ch0_24_05_16_09_FRET_representative-labels_ΔF.csv') %>%
  mutate(roi = as.factor(roi),
         rel_time = time * 2 - 60)
df.ch3 <- read.csv('./fret_spine/24_05_16_09_ch3_24_05_16_09_FRET_representative-labels_ΔF.csv') %>%
  mutate(roi = as.factor(roi),
         rel_time = time * 2 - 60)
df.eapp <- read.csv('./fret_spine/24_05_16_09_Eapp_24_05_16_09_FRET_representative-labels_abs.csv') %>%
  mutate(roi = as.factor(roi),
         rel_time = time * 2 - 60)

plot.eapp <- ggplot(data = df.eapp,
       aes(x = rel_time, y = int, color = roi)) +
  geom_hline(yintercept = 0, lty = 2) +
  annotate('rect',
           xmin = 0, xmax = 60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red', size = 0) +
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('4' = 'green2',
                                '1' = 'red',
                                '3' = 'magenta',
                                '5' = 'green4',
                                '2' = 'red4')) +
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
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('4' = 'green2',
                                '1' = 'red',
                                '3' = 'magenta',
                                '5' = 'green4',
                                '2' = 'red4')) +
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
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('4' = 'green2',
                                '1' = 'red',
                                '3' = 'magenta',
                                '5' = 'green4',
                                '2' = 'red4')) +
  labs(x = 'Time, s',
       y = expression(ΔF/F[0])) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size=font.size, family=font.fam, face="bold"))

draw.ch0 <- ggdraw(plot.ch0) +
  draw_plot_label(c("E"),
                  c(0),
                  c(1),
                  size = font.size + 5)
draw.ch3 <- ggdraw(plot.ch3) +
  draw_plot_label(c("F"),
                  c(0),
                  c(1),
                  size = font.size + 5)
draw.eapp <- ggdraw(plot.eapp) +
  draw_plot_label(c("G"),
                  c(0),
                  c(1),
                  size = font.size + 5)

draw.fret <- plot_grid(draw.ch0, draw.ch3, draw.eapp, nrow = 1)

save_plot('0_plot_fret_rep.png', draw.fret, base_width = 20, base_height = 4)
