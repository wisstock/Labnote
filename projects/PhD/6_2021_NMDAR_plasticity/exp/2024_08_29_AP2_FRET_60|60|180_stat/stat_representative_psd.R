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

df.spine <- read.csv('./psd_sasha/psd_representative/24_06_06_cell1_ch0_24_06_06_cell1_spine1_ΔF.csv') %>%
  mutate(roi = as.factor(roi),
         rel_time = time - 25)

plot.spine <- ggplot(data = df.spine,
       aes(x = rel_time, y = int, color = roi)) +
  geom_hline(yintercept = 0, lty = 2) +
  # geom_segment(aes(x = 0, xend = 15,
  #                  y = 0.47, yend = 0.47),
  #              color = 'black', size = 2) +
  annotate('rect',
           xmin = 0, xmax = 15,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red', size = 0) +
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('4' = 'red',
                                '1' = 'magenta2',
                                '3' = 'blue2',
                                '5' = 'yellow3',
                                '2' = 'green4')) +
  labs(x = 'Time, s',
       y = expression(ΔF/F[0])) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size=font.size, family=font.fam, face="bold"))

draw.spine <- ggdraw(plot.spine) + draw_plot_label(c("D"),
                                  c(0),
                                  c(1),
                                  size = font.size + 5)

save_plot('0_plot_spine.png', draw.spine, base_width = 10, base_height = 4)
