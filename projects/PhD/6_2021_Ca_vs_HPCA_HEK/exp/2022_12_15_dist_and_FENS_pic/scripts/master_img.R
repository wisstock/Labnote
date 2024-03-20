# HPCA-TagRFP + Fluo-4 profiles data analysis, img for FENS forum 2022
# Copyright © 2022 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/astria/bio/institute/my_pub/conf/2022_FENS/HPCA/poster/scripts')

##### GLOBAL PLOT OPTIONS #####
font_family <- "ubuntu mono"
font_size <- 30

line_size <- 1
line_size_light <- 1.5


##### DATA PREPROCESSING #####
# profile section
df.total <- read.csv('master_px_ca.csv') %>%
  mutate(ID = as.factor(ID),
         frame = as.factor(frame))


ca.blue <- c('#0000FF', '#6767FF', '#00AEFF', '#366BFF')

ca.frame <- ggplot(df.total,
       aes(x = delta, fill = frame)) +
  geom_density(alpha = .5) +
  scale_fill_manual(name = NULL,
                    values = ca.blue,
                    labels = c('Native','After 1st stimul','After 3d stimul')) +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     limits = c(-0.2, 3),
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'Density',
                     breaks = c(),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.7, 1),
        legend.justification = c("left", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))


ggsave('ca_frame.png', ca.frame, 
       width = 32.5, height = 17.7, units = 'cm', dpi = 300)

  