# NMDA ionophoresis, EPSCs record plot for bio-protocol
# Copyright © 2025 Borys Olifirov

require(dplyr)
require(tidyr)

require(ggplot2)


setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/epsc')

font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6

base.01 <- read.csv('traces/base_01.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))
base.02 <- read.csv('traces/base_02.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))

wo.01 <- read.csv('traces/wo_01.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))
wo.02 <- read.csv('traces/wo_01.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))


df.rep <- rbind(base.01 %>% slice(1:31000) %>% mutate(time = 'base'),
                wo.02 %>% slice(1:31000) %>% mutate(time = 'wo')) %>%
  mutate(time = as.factor(time)) %>%
  group_by(time) %>%
  mutate(index = seq(1, 31000),
         pa = if_else(time == 'base', pa + 80, pa - 80))
  


prof_epsc <- ggplot(data = df.rep,
       aes(x = index, y = pa, group = time)) +
  geom_line() +
  annotate(geom = 'rect',  # 20 pA scale bar
           xmin = 0, xmax = 280,
           ymin = -44, ymax = -24,
           fill = 'black') +
  annotate(geom = 'rect',  # 10 s scale bar
           xmin = 0, xmax = 3300,
           ymin = -44, ymax = -40,
           fill = 'black') +
  theme_void() +
  theme(legend.position = 'none')

prof_epsc
save_plot('0_pic_prof_epsc.png', prof_epsc, base_width = 7, base_height = 4.4, dpi = 300)
remove(prof_epsc)