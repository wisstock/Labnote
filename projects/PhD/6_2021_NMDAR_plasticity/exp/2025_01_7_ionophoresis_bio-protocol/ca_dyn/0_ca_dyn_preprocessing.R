# NMDA ionophoresis, ionophoresis ca dynamic data preprocessing for bio-protocol
# Copyright © 2025 Borys Olifirov

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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/ca_dyn')

df.max <- rbind(read.csv('24_05_08_cell7_ch0_24_05_08_cell7_bottom_cell_max_mask_ΔF.csv'),
                read.csv('24_05_08_cell7_ch0_24_05_08_cell7_up_cell_max_mask_ΔF.csv') %>%
                  mutate(id = '24_05_08_cell7_ch0_up'),
                read.csv('24_05_9_cell01_ch0_24_05_9_cell01_max_mask_ΔF.csv'),
                read.csv('24_05_9_cell02_ch0_24_05_9_cell02_ch0_max_mask_ΔF.csv')) %>%
  mutate(roi_type = 'Max')

df.mid <- rbind(read.csv('24_05_08_cell7_ch0_24_05_08_cell7_bottom_cell_mid_mask_ΔF.csv'),
                read.csv('24_05_08_cell7_ch0_24_05_08_cell7_up_cell_mid_mask_ΔF.csv') %>%
                  mutate(id = '24_05_08_cell7_ch0_up'),
                read.csv('24_05_9_cell01_ch0_24_05_9_cell01_mid_mask_ΔF.csv'),
                read.csv('24_05_9_cell02_ch0_24_05_9_cell02_ch0_mid_mask_ΔF.csv')) %>%
  mutate(roi_type = 'Mid')

df.min <- rbind(read.csv('24_05_08_cell7_ch0_24_05_08_cell7_bottom_cell_min_mask_ΔF.csv'),
                read.csv('24_05_08_cell7_ch0_24_05_08_cell7_up_cell_min_mask_ΔF.csv') %>%
                  mutate(id = '24_05_08_cell7_ch0_up'),
                read.csv('24_05_9_cell01_ch0_24_05_9_cell01_min_mask_ΔF.csv'),
                read.csv('24_05_9_cell02_ch0_24_05_9_cell02_ch0_min_mask_ΔF.csv')) %>%
  mutate(roi_type = 'Min')

df <- rbind(df.max, df.mid, df.min) %>%
  mutate(roi = as.factor(roi)) %>%
  select(-X) %>%
  filter(index <= 60)
remove(df.max, df.mid, df.min)

ggplot() +
  geom_line(data = df %>% group_by(id, index, roi) %>% mutate(int_mean = mean(int)) %>% distinct(),
            aes(x = index-5, y = int_mean, color = id)) +
  facet_wrap(~roi, ncol = 3) +
  annotate('rect', xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black')

write.csv(df, '0_df_ca_dyn.csv')
