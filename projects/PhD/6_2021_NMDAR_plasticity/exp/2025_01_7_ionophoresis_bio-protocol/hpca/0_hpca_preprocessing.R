# NMDA ionophoresis, ionophoresis hpca translocation data preprocessing for bio-protocol
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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/hpca')

df.05 <- rbind(read.csv('03_04_24_cell2/cell2_0.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell2'),
               read.csv('03_04_24_cell3/cell3_0.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell3'),
               read.csv('03_04_24_cell4/cell4_0.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell4'),
               read.csv('03_04_24_cell5/cell5_0.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell5'),
               read.csv('01_04_24_cell3/cell3_0.5s_ch1_01.01.2024.csv') %>% mutate(cell_id = '01_04_24_cell3')) %>%
  mutate(app_time = '0.5')

df.25 <- rbind(read.csv('03_04_24_cell2/cell2_2.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell2'),
               read.csv('03_04_24_cell3/cell3_2.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell3'),
               read.csv('03_04_24_cell4/cell4_2.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell4'),
               read.csv('03_04_24_cell5/cell5_2.5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell5'),
               read.csv('01_04_24_cell3/cell3_2.5s_ch1_01.01.2024.csv') %>% mutate(cell_id = '01_04_24_cell3')) %>%
  mutate(app_time = '2.5')

df.50 <- rbind(read.csv('03_04_24_cell2/cell2_5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell2'),
               read.csv('03_04_24_cell3/cell3_5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell3'),
               read.csv('03_04_24_cell4/cell4_5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell4'),
               read.csv('03_04_24_cell5/cell5_5s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell5'),
               read.csv('01_04_24_cell3/cell3_5s_ch1_01.01.2024.csv') %>% mutate(cell_id = '01_04_24_cell3')) %>%
  mutate(app_time = '5')

df.60 <- rbind(read.csv('03_04_24_cell2/cell2_60s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell2'),
               read.csv('03_04_24_cell3/cell3_60s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell3'),
               read.csv('03_04_24_cell4/cell4_60s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell4'),
               read.csv('03_04_24_cell5/cell5_60s_ch1_03.04.2024.csv') %>% mutate(cell_id = '03_04_24_cell5')) %>%
  mutate(app_time = '60')

df <- rbind(df.05, df.25, df.50, df.60) %>%
  mutate(roi = as.factor(roi), app_time = as.factor(app_time)) %>%
  select(-X)
remove(df.05, df.25, df.50, df.60)

ggplot() +
  geom_line(data = df %>% group_by(id, index) %>% mutate(int_mean = mean(int)) %>% distinct(),
            aes(x = index-5, y = int_mean, color = app_time)) +
  facet_wrap(~id)

write.csv(df, '0_df_hpca.csv')
