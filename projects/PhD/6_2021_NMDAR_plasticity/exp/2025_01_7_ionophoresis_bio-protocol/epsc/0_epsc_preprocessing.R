# NMDA ionophoresis, EPSCs data preprocessing for bio-protocol
# Copyright © 2025 Borys Olifirov


require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)
require(minpack.lm)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/epsc')

base_70 <- c(0,150)
base_40 <- c(150,170)
iono_40 <- c(170, 230)
post_40 <- c(230, 250)
post_70 <- c(250, 400)
ends_70 <- c(400, 550)

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

df <- rbind(read.csv('24_10_17_3_60s_100mm.csv') %>% mutate(id = '24_10_17_3'),
            read.csv('24_10_16_11_60s_100mm.csv') %>% mutate(id = '24_10_16_11'),
            read.csv('24_10_22_6_60s_100mm.csv') %>% mutate(id = '24_10_22_6'),
            read.csv('24_10_22_8_60s_100mm.csv') %>% mutate(id = '24_10_22_8'),
            read.csv('24_10_22_10_60s_100mm.csv') %>% mutate(id = '24_10_22_10'),
            read.csv('24_10_23_2_60s_100mm.csv') %>% mutate(id = '24_10_23_2'),
            read.csv('24_10_24_7_60s_100mm.csv') %>% mutate(id = '24_10_24_7'),
            read.csv('24_10_24_17_60s_100mm.csv') %>% mutate(id = '24_10_24_17'),
            read.csv('24_10_24_19_60s_100mm.csv') %>% mutate(id = '24_10_24_19')) %>%
  mutate_at(colnames(df)[1:10], as.numeric) %>%
  mutate(id = as.factor(id),
         hold_interval = as.factor(case_when(between(TimRef, base_70[1], base_70[2]) ~ 'base_70',
                                             between(TimRef, base_40[1], base_40[2]) ~ 'base_40',
                                             between(TimRef, iono_40[1], iono_40[2]) ~ 'iono_40',
                                             between(TimRef, post_40[1], post_40[2]) ~ 'post_40',
                                             between(TimRef, post_70[1], post_70[2]) ~ 'post_70',
                                             between(TimRef, ends_70[1], ends_70[2]) ~ 'ends_70',
                                             .default = 'out')),
         time_interval = as.factor(case_when(between(TimRef, 0, 25) ~ '0',
                                             between(TimRef, 25, 50) ~ '1',
                                             between(TimRef, 50, 75) ~ '2',
                                             between(TimRef, 75, 100) ~ '3',
                                             between(TimRef, 100, 125) ~ '4',
                                             between(TimRef, 125, 150) ~ '5',
                                             between(TimRef, 150, 175) ~ '6',
                                             between(TimRef, 175, 200) ~ '7',
                                             between(TimRef, 200, 225) ~ '8',
                                             between(TimRef, 225, 250) ~ '9',
                                             between(TimRef, 250, 275) ~ '10',
                                             between(TimRef, 275, 300) ~ '11',
                                             between(TimRef, 300, 325) ~ '12',
                                             between(TimRef, 325, 350) ~ '13',
                                             between(TimRef, 350, 375) ~ '14',
                                             between(TimRef, 375, 400) ~ '15',
                                             between(TimRef, 400, 425) ~ '16',
                                             between(TimRef, 425, 450) ~ '17',
                                             between(TimRef, 450, 475) ~ '18',
                                             between(TimRef, 475, 500) ~ '19',
                                             between(TimRef, 500, 525) ~ '20',
                                             between(TimRef, 525, 550) ~ '21',
                                             .default = 'out'))) %>%
  group_by(id) %>%
  mutate(Amp_base = Amp / median(Amp[hold_interval == 'base_70']),
         Amp_norm = Amp / median(Amp[time_interval == '0']),
         Amp_delta = normalize(Amp)) %>%
  ungroup() %>%
  filter(Amp < 150)

ggplot() +
  geom_point(data = df,
             aes(x = TimPeak, y = Amp_norm, color = id), alpha = .15) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1) +
  facet_wrap(facets = vars(id), ncol = nlevels(df$id)) 
  
write.csv(df, 'epsc_olf.csv')
