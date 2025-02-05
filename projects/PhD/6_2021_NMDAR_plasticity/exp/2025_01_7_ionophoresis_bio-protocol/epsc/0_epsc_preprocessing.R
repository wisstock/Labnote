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

## 100 mM
df.100 <- rbind(read.csv('24_10_17_3_60s_100mm.csv') %>% mutate(id = '24_10_17_3'),
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
         mid_interval = as.factor(case_when(between(TimRef, 0, 75) ~ 't-75',
                                            between(TimRef, 75, 150) ~ 't0',
                                            between(TimRef, 250, 325) ~ 't75',
                                            between(TimRef, 325, 400) ~ 't150',
                                            between(TimRef, 400, 475) ~ 't225',
                                            between(TimRef, 475, 550) ~ 't300',
                                            .default = '-40'))) %>%
  group_by(id) %>%
  mutate(Amp_base = Amp / median(Amp[hold_interval == 'base_70']),
         Amp_mid = Amp / median(Amp[mid_interval == 't-75'])) %>%
  ungroup() %>%
  filter(Amp < 150)

ggplot() +
  geom_point(data = df.100,
             aes(x = TimPeak, y = Amp_mid, group = id, color = mid_interval), alpha = .15) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1) +
  facet_wrap(facets = vars(id), ncol = nlevels(df.100$id)) 
  
# write.csv(df.100, 'epsc_100.csv')

## 100 mM + APV
df.100.apv <- rbind(read.csv('25_01_30_15_60s_100mm_apv.csv') %>% mutate(id = '25_01_30_15'),
            read.csv('25_01_30_18_60s_100mm_apv.csv') %>% mutate(id = '25_01_30_18'),
            read.csv('25_01_31_3_60s_100mm_apv.csv') %>% mutate(id = '25_01_31_3'),
            read.csv('25_01_31_6_60s_100mm_apv.csv') %>% mutate(id = '25_01_31_6'),
            read.csv('25_02_3_12_60s_100mm_apv.csv') %>% mutate(id = '25_02_3_12'),
            read.csv('25_02_3_15_60s_100mm_apv.csv') %>% mutate(id = '25_02_3_15'),
            read.csv('25_02_4_2_60s_100mm_apv.csv') %>% mutate(id = '25_02_4_2'),
            read.csv('25_02_4_12_60s_100mm_apv.csv') %>% mutate(id = '25_02_4_12')) %>%
  mutate_at(colnames(df)[1:10], as.numeric) %>%
  mutate(id = as.factor(id),
         hold_interval = as.factor(case_when(between(TimRef, base_70[1], base_70[2]) ~ 'base_70',
                                             between(TimRef, base_40[1], base_40[2]) ~ 'base_40',
                                             between(TimRef, iono_40[1], iono_40[2]) ~ 'iono_40',
                                             between(TimRef, post_40[1], post_40[2]) ~ 'post_40',
                                             between(TimRef, post_70[1], post_70[2]) ~ 'post_70',
                                             between(TimRef, ends_70[1], ends_70[2]) ~ 'ends_70',
                                             .default = 'out')),
         mid_interval = as.factor(case_when(between(TimRef, 0, 75) ~ 't-75',
                                            between(TimRef, 75, 150) ~ 't0',
                                            between(TimRef, 250, 325) ~ 't75',
                                            between(TimRef, 325, 400) ~ 't150',
                                            between(TimRef, 400, 475) ~ 't225',
                                            between(TimRef, 475, 550) ~ 't300',
                                            .default = '-40'))) %>%
  group_by(id) %>%
  mutate(Amp_base = Amp / median(Amp[hold_interval == 'base_70']),
         Amp_mid = Amp / median(Amp[mid_interval == 't-75'])) %>%
  ungroup() %>%
  filter(Amp < 150)

ggplot() +
  geom_point(data = df.100.apv,
             aes(x = TimPeak, y = Amp_mid, group = id, color = mid_interval), alpha = .15) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1) +
  facet_wrap(facets = vars(id), ncol = nlevels(df.100.apv$id)) 

# write.csv(df, 'epsc_olf.csv')

## 15 mM
df.15 <- rbind(read.csv('24_10_14_4_60s_15mm.csv') %>% mutate(id = '24_10_14_4'),
               read.csv('25_01_2_3_60s_15mm.csv') %>% mutate(id = '25_01_2_3'),
               read.csv('25_01_7_1_60s_15mm.csv') %>% mutate(id = '25_01_7_1'),
               read.csv('25_01_8_4_filtered_60s_15mm.csv') %>% mutate(id = '25_01_8_4'),
               read.csv('25_01_8_9_60s_15mm.csv') %>% mutate(id = '25_01_8_9'),
               read.csv('25_01_9_2_60s_15mm.csv') %>% mutate(id = '25_01_9_2')) %>%
  mutate_at(colnames(df)[1:10], as.numeric) %>%
  mutate(id = as.factor(id),
         hold_interval = as.factor(case_when(between(TimRef, base_70[1], base_70[2]) ~ 'base_70',
                                                      between(TimRef, base_40[1], base_40[2]) ~ 'base_40',
                                                      between(TimRef, iono_40[1], iono_40[2]) ~ 'iono_40',
                                                      between(TimRef, post_40[1], post_40[2]) ~ 'post_40',
                                                      between(TimRef, post_70[1], post_70[2]) ~ 'post_70',
                                                      between(TimRef, ends_70[1], ends_70[2]) ~ 'ends_70',
                                                      .default = 'out')),
                  mid_interval = as.factor(case_when(between(TimRef, 0, 75) ~ 't-75',
                                                     between(TimRef, 75, 150) ~ 't0',
                                                     between(TimRef, 250, 325) ~ 't75',
                                                     between(TimRef, 325, 400) ~ 't150',
                                                     between(TimRef, 400, 475) ~ 't225',
                                                     between(TimRef, 475, 550) ~ 't300',
                                                     .default = '-40'))) %>%
  group_by(id) %>%
  mutate(Amp_base = Amp / median(Amp[hold_interval == 'base_70']),
         Amp_mid = Amp / median(Amp[mid_interval == 't-75'])) %>%
  ungroup() %>%
  filter(Amp < 150)

ggplot() +
  geom_point(data = df.15,
             aes(x = TimPeak, y = Amp_mid, group = id, color = mid_interval), alpha = .15) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1) +
  facet_wrap(facets = vars(id), ncol = nlevels(df.15$id)) 

# write.csv(df.15, 'epsc_olf_15mm.csv')

df.output <- rbind(df.100 %>% mutate(conc = '100mM'),
                   df.100.apv %>% mutate(conc = '100mM+APV'),
                   df.15 %>% mutate(conc = '15mM')) %>%
  mutate(conc = as.factor(conc))

write.csv(df.output, 'df_epsc_raw.csv')

# mutate(id = as.factor(id),
#        hold_interval = as.factor(case_when(between(TimRef, base_70[1], base_70[2]) ~ 'base_70',
#                                            between(TimRef, base_40[1], base_40[2]) ~ 'base_40',
#                                            between(TimRef, iono_40[1], iono_40[2]) ~ 'iono_40',
#                                            between(TimRef, post_40[1], post_40[2]) ~ 'post_40',
#                                            between(TimRef, post_70[1], post_70[2]) ~ 'post_70',
#                                            between(TimRef, ends_70[1], ends_70[2]) ~ 'ends_70',
#                                            .default = 'out')),
#        mid_interval = as.factor(case_when(between(TimRef, 0, 75) ~ 't-75',
#                                           between(TimRef, 75, 150) ~ 't0',
#                                           between(TimRef, 250, 325) ~ 't75',
#                                           between(TimRef, 325, 400) ~ 't150',
#                                           between(TimRef, 400, 475) ~ 't225',
#                                           between(TimRef, 475, 550) ~ 't300',
#                                            .default = '-40')),
#        time_interval = as.factor(case_when(between(TimRef, 0, 25) ~ '0',
#                                            between(TimRef, 25, 50) ~ '1',
#                                            between(TimRef, 50, 75) ~ '2',
#                                            between(TimRef, 75, 100) ~ '3',
#                                            between(TimRef, 100, 125) ~ '4',
#                                            between(TimRef, 125, 150) ~ '5',
#                                            between(TimRef, 150, 175) ~ '6',
#                                            between(TimRef, 175, 200) ~ '7',
#                                            between(TimRef, 200, 225) ~ '8',
#                                            between(TimRef, 225, 250) ~ '9',
#                                            between(TimRef, 250, 275) ~ '10',
#                                            between(TimRef, 275, 300) ~ '11',
#                                            between(TimRef, 300, 325) ~ '12',
#                                            between(TimRef, 325, 350) ~ '13',
#                                            between(TimRef, 350, 375) ~ '14',
#                                            between(TimRef, 375, 400) ~ '15',
#                                            between(TimRef, 400, 425) ~ '16',
#                                            between(TimRef, 425, 450) ~ '17',
#                                            between(TimRef, 450, 475) ~ '18',
#                                            between(TimRef, 475, 500) ~ '19',
#                                            between(TimRef, 500, 525) ~ '20',
#                                            between(TimRef, 525, 550) ~ '21',
#                                            .default = 'out'))) %>%
# group_by(id) %>%
# mutate(Amp_base = Amp / median(Amp[hold_interval == 'base_70']),
#        Amp_norm = Amp / median(Amp[time_interval == '0']),
#        Amp_delta = normalize(Amp),
#        Amp_mid = Amp / median(Amp[mid_interval == 't-75'])) %>%
# ungroup() %>%
# filter(Amp < 150)
