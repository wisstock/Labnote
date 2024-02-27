# SEP calibraion results analysis
# Copyright © 2024 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio/note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_02_18_SEP_calibration')

##### 01.31 #####
df.speed <- read.csv('24_01_31_02_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(roi = as.factor(roi))

ggplot(data = df.speed,
       aes(x = time, y = int, color = roi)) +
  annotate('rect', xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = 'red') +
  geom_line() +
  geom_point() +
  scale_y_continuous(name = 'ΔF/F0')

##### 02.16 #####
df.01 <- read.csv('24_02_16_01_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('01'), ch = as.factor('SEP'), exp = as.factor(1),
         ΔH = 11, Δh = -4, roi = as.factor(roi))
df.02 <- read.csv('24_02_16_02_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('02'), ch = as.factor('Alexa'), exp = as.factor(1),
         ΔH = 11, Δh = -4, roi = as.factor(roi))

df.03 <- read.csv('24_02_16_03_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('03'), ch = as.factor('SEP'), exp = as.factor(2),
         ΔH = 11, Δh = -4, roi = as.factor(roi))
df.04 <- read.csv('24_02_16_04_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('04'), ch = as.factor('Alexa'), exp = as.factor(2),
         ΔH = 11, Δh = -4, roi = as.factor(roi))

df.05 <- read.csv('24_02_16_05_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('05'), ch = as.factor('Alexa'), exp = as.factor(3),
         ΔH = 5, Δh = -4, roi = as.factor(roi))
df.06 <- read.csv('24_02_16_06_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('06'), ch = as.factor('SEP'), exp = as.factor(3),
         ΔH = 5, Δh = -4, roi = as.factor(roi))

df.07 <- read.csv('24_02_16_07_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('07'), ch = as.factor('SEP'), exp = as.factor(4),
         ΔH = 5, Δh = 0, roi = as.factor(roi))
df.08 <- read.csv('24_02_16_08_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor('08'), ch = as.factor('Alexa'), exp = as.factor(4),
         ΔH = 5, Δh = 0, roi = as.factor(roi))

df.full <- bind_rows(df.01, df.02, df.03, df.04, df.05, df.06, df.07, df.08)
remove(df.01, df.02, df.03, df.04, df.05, df.06, df.07, df.08)



ggplot(data = df.full %>% filter(time <= 10, exp != '1', ch == 'SEP') %>% group_by(id) %>% mutate(time = time - 4.0),
       aes(x = time, y = int, color = exp)) +
  annotate('rect', xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = 'red') +
  annotate('segment', x = -2, xend = -2, y = -Inf, yend = Inf,
           colour = 'black') +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = vars(roi), nrow = 5, strip.position = 'right') +
  scale_y_continuous(name = 'ΔF/F0')


ggplot(data = df.full %>% filter(exp == '2'),
       aes(x = time, y = int, color = roi)) +
  annotate('rect', xmin = 2, xmax = 22, ymin = -Inf, ymax = Inf,
           alpha = 0.075, fill = 'blue') +
  annotate('rect', xmin = 4, xmax = 5, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = 'red') +
  annotate('rect', xmin = 19, xmax = 20, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = 'red') +
  geom_line() +
  geom_point() +
  facet_wrap(facets = vars(ch), scales = 'free_y', nrow = 2) +
  scale_y_continuous(name = 'ΔF/F0')
