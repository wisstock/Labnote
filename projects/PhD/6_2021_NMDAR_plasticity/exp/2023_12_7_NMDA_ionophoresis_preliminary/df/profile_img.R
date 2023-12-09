# NMDA ionophoresis preliminary results analysis
# Copyright © 2023 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio/note/projects/PhD/6_2021_NMDAR_plasticity/exp/2023_12_7_NMDA_ionophoresis_preliminary/df')

##### DATA PREPROCESSING #####
index.shift <- 4

df.2 <- read.csv('23_11_28_07_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = '23_11_28_07') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(index = index - index.shift) %>%
  filter(index > 0)
df.3 <- read.csv('23_11_28_07_ch3_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = '23_11_28_07') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'TagRFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(index = index - index.shift) %>%
  filter(index > 0)

df.4 <- read.csv('23_11_28_08_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = '23_11_28_08') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(index = index - index.shift) %>%
  filter(index > 0)
df.5 <- read.csv('23_11_28_08_ch3_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = '23_11_28_08') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'TagRFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(index = index - index.shift) %>%
  filter(index > 0)

df.6 <- read.csv('Fluorescence 435nm2_ch0_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = 'sh') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'CFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi))
df.7 <- read.csv('Fluorescence 505nm2_lab_prof_dF.csv') %>%
  select(-X) %>%
  mutate(id = 'sh') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi))

df.full <- bind_rows(df.2, df.3, df.4, df.5, df.6, df.7)

ggplot(data = df.full %>% filter(index < 15)) +
  geom_line(aes(x = index, y = int, color = roi, linetype = id)) +
  geom_point(aes(x = index, y = int, color = roi)) +
  facet_grid(rows = vars(hpca))

df.full <- df.full %>%
  filter(index <= 15, id == 'sh') %>%
  mutate(hpca = factor(hpca, levels = c('WT', 'N75K'))) %>%
  ungroup()

##### AVERAGE PROFILES #####
green.index <- 7
blue.index <- 12

df.avg <- df.full %>%
  ungroup %>%
  select(id, roi, int, hpca, index) %>%
  group_by(hpca, index) %>%
  mutate(int_mean = mean(int), int_se = sd(int)/sqrt(n())) %>%
  select(-int, -roi) %>%
  unique()

prof.avg <- ggplot(data = df.avg) +
  geom_line(aes(x = index, y = int_mean, color = hpca)) +
  geom_point(aes(x = index, y = int_mean, color = hpca)) +
  geom_errorbar(aes(x = index,
                    ymax = int_mean + int_se,
                    ymin = int_mean - int_se,
                    color = hpca),
                width = 0.2) +
  scale_colour_manual(name = NULL,
                    values = c('black', 'red')) +
  geom_segment(aes(x = blue.index, xend = blue.index,
                   y = 0.21, yend = 0.22),
               color = 'blue',
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = green.index, xend = green.index,
                   y = 0.025, yend = 0.035),
               color = 'green',
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = 5, linetype='dashed') +
  scale_x_continuous(name = 'Time, s',
                     limits = c(1, 15),
                     breaks = seq(1, 15, 1),
                     labels = seq(0, 74, 5)) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.05, 0.3),
                     breaks = seq(-1, 1, 0.025),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.9, .18),
        text = element_text(family='arial', face="bold", size=9))

ggsave('prof_avg_sh.png', prof.avg, 
       width = 12, height = 8.5, units = 'cm', dpi = 300)


##### GREEN BOXPLOT #####
df.green <- df.full %>%
  filter(index == green.index) %>%
  select(id, roi, int, hpca)

df.green.avg <- df.green %>%
  group_by(id, hpca) %>%
  mutate(int_mean = mean(int)) %>%
  select(-roi, -int) %>%
  unique()

df.green.stat <- df.green %>%
  select(int, hpca) %>%
  wilcox_test(int ~ hpca) %>%
  add_significance() %>%
  add_xy_position(fun = 'median_iqr', scales = 'free') %>%
  mutate(y.position = y.position)
  
box.green <- ggplot() +
  geom_boxplot(data = df.green,
               aes(x = hpca, y = int, fill = hpca), alpha = .6) +
  geom_point(data = df.green.avg,
             aes(x = hpca, y = int_mean)) +
  geom_line(data = df.green.avg,
            aes(x = hpca, y = int_mean, group = id),
            linetype = 'dashed') +
  stat_pvalue_manual(df.green.stat, label = 'p.signif',
                     hide.ns = FALSE, size = 5, family='arial') +
  scale_fill_manual(name = NULL,
                      values = c('black', 'red')) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.01, 0.35),
                     breaks = seq(0, 0.5, 0.05)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        text = element_text(family='arial', face="bold", size=8),
        panel.border = element_rect(color = "green", 
                                    fill = NA, 
                                    size = 3))
ggsave('box_green_sh.png', box.green, 
       width = 5, height = 10, units = 'cm', dpi = 300)

##### BLUE BOXPLOT #####
df.blue <- df.full %>%
  filter(index == blue.index) %>%
  select(id, roi, int, hpca)

df.blue.avg <- df.blue %>%
  group_by(id, hpca) %>%
  mutate(int_mean = mean(int)) %>%
  select(-roi, -int) %>%
  unique()

df.blue.stat <- df.blue %>%
  select(int, hpca) %>%
  wilcox_test(int ~ hpca) %>%
  add_significance() %>%
  add_xy_position(fun = 'median_iqr', scales = 'free') %>%
  mutate(y.position = y.position + 0.02)

box.blue <- ggplot() +
  geom_boxplot(data = df.blue,
               aes(x = hpca, y = int, fill = hpca), alpha = .6) +
  geom_point(data = df.blue.avg,
             aes(x = hpca, y = int_mean)) +
  geom_line(data = df.blue.avg,
            aes(x = hpca, y = int_mean, group = id),
            linetype = 'dashed') +
  stat_pvalue_manual(df.blue.stat, label = 'p.signif',
                     hide.ns = FALSE, size = 5, family='arial') +
  scale_fill_manual(name = NULL,
                    values = c('black', 'red')) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.01, 0.35),
                     breaks = seq(0, 0.5, 0.05)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        text = element_text(family='arial', face="bold", size=8),
        panel.border = element_rect(color = "blue", 
                                    fill = NA, 
                                    size = 3))
ggsave('box_blue_sh.png', box.blue, 
       width = 5, height = 10, units = 'cm', dpi = 300)
