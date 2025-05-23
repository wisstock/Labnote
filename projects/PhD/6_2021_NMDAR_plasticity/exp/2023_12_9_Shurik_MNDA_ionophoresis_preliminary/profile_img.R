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

setwd('/home/wisstock/bio/note/projects/PhD/6_2021_NMDAR_plasticity/exp/2023_12_9_Shurik_MNDA_ionophoresis_preliminary')

##### DATA PREPROCESSING #####
df.start <- 2
df.end <- 5

index.shift <- 3

df.0 <- read.csv('00_435_raw.csv') %>%
  select(-X) %>%
  mutate(id = '00') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'CFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  mutate(index = index - index.shift) %>%
  filter(index > 0) %>%
  select(-f0)
df.1 <- read.csv('00_505_raw.csv') %>%
  select(-X) %>%
  mutate(id = '00') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  mutate(index = index - index.shift) %>%
  filter(index > 0) %>%
  select(-f0)

df.2 <- read.csv('02_435_raw.csv') %>%
  select(-X) %>%
  mutate(id = '02') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'CFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)
df.3 <- read.csv('02_505_raw.csv') %>%
  select(-X) %>%
  mutate(id = '02') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)

df.4 <- read.csv('03_435_raw.csv') %>%
  select(-X) %>%
  mutate(id = '03') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'CFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)
df.5 <- read.csv('03_505_raw.csv') %>%
  select(-X) %>%
  mutate(id = '03') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)

df.6 <- read.csv('04_435_raw.csv') %>%
  select(-X) %>%
  mutate(id = '04') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'CFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)
df.7 <- read.csv('04_505_raw.csv') %>%
  select(-X) %>%
  mutate(id = '04') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)

df.8 <- read.csv('05_435_raw.csv') %>%
  select(-X) %>%
  mutate(id = '05') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'CFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)
df.9 <- read.csv('05_505_raw.csv') %>%
  select(-X) %>%
  mutate(id = '05') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi)) %>%
  mutate(f0 = mean(int[index >= df.start & index < df.end])) %>%
  mutate(df = (int-f0)/ f0) %>%
  select(-f0)

df.full <- bind_rows(df.0, df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8, df.9)

df.full <- df.full %>%
  filter(index <= 8) %>%
  mutate(hpca = factor(hpca, levels = c('WT', 'N75K'))) %>%
  ungroup()

ggplot(data = df.full %>% filter(index < 15)) +
  geom_line(aes(x = index, y = df, color = roi, linetype = id)) +
  geom_point(aes(x = index, y = df, color = roi, linetype = id)) +
  facet_grid(rows = vars(hpca))

##### FULL REPRESENTATIVE PROFILE #####
df.rep.0 <- read.csv('00_435_dF.csv') %>%
  select(-X) %>%
  mutate(id = '00') %>%
  mutate(hpca = 'WT') %>%
  mutate(fp = 'CFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi))
df.rep.1 <- read.csv('00_505_dF.csv') %>%
  select(-X) %>%
  mutate(id = '00') %>%
  mutate(hpca = 'N75K') %>%
  mutate(fp = 'EYFP') %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(index = row_number(roi))

df.rep <- bind_rows(df.rep.0, df.rep.1) %>%
  select(roi, int, hpca, index, time) %>%
  group_by(hpca, index) %>%
  mutate(int_mean = mean(int), int_se = sd(int)/sqrt(n())) %>%
  mutate(hpca = factor(hpca, levels = c('WT', 'N75K'))) %>%
  mutate(time = time - 15.0) %>%
  mutate(`ΔF/F0` = int_mean)

prof.rep <- ggplot(data = df.rep) +
  geom_line(aes(x = time, y = `ΔF/F0`, color = hpca)) +
  geom_point(aes(x = time, y = `ΔF/F0`, color = hpca)) +
  geom_ribbon(aes(x = time,
                  ymax = `ΔF/F0` + int_se,
                  ymin = `ΔF/F0` - int_se,
                  color = hpca,
                  fill = hpca),
              alpha=.15,
              size=0.25) +
  geom_segment(aes(x = 15, xend = 315, y = 0.21, yend = 0.21), size = 1) +
  scale_colour_manual(name = NULL,
                      values = c('black', 'red')) +
  scale_fill_manual(name = NULL,
                      values = c('black', 'red')) +
  scale_x_continuous(name = 'Time, s',
                     limits = c(0, 400),
                     breaks = seq(0, 1000, 50)) +
  scale_y_continuous(limits = c(-0.02, 0.22),
                     breaks = seq(-1, 1, 0.05),
                     expand = c(0,0)) +
#  scale_x_break(c(75, 270)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        axis.title.x = element_text(margin = margin(t = -0.1)),
        text = element_text(family='arial', face="bold", size=9))

ggsave('prof_rep.png', prof.rep, 
       width = 10.5, height = 6, units = 'cm', dpi = 300)

##### RISE AVERAGE BOX #####
df.full.stat <- df.full %>%
  group_by(index) %>%
  pairwise_wilcox_test(df ~ hpca) %>%  # adjust_pvalue(method = "bonferroni")
  add_significance("p.adj") %>%
  add_y_position(fun = "mean_sd")  

box.avg <- ggplot() +
  geom_boxplot(data = df.full,
               aes(y = df, x = as.factor(index), fill = hpca,
                   group = interaction(index, hpca)),
               alpha = .5,
               width = 8) +
  stat_summary(data = df.full,
               fun = median,
               geom = 'line',
               aes(y = df, x = index, color = hpca,
                   group = hpca),
               position = position_dodge(width = 1)) +
  stat_summary(data = df.full,
               fun = median,
               geom = 'point',
               aes(y = df, x = index, color = hpca,
                   group = hpca),
               position = position_dodge(width = 1)) +
  geom_text(data = df.full.stat,
            aes(x = index, y = y.position, label = p.adj.signif),
            size = 3) +
  scale_colour_manual(name = NULL,
                      values = c('black', 'red')) +
  scale_fill_manual(name = NULL,
                      values = c('black', 'red')) +
  geom_vline(xintercept = 4.52, linetype='dashed', size = 0.4) +
  scale_x_discrete(name = 'Time, s',
                     breaks = seq(1, 15, 1),
                     labels = seq(0, 74, 5)) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.1, 0.17),
                     breaks = seq(-1, 1, 0.05),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.05,
                                        colour = "black"),
        legend.position = c(.3, .9),
        text = element_text(family='arial', face="bold", size=9))

ggsave('box_avg.png', box.avg, 
       width = 15, height = 7, units = 'cm', dpi = 300)

##### RISE AVERAGE PROFILES #####
magenta.index <- 3
green.index <- 5
blue.index <- 6

df.avg <- df.full %>%
  select(id, roi, df, hpca, index) %>%
  group_by(hpca, index) %>%
  mutate(int_mean = mean(df), int_se = sd(df)/sqrt(n())) %>%
  select(-df, -roi) %>%
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
  geom_segment(aes(x = magenta.index, xend = magenta.index,
                   y = -0.016, yend = -0.006),
               color = 'magenta',
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = green.index, xend = green.index,
                   y = -0.025, yend = -0.015),
               color = 'green',
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = blue.index, xend = blue.index,
                   y = 0.06, yend = 0.07),
               color = 'blue',
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = 4.95, linetype='dashed', size = 0.4) +
  scale_x_continuous(name = 'Time, s',
                     limits = c(1, 7),
                     breaks = seq(1, 15, 1),
                     labels = seq(0, 74, 5)) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.03, 0.11),
                     breaks = seq(-1, 1, 0.025),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.3, .9),
        text = element_text(family='arial', face="bold", size=9))

ggsave('prof_avg.png', prof.avg, 
       width = 4, height = 6, units = 'cm', dpi = 300)

##### MAGENTA BOXPLOT #####
df.magenta <- df.full %>%
  filter(index == magenta.index) %>%
  select(df, roi, id, hpca)

df.magenta.avg <- df.magenta %>%
  group_by(id, hpca) %>%
  mutate(int_mean = mean(df)) %>%
  select(-roi, -df) %>%
  unique()

df.magenta.stat <- df.magenta %>%
  select(df, hpca) %>%
  wilcox_test(df ~ hpca) %>%
  add_significance() %>%
  add_xy_position(fun = 'median_iqr', scales = 'free') %>%
  mutate(y.position = y.position)

box.magenta <- ggplot() +
  geom_boxplot(data = df.magenta,
               aes(x = hpca, y = df, fill = hpca), alpha = .6) +
  geom_point(data = df.magenta.avg,
             aes(x = hpca, y = int_mean)) +
  geom_line(data = df.magenta.avg,
            aes(x = hpca, y = int_mean, group = id),
            linetype = 'dashed') +
  stat_pvalue_manual(df.magenta.stat, label = 'p.signif',
                     hide.ns = FALSE, size = 3, family='arial') +
  scale_fill_manual(name = NULL,
                    values = c('black', 'red')) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.05, 0.175),
                     breaks = seq(-0.5, 0.5, 0.025)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        text = element_text(family='arial', face="bold", size=9),
        panel.border = element_rect(color = "magenta", 
                                    fill = NA, 
                                    size = 2))

ggsave('box_magenta.png', box.magenta, 
       width = 3, height = 6, units = 'cm', dpi = 300)

##### GREEN BOXPLOT #####
df.green <- df.full %>%
  filter(index == green.index) %>%
  select(df, roi, id, hpca)

df.green.avg <- df.green %>%
  group_by(id, hpca) %>%
  mutate(int_mean = mean(df)) %>%
  select(-roi, -df) %>%
  unique()

df.green.stat <- df.green %>%
  select(df, hpca) %>%
  wilcox_test(df ~ hpca) %>%
  add_significance() %>%
  add_xy_position(fun = 'median_iqr', scales = 'free') %>%
  mutate(y.position = y.position)
  
box.green <- ggplot() +
  geom_boxplot(data = df.green,
               aes(x = hpca, y = df, fill = hpca), alpha = .6) +
  geom_point(data = df.green.avg,
             aes(x = hpca, y = int_mean)) +
  geom_line(data = df.green.avg,
            aes(x = hpca, y = int_mean, group = id),
            linetype = 'dashed') +
  stat_pvalue_manual(df.green.stat, label = 'p.signif',
                     hide.ns = FALSE, size = 3, family='arial') +
  scale_fill_manual(name = NULL,
                      values = c('black', 'red')) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.05, 0.175),
                     breaks = seq(-0.5, 0.5, 0.025)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        text = element_text(family='arial', face="bold", size=9),
        panel.border = element_rect(color = "green", 
                                    fill = NA, 
                                    size = 2))
ggsave('box_green.png', box.green, 
       width = 3, height = 6, units = 'cm', dpi = 300)

##### BLUE BOXPLOT #####
df.blue <- df.full %>%
  filter(index == blue.index) %>%
  select(id, roi, df, hpca)

df.blue.avg <- df.blue %>%
  group_by(id, hpca) %>%
  mutate(int_mean = mean(df)) %>%
  select(-roi, -df) %>%
  unique()

df.blue.stat <- df.blue %>%
  select(df, hpca) %>%
  wilcox_test(df ~ hpca) %>%
  add_significance() %>%
  add_xy_position(fun = 'median_iqr', scales = 'free') %>%
  mutate(y.position = y.position + 0.02)

box.blue <- ggplot() +
  geom_boxplot(data = df.blue,
               aes(x = hpca, y = df, fill = hpca), alpha = .6) +
  geom_point(data = df.blue.avg,
             aes(x = hpca, y = int_mean)) +
  geom_line(data = df.blue.avg,
            aes(x = hpca, y = int_mean, group = id),
            linetype = 'dashed') +
  stat_pvalue_manual(df.blue.stat, label = 'p.signif',
                     hide.ns = FALSE, size = 3, family='arial') +
  scale_fill_manual(name = NULL,
                    values = c('black', 'red')) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.05, 0.175),
                     breaks = seq(-0.5, 0.5, 0.025)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        text = element_text(family='arial', face="bold", size=9),
        panel.border = element_rect(color = "blue", 
                                    fill = NA, 
                                    size = 2))
ggsave('box_blue.png', box.blue, 
       width = 3, height = 6, units = 'cm', dpi = 300)

##### AREA BOXPLOT #####
df.area <- data.frame(id = c('00', '02', '03', '04', '05'),
                      hpca = c('N75K', 'N75K','N75K','N75K','N75K'),
                      wt = c(4903, 1087, 1133, 1075, 538),
                      mut = c(2909, 202, 199, 115, 64),
                      rel = c(0.5933, 0.1858, 0.1756, 0.107, 0.119))

df.area.stat <- df.area %>%
  wilcox_test(rel ~ 1, mu = 1, alternative = 'less') %>%
  add_significance()
  
box.area <- ggplot() +
  geom_boxplot(data = df.area,
               aes(x = hpca, y = rel),
               fill = 'red', alpha = .6) +
  geom_point(data = df.area,
             aes(x = hpca, y = rel)) +
  geom_text(data = df.area,
            aes(x = hpca, y = 0.75), label = '*', size = 5) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = 'N75K area / WT area',
                     limits = c(0, 1),
                     breaks = seq(-2, 2, 0.25)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        text = element_text(family='arial', face="bold", size=9))

ggsave('box_area.png', box.area, 
       width = 3, height = 6, units = 'cm', dpi = 300)
