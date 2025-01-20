# NMDA ionophoresis, EPSCs data analysis for bio-protocol
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
require(wesanderson)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/epsc')

font.size <- 19
font.fam <- 'Arial'
box.alpha <- 0.6

nmda.100.color <- 'magenta' 
nmda.15.color <- 'green3'

base_70 <- c(0,150)
base_40 <- c(150,170)
iono_40 <- c(170, 230)
post_40 <- c(230, 250)
post_70 <- c(250, 400)
ends_70 <- c(400, 550)

df <- read.csv('epsc_olf.csv') %>%
  select(-X) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(time_interval = as.factor(time_interval)) %>%
  mutate(hold_interval = factor(hold_interval, c('base_70', 'base_40', 'iono_40', 'post_40', 'post_70', 'ends_70'), ordered = TRUE),
         mid_interval = factor(mid_interval, c('t-75', 't0', '-40', 't75', 't150', 't225', 't300'), ordered = TRUE)) %>%
  filter(Amp < 75)

df.ltd <- df %>%
  filter(id %in% c('24_10_17_3',
                   '24_10_22_6',
                   '24_10_23_2',
                   '24_10_16_11',
                   '24_10_24_19'))

df.no.ltd <- df %>%
  filter(id %in% c('24_10_16_11',
                   '24_10_22_10',
                   '24_10_22_8',
                   '24_10_24_7'))


df.15 <- read.csv('epsc_olf_15mm.csv') %>%
  select(-X) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(hold_interval = factor(hold_interval, c('base_70', 'base_40', 'iono_40', 'post_40', 'post_70', 'ends_70'), ordered = TRUE),
         mid_interval = factor(mid_interval, c('t-75', 't0', '-40', 't75', 't150', 't225', 't300'), ordered = TRUE)) %>%
  filter(Amp < 75)
  #        id %in% c('24_10_14_4',
  #                  '25_01_7_1',
  #                  '25_01_8_9',
  #                  '25_01_9_2'))

##### HIST #####
# ggplot(data = df, aes(x = Amp_base, color = hold_interval, group = hold_interval)) +
#   stat_ecdf() +
#   facet_wrap(~id, nrow = nlevels(df$id))
# 
# ggplot(data = df, aes(x = Amp, fill = id, group = id, color = id)) +
#   geom_histogram(alpha = .75, bins = 100) +
#   facet_wrap(~hold_interval, ncol = 5)


###### STAT OVERVIEW #####
## 100 mM
df.hold.stat <- df %>%
  group_by(id) %>%
  pairwise_wilcox_test(Amp_mid ~ mid_interval,
                       p.adjust.method = 'BH', detailed = TRUE) %>%
  add_significance() %>%
  add_xy_position()

ggplot(data = df %>% filter(mid_interval != '-40'),
       aes(x = mid_interval, y = Amp_mid, fill = id)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_boxplot(aes(group = mid_interval)) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5) +
  stat_pvalue_manual(df.hold.stat, label = 'p.adj.signif', hide.ns = TRUE,
                     tip.length = 0.01) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.7)) +
  labs(title = '100 mM') +
  facet_wrap(~id, ncol = nlevels(df$id))

## 15 mM
df.15.hold.stat <- df.15 %>%
  filter(mid_interval != '-40') %>%
  group_by(id) %>%
  pairwise_wilcox_test(Amp_mid ~ mid_interval,
                       p.adjust.method = 'BH', detailed = TRUE) %>%
  add_significance() %>%
  add_xy_position()

ggplot(data = df.15 %>% filter(mid_interval != '-40'),
       aes(x = mid_interval, y = Amp_mid, fill = id)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_boxplot(aes(group = mid_interval)) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5) +
  stat_pvalue_manual(df.15.hold.stat, label = 'p.adj.signif', hide.ns = TRUE,
                     tip.length = 0.01) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.7)) +
  labs(title = '15 mM') +
  facet_wrap(~id, ncol = nlevels(df.15$id))


remove(df.15.hold.stat, df.hold.stat)

## LTD STAT
ggplot(data = df.ltd, aes(x = mid_interval, y = Amp_base,
                          color = id)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .3) +
  scale_x_discrete(labels=c("base_70" = "B -70 mV",
                            "base_40" = "B -40 mV",
                            "iono_40" = "I -40 mV",
                            "post_40" = "WO -40 mV",
                            "post_70" = "WO -70 mV 250-400s",
                            "ends_70" = "WO -70 mV 400-550s")) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, vjust = 0.5))




###### AMP CDF STAT #####
df.ltd.cdf <- df.ltd %>%
  filter(mid_interval %in% c('t-75', 't0', 't150', 't300')) %>%
  droplevels()

df.ltd.cdf %>%
  select(-id) %>%
  group_by(mid_interval) %>%
  summarise(med = median(Amp), iqr = IQR(Amp))

df.ltd.cdf %>%
  select(-id) %>%
  group_by(mid_interval) %>%
  summarise(med = median(IEI), iqr = IQR(IEI))


df.15.cdf <- df.15 %>%
  filter(mid_interval %in% c('t-75', 't0', 't150', 't300')) %>%
  droplevels()

df.15.cdf %>%
  select(-id) %>%
  group_by(mid_interval) %>%
  summarise(med = median(Amp), iqr = IQR(Amp))

df.15.cdf %>%
  select(-id) %>%
  group_by(mid_interval) %>%
  summarise(med = median(IEI), iqr = IQR(IEI))


set.seed(1)
ks.test(sample(df.ltd.cdf$Amp[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$Amp[df.15.cdf$mid_interval == 't-75'], size = samp.size))

## AMP 100
samp.size <- 500

set.seed(1)
ks.test(sample(df.ltd.cdf$Amp[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.ltd.cdf$Amp[df.ltd.cdf$mid_interval == 't0'], size = samp.size))

set.seed(1)
ks.test(sample(df.ltd.cdf$Amp[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.ltd.cdf$Amp[df.ltd.cdf$mid_interval == 't150'], size = samp.size))

set.seed(1)
ks.test(sample(df.ltd.cdf$Amp[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.ltd.cdf$Amp[df.ltd.cdf$mid_interval == 't300'], size = samp.size))


cdf_amp_100 <- ggplot(data = df.ltd.cdf,
       aes(x = Amp, color = mid_interval, group = mid_interval)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  stat_ecdf(geom = 'point', alpha = .75) +
  annotate("text", label = 'KS-test t-75~t0 p=0.129', hjust = 0,
           size = font.size - 14, x = 24, y = 0.11) +
  annotate("text", label = 'KS-test t-75~t150 p<0.001', hjust = 0,
           size = font.size - 14, x = 24, y = 0.06) +
  annotate("text", label = 'KS-test t-75~t300 p<0.001', hjust = 0,
           size = font.size - 14, x = 24, y = 0.01) +
  # annotate("text", label = 'NMDA 100 mM', hjust = 0.5, fontface = 2,
  #          size = font.size - 12, x = 0, y = 1, color = nmda.100.color) +
  guides(color=guide_legend('Time interval')) +
  labs(subtitle = 'NMDA 100 mM',
       caption = 'n = 5/5 (cultures/cells)',
       x = 'Amplitude, pA',
       y = 'Cumulative probability') +
  theme_classic() +
  theme(text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        plot.subtitle = element_text(hjust = 0.5, size = font.size - 5,
                                     face = 'bold', family = font.fam, color = nmda.100.color),
        legend.position = c(0.78, 0.69),
        legend.title = element_text(size = font.size-3),
        panel.border = element_rect(color = nmda.100.color, fill = NA,  size = 3)) +
  scale_color_manual(values = c('t-75' = 'black', 't0' = 'coral',
                                't150' = 'deepskyblue', 't300' = 'deeppink2')) +
  scale_x_continuous(breaks = seq(0, 100, 10),
                     limits = c(0,75))

cdf_amp_100
save_plot('0_pic_cdf_amp_100.png', cdf_amp_100, base_width = 4.5, base_height = 5, dpi = 300)
remove(cdf_amp_100)


## AMP 15
samp.size <- 500

set.seed(1)
ks.test(sample(df.15.cdf$Amp[df.15.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$Amp[df.15.cdf$mid_interval == 't0'], size = samp.size))

set.seed(1)
ks.test(sample(df.15.cdf$Amp[df.15.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$Amp[df.15.cdf$mid_interval == 't150'], size = samp.size))

set.seed(1)
ks.test(sample(df.15.cdf$Amp[df.15.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$Amp[df.15.cdf$mid_interval == 't300'], size = samp.size))


cdf_amp_15 <- ggplot(data = df.15.cdf,
                     aes(x = Amp, color = mid_interval, group = mid_interval)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  stat_ecdf(geom = 'point', alpha = .75) +
  annotate("text", label = 'KS-test t-75~t0 p<0.001', hjust = 0,
           size = font.size - 14, x = 24, y = 0.11) +
  annotate("text", label = 'KS-test t-75~t150 p<0.05', hjust = 0,
           size = font.size - 14, x = 24, y = 0.06) +
  annotate("text", label = 'KS-test t-75~t300 p<0.001', hjust = 0,
           size = font.size - 14, x = 24, y = 0.01) +
  guides(color=guide_legend('Time interval')) +
  labs(subtitle = 'NMDA 15 mM',
       caption = 'n = 5/6 (cultures/cells)',
       x = 'Amplitude, pA',
       y = 'Cumulative probability') +
  theme_classic() +
  theme(text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        plot.subtitle = element_text(hjust = 0.5, size = font.size - 5,
                                     face = 'bold', family = font.fam, color = nmda.15.color),
        legend.position = 'none',
        legend.title = element_text(size = font.size-2),
        panel.border = element_rect(color = nmda.15.color, fill = NA,  size = 3)) +
  scale_color_manual(values = c('t-75' = 'black', 't0' = 'coral',
                                't150' = 'deepskyblue', 't300' = 'deeppink2')) +
  scale_x_continuous(breaks = seq(0, 100, 10),
                     limits = c(0,75))

cdf_amp_15
save_plot('0_pic_cdf_amp_15.png', cdf_amp_15, base_width = 4.5, base_height = 5, dpi = 300)
remove(cdf_amp_15)



##### IEI CDF STAT #####

set.seed(1)
ks.test(sample(df.ltd.cdf$IEI[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$IEI[df.15.cdf$mid_interval == 't-75'], size = samp.size))

## IEI 100
samp.size <- 500

set.seed(1)
ks.test(sample(df.ltd.cdf$IEI[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.ltd.cdf$IEI[df.ltd.cdf$mid_interval == 't0'], size = samp.size))

set.seed(1)
ks.test(sample(df.ltd.cdf$IEI[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.ltd.cdf$IEI[df.ltd.cdf$mid_interval == 't150'], size = samp.size))

set.seed(1)
ks.test(sample(df.ltd.cdf$IEI[df.ltd.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.ltd.cdf$IEI[df.ltd.cdf$mid_interval == 't300'], size = samp.size))


cdf_iei_100 <- ggplot(data = df.ltd.cdf,
       aes(x = IEI, color = mid_interval, group = mid_interval)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  stat_ecdf(geom = 'point', alpha = .75) +
  annotate("text", label = 'KS-test t-75~t0 p=0.863', hjust = 0,
           size = font.size - 14, x = 0.31, y = 0.11) +
  annotate("text", label = 'KS-test t-75~t150 p=0.095', hjust = 0,
           size = font.size - 14, x = 0.31, y = 0.06) +
  annotate("text", label = 'KS-test t-75~t300 p=0.082', hjust = 0,
           size = font.size - 14, x = 0.31, y = 0.01) +
  guides(color=guide_legend('Time interval')) +
  labs(subtitle = 'NMDA 100 mM',
       caption = 'n = 5/5 (cultures/cells)',
       x = 'IEI, s',
       y = 'Cumulative probability') +
  theme_classic() +
  theme(text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        plot.subtitle = element_text(hjust = 0.5, size = font.size - 5,
                                     face = 'bold', family = font.fam, color = nmda.100.color),
        legend.position = 'none',
        legend.title = element_text(size = font.size-2),
        panel.border = element_rect(color = nmda.100.color, fill = NA,  size = 3)) +
  scale_color_manual(values = c('t-75' = 'black', 't0' = 'coral',
                                't150' = 'deepskyblue', 't300' = 'deeppink2')) +
  scale_x_continuous(breaks = seq(0, 100, 0.25),
                     limits = c(0,1))

cdf_iei_100
save_plot('0_pic_cdf_iei_100.png', cdf_iei_100, base_width = 4.5, base_height = 5, dpi = 300)
remove(cdf_iei_100)


## IEI 15
samp.size <- 500

set.seed(1)
ks.test(sample(df.15.cdf$IEI[df.15.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$IEI[df.15.cdf$mid_interval == 't0'], size = samp.size))

set.seed(1)
ks.test(sample(df.15.cdf$IEI[df.15.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$IEI[df.15.cdf$mid_interval == 't150'], size = samp.size))

set.seed(1)
ks.test(sample(df.15.cdf$IEI[df.15.cdf$mid_interval == 't-75'], size = samp.size),
        sample(df.15.cdf$IEI[df.15.cdf$mid_interval == 't300'], size = samp.size))


cdf_iei_15 <- ggplot(data = df.15.cdf,
                      aes(x = IEI, color = mid_interval, group = mid_interval)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  stat_ecdf(geom = 'point', alpha = .75) +
  annotate("text", label = 'KS-test t-75~t0 p<0.001', hjust = 0,
           size = font.size - 14, x = 0.31, y = 0.11) +
  annotate("text", label = 'KS-test t-75~t150 p=0.770', hjust = 0,
           size = font.size - 14, x = 0.31, y = 0.06) +
  annotate("text", label = 'KS-test t-75~t300 p<0.001', hjust = 0,
           size = font.size - 14, x = 0.31, y = 0.01) +
  guides(color=guide_legend('Time interval')) +
  labs(subtitle = 'NMDA 15 mM',
       caption = 'n = 5/6 (cultures/cells)',
       x = 'IEI, s',
       y = 'Cumulative probability') +
  theme_classic() +
  theme(text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        plot.subtitle = element_text(hjust = 0.5, size = font.size - 5,
                                     face = 'bold', family = font.fam, color = nmda.15.color),
        legend.position = 'none',
        legend.title = element_text(size = font.size-2),
        panel.border = element_rect(color = nmda.15.color, fill = NA,  size = 3)) +
  scale_color_manual(values = c('t-75' = 'black', 't0' = 'coral',
                                't150' = 'deepskyblue', 't300' = 'deeppink2')) +
  scale_color_manual(values = c('t-75' = 'black', 't0' = 'coral',
                                't150' = 'deepskyblue', 't300' = 'deeppink2')) +
  scale_x_continuous(breaks = seq(0, 100, 0.25),
                     limits = c(0,1))

cdf_iei_15
save_plot('0_pic_cdf_iei_15.png', cdf_iei_15, base_width = 4.5, base_height = 5, dpi = 300)
remove(cdf_iei_15)


## Area
ks.test(df.ltd.cdf$Area[df.ltd.cdf$hold_interval == 'base_70' & df.ltd.cdf$IEI < 1],
        df.ltd.cdf$Area[df.ltd.cdf$hold_interval == 'post_70' & df.ltd.cdf$IEI < 1])

ks.test(df.ltd.cdf$Area[(df.ltd.cdf$hold_interval == 'base_70') & (df.ltd.cdf$IEI < 1)],
        df.ltd.cdf$Area[(df.ltd.cdf$hold_interval == 'ends_70') & (df.ltd.cdf$IEI < 1)])

ggplot(data = df.ltd.cdf,
       aes(x = Area, color = hold_interval, group = hold_interval)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  stat_ecdf(geom = 'point', alpha = .5) +
  annotate("text", label = 'KS-test t0.0~t2.5 p=0.497', x = 700, y = 0.87) +
  annotate("text", label = 'KS-test t0.0~t5.0 p<<0.001', x = 700, y = 0.8) +
  theme_classic() +
  theme(text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.position = c(0.9, 0.25)) +
  scale_color_manual(values = c('base_70' = 'black', 'post_70' = 'coral', 'ends_70' = 'deepskyblue'),
                     labels = c('base_70' = 't0.0', 'post_70' = 't2.5', 'ends_70' = 't5.0')) +
  scale_x_continuous(breaks = seq(0, 1000, 100)) +
  guides(color=guide_legend('Time interval')) +
  labs(caption = 'n = 3/4 (cultures/cells)',
       x = 'sEPSCs area, fC',
       y = 'Cumulative probability')

##### AGGREGATED STAT #####
## AMP
df.100.med <- df.ltd %>%
  # filter(mid_interval %in% c('t-75', 't0', 't150', 't300')) %>%
  filter(mid_interval != '-40') %>%
  group_by(mid_interval, id) %>%
  mutate(Amp_med = median(Amp_mid)) %>%
  select(id, mid_interval, Amp_med) %>%
  distinct() %>%
  ungroup() %>%
  droplevels()

df.15.med <- df.15 %>%
  filter(mid_interval != '-40', id != '25_01_8_4') %>%
  group_by(mid_interval, id) %>%
  mutate(Amp_med = median(Amp_mid)) %>%
  select(id, mid_interval, Amp_med) %>%
  distinct() %>%
  ungroup() %>%
  droplevels()

kruskal.test(Amp_med ~ mid_interval, data = df.15.med)

df.median <- rbind(df.100.med %>% mutate(conc = '100'),
                   df.15.med %>% mutate(conc = '15')) %>%
  mutate(conc = as.factor(conc))
remove(df.100.med, df.15.med)

df.median %>%
  select(-id) %>%
  group_by(mid_interval, conc) %>%
  summarise(med = median(Amp_med), iqr = IQR(Amp_med))


point.shift <- .075

df.median.stat <- df.median %>%
  group_by(conc) %>%
  mutate(mid_interval = factor(mid_interval, ordered = FALSE)) %>%
  pairwise_wilcox_test(Amp_med ~ mid_interval, p.adjust.method = 'BH',
                       detailed = TRUE, ref.group = 't-75') %>%
  add_significance() %>%
  add_xy_position(fun = 'median') %>%
  mutate(y.position = c(0, 0.84, 0.79, 0.78, 0.75,
                        0,0,0.93,0,0),
         xmax = if_else(conc == '15', xmax + point.shift + 0.17, xmax - point.shift - 0.17))

df.median.plot <- df.median %>%
  mutate(mid_interval = as.numeric(mid_interval)) %>%
  mutate(mid_interval = if_else(conc == '100', mid_interval - .075, mid_interval + .075))

prof_amp <- ggplot(data = df.median.plot, aes(x = mid_interval, y = Amp_med,
                                  color = conc, group = conc)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75) +
  stat_pvalue_manual(df.median.stat, label = 'p.adj.signif', size = font.size - 10,
                     hide.ns = TRUE, remove.bracket = TRUE) +
  annotate("text", size = font.size - 14, label = 'H-test 15 mM p=0.386', hjust = 0, x = 4.5, y = 1.015) +
  scale_color_manual(values = c('100' = nmda.100.color, '15' = nmda.15.color)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c('t-75', 't0', 't75', 't150', 't225', 't300')) +
  theme_classic() +
  theme(legend.position = c(0.12, 0.18),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-2)) +
  labs(linetype = 'NMDA, mM',
       color = 'NMDA, mM',
       caption = '100 mM n = 4/5, 15 mM n = 5/5 (cultures/cells)',
       x = 'Time interval',
       y = 'Norm. amplitude')

prof_amp
save_plot('0_pic_prof_amp.png',prof_amp, base_width = 7, base_height = 5, dpi = 300)
remove(prof_amp)

## IEI
df.100.med.iei <- df.ltd %>%
  filter(mid_interval != '-40') %>%
  group_by(mid_interval, id) %>%
  mutate(IEI_med = median(IEI)) %>%
  select(id, mid_interval, IEI_med) %>%
  distinct() %>%
  ungroup() %>%
  droplevels()

df.15.med.iei <- df.15 %>%
  filter(mid_interval != '-40', id != '25_01_8_4') %>%
  group_by(mid_interval, id) %>%
  mutate(IEI_med = median(IEI)) %>%
  select(id, mid_interval, IEI_med) %>%
  distinct() %>%
  ungroup() %>%
  droplevels()

df.median.iei <- rbind(df.100.med.iei %>% mutate(conc = '100'),
                       df.15.med.iei %>% mutate(conc = '15')) %>%
  mutate(conc = as.factor(conc))
remove(df.100.med.iei, df.15.med.iei)


df.median.iei %>%
  select(-id) %>%
  group_by(mid_interval, conc) %>%
  summarise(med = median(IEI_med), iqr = IQR(IEI_med))


df.median.iei.stat <- df.median.iei %>%
  group_by(conc) %>%
  kruskal_test(IEI_med ~ mid_interval) %>%
  add_significance()

df.median.iei.plot <- df.median.iei %>%
  mutate(mid_interval = as.numeric(mid_interval)) %>%
  mutate(mid_interval = if_else(conc == '100', mid_interval - .075, mid_interval + .075))

prof_iei <- ggplot(data = df.median.iei.plot, aes(x = mid_interval, y = IEI_med,
                                      color = conc, group = conc)) +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75) +
  annotate("text", size = font.size - 14, label = 'H-test 100 mM p=0.697', hjust = 0, x = 1, y = 1.5) +
  annotate("text", size = font.size - 14, label = 'H-test 15 mM p=0.989', hjust = 0, x = 1, y = 1.35) +
  scale_color_manual(values = c('100' = nmda.100.color, '15' = nmda.15.color)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c('t-75', 't0', 't75', 't150', 't225', 't300')) +
  scale_y_continuous(breaks = seq(0, 10, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-2)) +
  labs(color = 'NMDA, mM',
       caption = '100 mM n = 4/5, 15 mM n = 5/5 (cultures/cells)',
       x = 'Time interval',
       y = 'IEI, s')

prof_iei
save_plot('0_pic_prof_iei.png', prof_iei, base_width = 7, base_height = 5, dpi = 300)
remove(prof_iei)
