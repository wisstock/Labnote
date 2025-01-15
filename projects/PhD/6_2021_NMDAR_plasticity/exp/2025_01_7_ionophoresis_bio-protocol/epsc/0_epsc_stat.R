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

font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6

base_70 <- c(0,150)
base_40 <- c(150,170)
iono_40 <- c(170, 230)
post_40 <- c(230, 250)
post_70 <- c(250, 400)
ends_70 <- c(400, 550)

df.15 <- read.csv('24_10_14_4_60s_15mm.csv') %>%
  select(-X) %>%
  mutate_at(colnames(df)[1:10], as.numeric) %>%
  mutate(id = as.factor('24_10_14_4'),
         hold_interval = as.factor(case_when(between(TimRef, base_70[1], base_70[2]) ~ 'base_70',
                                             between(TimRef, base_40[1], base_40[2]) ~ 'base_40',
                                             between(TimRef, iono_40[1], iono_40[2]) ~ 'iono_40',
                                             between(TimRef, post_40[1], post_40[2]) ~ 'post_40',
                                             between(TimRef, post_70[1], post_70[2]) ~ 'post_70',
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
         Amp_norm = Amp / median(Amp[time_interval == '0'])) %>%
  ungroup() %>%
  filter(Amp < 150)

df <- read.csv('epsc_olf.csv') %>%
  select(-X) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(time_interval = as.factor(time_interval)) %>%
  mutate(hold_interval = factor(hold_interval, c('base_70', 'base_40', 'iono_40', 'post_40', 'post_70', 'ends_70'), ordered = TRUE)) %>%
  filter(Amp < 75)




##### HIST #####
ggplot(data = df, aes(x = Amp_base, color = hold_interval, group = hold_interval)) +
  stat_ecdf() +
  facet_wrap(~id, nrow = nlevels(df$id))

ggplot(data = df, aes(x = Amp, fill = id, group = id, color = id)) +
  geom_histogram(alpha = .75, bins = 100) +
  facet_wrap(~hold_interval, ncol = 5)


###### STAT OVERVIEW #####
df.hold.stat <- df %>%
  group_by(id) %>%
  pairwise_wilcox_test(Amp_base ~ hold_interval,
                       p.adjust.method = 'BH', detailed = TRUE) %>%
  add_significance() %>%
  add_xy_position()

ggplot(data = df, aes(x = hold_interval, y = Amp_base,
                      fill = id)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_boxplot(aes(group = hold_interval)) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5) +
  stat_pvalue_manual(df.hold.stat, label = 'p.adj.signif', hide.ns = TRUE,
                     tip.length = 0.01) +
  # stat_summary(fun.min = function(z) { quantile(z,0.25) },
  #              fun.max = function(z) { quantile(z,0.75) },
  #              fun = median,
  #              geom = 'errorbar', width = .3) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.7)) +
  facet_wrap(~id, ncol = nlevels(df$id))

##### LTD STAT #####
df.ltd <- df %>%
  filter(id %in% c('24_10_17_3',
                   '24_10_22_6',
                   '24_10_23_2',
                   '24_10_23_2',
                   '24_10_24_19'))

ggplot(data = df.ltd, aes(x = hold_interval, y = Amp_base,
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



###### CDF STAT #####
df.ltd.cdf <- df.ltd %>%
  filter(hold_interval %in% c('base_70', 'post_70', 'ends_70'))

ks.test(df.ltd.cdf$Amp[df.ltd.cdf$hold_interval == 'base_70'],
        df.ltd.cdf$Amp[df.ltd.cdf$hold_interval == 'post_70'])

ks.test(df.ltd.cdf$Amp[df.ltd.cdf$hold_interval == 'base_70'],
        df.ltd.cdf$Amp[df.ltd.cdf$hold_interval == 'ends_70'])

ggplot(data = df.ltd.cdf,
       aes(x = Amp, color = hold_interval, group = hold_interval)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  stat_ecdf(geom = 'point', alpha = .75) +
  annotate("text", label = 'KS-test t0~t2 p << 0.001', x = 7, y = 0.97) +
  annotate("text", label = 'KS-test t0~t5 p << 0.001', x = 7, y = 0.9) +
  theme_classic() +
  theme(text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.position = c(0.9, 0.25)) +
  scale_color_manual(values = c('base_70' = 'black', 'post_70' = 'coral', 'ends_70' = 'deepskyblue'),
                     labels = c('base_70' = 't0', 'post_70' = 't2', 'ends_70' = 't5')) +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  guides(color=guide_legend('Time interval')) +
  labs(caption = 'n = 3/4 (cultures/cells)',
       x = 'mEPSC amplitude, pA',
       y = 'Cumulative probability')

##### AGGREGATED STAT #####

df.ltd.med <- df.ltd %>%
  filter(hold_interval %in% c('base_70', 'post_70', 'ends_70')) %>%
  group_by(hold_interval, id) %>%
  mutate(Amp_base_med = median(Amp_base)) %>%
  select(id, hold_interval, Amp_base_med) %>%
  distinct() %>%
  ungroup()

ggplot(data = df.ltd.med, aes(x = hold_interval, y = Amp_base_med,
                              color = id)) +
  geom_line(aes(group = id), size = 0.2) +
  geom_point(size = 0.2) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  stat_summary(aes(group = -1),
               fun = median,
               geom = 'line', size = 0.75, color = 'black') +
  stat_summary(aes(group = -1),
               fun = median,
               geom = 'point', size = 1.5, color = 'black') +
  stat_summary(aes(group = -1),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .3, color = 'black') +
  scale_x_discrete(labels=c("base_70" = "t0",
                            "base_40" = "B -40 mV",
                            "iono_40" = "I -40 mV",
                            "post_40" = "WO -40 mV",
                            "post_70" = "t2",
                            "ends_70" = "t5")) +
  theme_classic() +
  theme(legend.position = 'none')
