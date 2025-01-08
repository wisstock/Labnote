# NMDA ionophoresis, ionophoresis cloud data analysis for bio-protocol
# Copyright © 2025 Borys Olifirov

require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)
require(broom)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol')

df.sweep.dF <- read.csv('cloud_sweeps_dF.csv') %>%
  select(-X) %>%
  mutate(id = as.factor(id), i_app = as.factor(i_app), roi = as.factor(roi)) %>%
  mutate(i_app = factor(i_app, c('25', '50', '75', '100'), ordered = TRUE))

df.sweep.abs <- read.csv('cloud_sweeps_abs.csv') %>%
  select(-X) %>%
  mutate(id = as.factor(id), i_app = as.factor(i_app), roi = as.factor(roi)) %>%
  mutate(i_app = factor(i_app, c('25', '50', '75', '100'), ordered = TRUE)) %>%
  group_by(i_app, id) %>%
  mutate(roi_type = as.factor(case_when(roi %in% c(5) ~ 'Max',
                                        roi %in% c(2,4) ~ 'Mid',
                                        roi %in% c(1,3) ~ 'Min',)),
         roi_position = as.factor(case_when(roi %in% c(3,4,5) ~ 'Up',
                                            roi %in% c(1,2) ~ 'Down')),
         roi_position = factor(roi_position, c('Up', 'Down'), ordered = TRUE),
         roi_name = case_when(roi == 1 ~ 'Min down',
                              roi == 2 ~ 'Mid down',
                              roi == 3 ~ 'Min up',
                              roi == 4 ~ 'Mid up',
                              roi == 5 ~ 'Max'),
         roi_name = factor(roi_name, c('Min up', 'Mid up', 'Max', 'Mid down', 'Min down'),
                           ordered = TRUE))


font.size <- 15
font.fam <- 'Arial'
box.alpha <- 0.6


###### PROF PLOT #####
# ROI
ggplot() +
  annotate('rect', xmin = 0.5, xmax = 17, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black') +
  stat_summary(data = df.sweep.abs %>% filter(i_app == '25'),
               aes(x = index_app-3, y = int,
                   color = roi_type, group = roi, linetype = roi_position),
               fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(data = df.sweep.abs %>% filter(i_app == '25'),
               aes(x = index_app-3, y = int,
                   color = roi_type, group = roi, shape = roi_position),
               fun = median,
               geom = 'point', size = 2) +
  stat_summary(data = df.sweep.abs %>% filter(i_app == '25'),
               aes(x = index_app-3, y = int,
                   color = roi_type, fill = roi_type, group = roi),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  scale_fill_manual(name = "ROI type",
                    values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  scale_color_manual(name = "ROI type",
                       values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  scale_shape_discrete(name = "ROI position") +
  scale_linetype_discrete(name = "ROI position") +
  theme_classic() +
  theme(legend.position = c(0.9, 0.7),
        text=element_text(size = font.size, family = font.fam)) +
  scale_x_continuous(breaks = seq(-5, 60, 5)) +
  labs(caption = 'n = 1/4 (cultures/cells)',
       x = 'Time, s',
       y = 'a.u.')

# I app
ggplot() +
  annotate('rect', xmin = 0.5, xmax = 17, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black') +
  stat_summary(data = df.sweep.abs %>% filter(roi_type == 'Max'),
               aes(x = index_app-3, y = int,
                   color = i_app, group = i_app),
               fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(data = df.sweep.abs %>% filter(roi_type == 'Max'),
               aes(x = index_app-3, y = int,
                   color = i_app, group = i_app, shape = i_app),
               fun = median,
               geom = 'point', size = 2) +
  stat_summary(data = df.sweep.abs %>% filter(roi_type == 'Max'),
               aes(x = index_app-3, y = int,
                   color = i_app, fill = i_app, group = i_app),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  scale_fill_manual(name = "I, nA",
                    values = c('25' = 'grey60', '50' = 'grey50', '75' = 'grey35', '100' = 'grey25')) +
  scale_color_manual(name = "I, nA",
                     values = c('25' = 'grey60', '50' = 'grey50', '75' = 'grey35', '100' = 'grey25')) +
  scale_shape_discrete(name = "I, nA") +
  # scale_linetype_discrete(name = "ROI position") +
  theme_classic() +
  theme(legend.position = c(0.9, 0.7),
        text=element_text(size = font.size, family = font.fam)) +
  scale_x_continuous(breaks = seq(-5, 60, 5)) +
  labs(caption = 'n = 1/4 (cultures/cells)',
       x = 'Time, s',
       y = 'a.u.')

##### SPATIAL STAT #####
df.abs.roi.box <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app == 15, i_app == '100') %>%
  ungroup() %>%
  select(int, roi_name, id, roi_type)

df.abs.roi.stat <- df.abs.roi.box %>%
  pairwise_wilcox_test(int ~ roi_name, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position()

ggplot(data = df.abs.roi.box,
       aes(x = roi_name, y = int)) +
  geom_boxplot(aes(fill = roi_type), alpha = .7) +
  stat_pvalue_manual(df.abs.roi.stat, label = 'p.adj.signif',
                     size = font.size - 10) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam)) +
  labs(caption = 'n = 1/4 (cultures/cells)',
       x = 'ROI',
       y = 'a.u.')


##### I STAT #####
df.abs.i.box <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app == 15, roi_position == 'Up') 

df.abs.i.stat <- df.abs.i.box %>%
  ungroup() %>%
  select(i_app, int, roi_type) %>%
  group_by(roi_type) %>%
  kruskal_test(int ~ i_app) %>%
  mutate(title = paste('H Test p=', p, sep=''),
         int = 3000,
         i_app = factor('50', c('25', '50', '75', '100'), ordered = TRUE)) %>%
  ungroup()

ggplot(data = df.abs.i.box, aes(x = i_app, y = int, fill = roi_type)) +
  geom_boxplot(alpha = .7) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  geom_text(data = df.abs.i.stat,
            aes(label = title)) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam)) +
  facet_wrap(~roi_type, ncol = 3) +
  labs(caption = 'n = 1/4 (cultures/cells)',
       x = 'I, nA',
       y = 'a.u.')


##### ROI TAU STAT #####
df.abs.decay.roi <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app %in% seq(17,50), i_app == '25') %>%
  mutate(t = index_app - 16) %>%
  ungroup() %>%
  select(-roi_position, -i_app, -roi, -index_app) %>%
  nest_by(roi_name, id) %>%
  mutate(fit = list(nls(int ~ SSasymp(t,a,a0,log_tau), data = data)))  %>%
  reframe(tidy(fit)) %>%
  filter(term == 'log_tau') %>%
  select(-term) %>%
  mutate(estimate = exp(estimate),
         tau = 1/estimate) %>%
  select(roi_name, id, tau) %>%
  mutate(roi_type = as.factor(c('Min','Min','Min','Min',
                           'Mid','Mid','Mid','Mid',
                           'Max','Max','Max','Max',
                           'Mid','Mid','Mid','Mid',
                           'Min','Min','Min','Min')))

ggplot(data = df.abs.decay.roi, aes(x = roi_name, y = tau, fill = roi_type)) +
  geom_boxplot() +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  # geom_text(data = df.abs.i.stat,
  #           aes(label = title)) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam)) +
  labs(caption = 'n = 1/4 (cultures/cells)',
       x = 'ROI',
       y = 'Tau, s')


df.abs.rise.roi <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app %in% seq(0,17), i_app == '100') %>%
  mutate(t = index_app) %>%
  ungroup() %>%
  select(-roi_type, -roi_position, -i_app, -roi, -index_app) %>%
  nest_by(roi_name, id) %>%
  mutate(fit = list(nls(int ~ SSasymp(t,a,a0,log_tau), data = data)))  %>%
  reframe(tidy(fit)) %>%
  filter(term == 'log_tau') %>%
  select(-term) %>%
  mutate(estimate = exp(estimate),
         tau = 1/estimate) %>%
  select(roi_name, id, tau)

ggplot(data = df.abs.rise.roi, aes(x = roi_name, y = tau, fill = roi_name)) +
  geom_boxplot()


##### I TAU STAT #####
# https://douglas-watson.github.io/post/2018-09_dplyr_curve_fitting/
df.abs.decay <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app %in% seq(17,50), roi_name == 'Max') %>%
  mutate(t = index_app - 16) %>%
  ungroup() %>%
  select(-roi_type, -roi_position, -roi_name, -roi, -index_app) %>%
  nest_by(i_app, id) %>%
  mutate(fit = list(nls(int ~ SSasymp(t,a,a0,log_tau), data = data)))  %>%
  reframe(tidy(fit)) %>%
  filter(term == 'log_tau') %>%
  select(-term) %>%
  mutate(estimate = exp(estimate),
         tau = 1/estimate) %>%
  select(i_app, id, tau)

ggplot(data = df.abs.decay, aes(x = i_app, y = tau, fill = i_app)) +
  geom_boxplot()


df.abs.rise <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app %in% seq(0,17), roi_name == 'Max') %>%
  mutate(t = index_app) %>%
  ungroup() %>%
  select(-roi_type, -roi_position, -roi_name, -roi, -index_app) %>%
  nest_by(i_app, id) %>%
  mutate(fit = list(nls(int ~ SSasymp(t,a,a0,log_tau), data = data)))  %>%
  reframe(tidy(fit)) %>%
  filter(term == 'log_tau') %>%
  select(-term) %>%
  mutate(estimate = exp(estimate),
         tau = 1/estimate) %>%
  select(i_app, id, tau)

ggplot(data = df.abs.rise, aes(x = i_app, y = tau, fill = i_app)) +
  geom_boxplot()
