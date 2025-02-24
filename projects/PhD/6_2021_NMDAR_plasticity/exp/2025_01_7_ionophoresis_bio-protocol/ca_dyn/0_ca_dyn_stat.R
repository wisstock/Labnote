# NMDA ionophoresis, ionophoresis cloud data analysis for bio-protocol
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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/ca_dyn')

df <- read.csv('0_df_ca_dyn.csv') %>%
  mutate(id = as.factor(id),
         roi = as.factor(roi),
         roi_type = as.factor(roi_type),
         index = index - 5) %>%
  select(-X)

font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6


##### AMP STAT #####
df.max.amp <- df %>%
  filter(index == 20)
  group_by(id, roi_type) %>%
  mutate(int_mean = mean(int)) %>%
  select(id, roi_type, int_mean) %>%
  distinct() %>%
  ungroup()

df.max.amp %>%
  group_by(roi_type) %>%
  summarise(median =  median(int), iqr = IQR(int))


df.amp.stat <- df.max.amp %>%
  filter(int < 2.6) %>%
  kruskal_test(int ~ roi_type) %>%
  add_significance() %>%
  mutate(title = paste('H Test p=', p, sep=''),
         roi_type = factor('Max', c('Max', 'Mid', 'Min')),
         int_mean = 2.0)

df.amp.roi.comp <- df %>%
  filter(index == 20, int < 2.6) %>%
  pairwise_wilcox_test(int ~ roi_type, p.adjust.method = 'BH') %>%
  add_significance()


boxplot_max_amp <- ggplot(data = df %>% filter(index == 20, int < 2.6),
                          aes(x = roi_type, y = int)) +
  geom_boxplot(aes(fill = roi_type), alpha = box.alpha) +
  geom_point(color = 'grey35', size = 1.5) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75, linetype = 'dashed', color = 'grey25') +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5, color = 'grey25') +
  geom_text(data = df.amp.roi.stat, aes(label = title)) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        # axis.text.x = element_text(angle = 90, vjust = 0.7),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank()
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 2/4/34 (cultures/cells/ROIs)',
       x = 'ROI type',
       y = expression(ΔF/F[0]))

boxplot_max_amp
save_plot('0_pic_boxplot_max_amp.png', boxplot_max_amp, base_width = 3.25, base_height = 4, dpi = 300)
remove(boxplot_max_amp)


##### RISE TAU #####
df.rise <- df %>%
  filter(index > 0 & index < 20)

ggplot(data = df.rise, aes(x = index, y = int_median, color = roi_type)) +
  geom_line() +
  facet_wrap(~id, ncol = 4)


df.rise.fit <- df.rise %>%
  select(id, roi_type, int, index, roi) %>%
  nest_by(roi_type, id, roi) %>%
  mutate(fit = list(nlsLM(int ~ A - R0 * exp(-B * index),
                          start = list(A = 1.5, R0 = 1, B = 5),
                          data = data)))  %>%
  reframe(tidy(fit)) %>%
  pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic, p.value)) %>%
  mutate(tau = 1/estimate_B)

df.rise.fit %>%
  group_by(roi_type) %>%
  summarise(median =  median(tau), iqr = IQR(tau))


df.rise.tau.stat <- df.rise.fit %>%
  kruskal_test(tau ~ roi_type) %>%
  mutate(title = paste('H Test p=', p, sep=''),
         roi_type = factor('Mid', c('Max', 'Mid', 'Min')),
         tau = 30)

boxplot_rise_tau <- ggplot(data = df.rise.fit, aes(x = roi_type, y = tau)) +
  geom_boxplot(aes(fill = roi_type), alpha = box.alpha) +
  geom_point(color = 'grey35', size = 1.5) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75, linetype = 'dashed', color = 'grey25') +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5, color = 'grey25') +
  geom_text(data = df.rise.tau.stat, aes(label = title)) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  scale_y_continuous(limits = c(0,30), breaks = seq(0,30,5)) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 2/4/35 (cultures/cells/ROIs)',
       x = 'ROI type',
       y = 'Rise \u2CA7, s')

boxplot_rise_tau
save_plot('0_pic_boxplot_rise_tau.png', boxplot_rise_tau, base_width = 3.25, base_height = 4, dpi = 300)
remove(boxplot_rise_tau)


##### REPRESENTATIVE PROFILES #####
# "24_05_08_cell7_ch0"    "24_05_08_cell7_ch0_up" "24_05_9_cell01_ch0"    "24_05_9_cell02_ch0"  
prof_rep <- ggplot(data = df %>% filter(id == '24_05_9_cell02_ch0'),
       aes(x = index, y = int, color = roi_type, fill = roi_type, group = roi_type)) +
  annotate('rect', xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(fun = median,
               geom = 'point', size = 2) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
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
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4)) +
  scale_x_continuous(breaks = c(-5, 0, 10, 20, 30, 40, 50)) +
  labs(x = 'Time, s',
       caption = 'n = 1/1/9 (cultures/cells/ROIs)',
       y = expression(ΔF/F[0]))

prof_rep
save_plot('0_pic_prof_rep.png', prof_rep, base_width = 7, base_height = 4, dpi = 300)
remove(prof_rep)

### all cells
ggplot(data = df,
       aes(x = index, y = int, color = roi_type, fill = roi_type, group = roi_type)) +
  annotate('rect', xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'black') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 0.75) +
  stat_summary(fun = median,
               geom = 'point', size = 2) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
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
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4)) +
  scale_x_continuous(breaks = c(-5, 0, 10, 20, 30, 40, 50)) +
  labs(x = 'Time, s',
      y = expression(ΔF/F[0])) +
  facet_wrap(~id, nrow=4, scales = "free_y", strip.position = 'right')
