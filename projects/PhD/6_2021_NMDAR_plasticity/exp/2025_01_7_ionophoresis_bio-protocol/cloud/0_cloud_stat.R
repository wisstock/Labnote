# NMDA ionophoresis, ionophoresis cloud data analysis for bio-protocol
# Copyright © 2025 Borys Olifirov

require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)
# require(broom)
require(minpack.lm)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/cloud')

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

font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6


###### PROF PLOT #####
# ROI
roi_profile <-ggplot() +
  annotate('rect', xmin = 0, xmax = 17, ymin = -Inf, ymax = Inf,
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
  theme(text=element_text(size = font.size, family = font.fam),
        # legend.title = element_text(size = font.size-4),
        # legend.text = element_text(size = font.size-4),
        plot.caption = element_text(size = font.size-4)) +
  scale_x_continuous(breaks = seq(-5, 60, 5)) +
  labs(caption = 'n = 1/4 (experimental days/replications)',
       x = 'Time, s',
       y = 'Intensity, a.u.')

roi_profile
save_plot('0_profile_roi.png', roi_profile, base_width = 12, base_height = 3, dpi = 300)
remove(roi_profile)


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

df.abs.roi.box %>%
  group_by(roi_name) %>%
  summarise(median =  median(int), iqr = IQR(int))

df.abs.roi.stat <- df.abs.roi.box %>%
  pairwise_wilcox_test(int ~ roi_name, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(step.increase = 0.075) %>%
  mutate(y.position = c(2000,3400,3523.828,3716.201,3100,4100.947,3750,3250,3550,2000))

roi_boxplot <- ggplot(data = df.abs.roi.box,
       aes(x = roi_name, y = int)) +
  geom_boxplot(aes(fill = roi_type), alpha = box.alpha) +
  stat_pvalue_manual(df.abs.roi.stat, label = 'p.adj.signif',
                     size = font.size - 10, hide.ns = TRUE, tip.length = 0.01) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        axis.text.x = element_text(angle = 30, vjust = 0.7),
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 1/4 (experimental days/replications)',
       x = 'ROI',
       y = 'Intensity, a.u.')

roi_boxplot
save_plot('0_boxplot_roi.png', roi_boxplot, base_width = 4, base_height = 4, dpi = 300)
remove(roi_boxplot)

##### I STAT #####
df.abs.i.box <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app == 15, roi_position == 'Up') 

df.abs.i.box %>%
  group_by(roi_name) %>%
  summarise(median =  median(int), iqr = IQR(int))

df.abs.i.stat <- df.abs.i.box %>%
  ungroup() %>%
  select(i_app, int, roi_type) %>%
  group_by(roi_type) %>%
  kruskal_test(int ~ i_app) %>%
  mutate(title = paste('H Test p=', p, sep=''),
         int = 3000,
         i_app = factor('50', c('25', '50', '75', '100'), ordered = TRUE)) %>%
  ungroup()

i_boxplot <- ggplot(data = df.abs.i.box,
       aes(x = i_app, y = int, color = roi_type)) +
  stat_summary(fun = median, size = 0.5) +
  stat_summary(aes(group = -1),
               fun = median,
               geom = 'line', size = 1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .3) +
  geom_point(color = 'grey35', size = 1.5) +
  geom_text(data = df.abs.i.stat,
            aes(label = title), color = 'black') +
  facet_wrap(~roi_type, ncol = 3) +
  scale_color_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  scale_x_discrete(labels=c("25" = "-25", "50" = "-50",
                            "75" = "-75", "100" = "-100")) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam)) +
  facet_wrap(~roi_type, ncol = 3) +
  labs(caption = 'n = 1/4 (cultures/cells)',
       x = 'I, nA',
       y = 'Intensity, a.u.')

save_plot('0_boxplot_i.png', i_boxplot, base_width = 6, base_height = 4, dpi = 300)
remove(i_boxplot)


##### ROI TAU STAT #####
### RISE
df.abs.rise.roi <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app %in% seq(0,17), i_app == '100') %>%
  mutate(t = index_app) %>%
  ungroup() %>%
  select(-roi_type, -roi_position, -i_app, -roi, -index_app)

df.abs.rise.fit <- df.abs.rise.roi %>%
  nest_by(roi_name, id) %>%
  mutate(fit = list(nlsLM(int ~ A - R0 * exp(-B * t),
                        start = list(A = 5000, R0 = 1000, B = 5),
                        data = data)))  %>%
  reframe(tidy(fit)) %>%
  pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic, p.value)) %>%
  mutate(tau = 1/estimate_B,
         roi_type = as.factor(c('Min','Min','Min','Min',
                                'Mid','Mid','Mid','Mid',
                                'Max','Max','Max','Max',
                                'Mid','Mid','Mid','Mid',
                                'Min','Min','Min','Min')))

df.abs.rise.fit %>%
  group_by(roi_name) %>%
  summarise(median =  median(tau), iqr = IQR(tau))


df.abs.rise.stat <- df.abs.rise.fit %>%
  pairwise_wilcox_test(tau ~ roi_name, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(step.increase = 0.075) %>%
  mutate(y.position = c(13.43980,13.95513,14.47047,14.98580,9.75,10.5,16.53180,0,9,0))

rise_boxplot <- ggplot(data = df.abs.rise.fit, aes(x = roi_name, y = tau)) +
  geom_boxplot(aes(fill = roi_type), alpha = box.alpha) +
  stat_pvalue_manual(df.abs.rise.stat, label = 'p.adj.signif',
                     size = font.size - 10, hide.ns = TRUE, tip.length = 0.01) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        axis.text.x = element_text(angle = 30, vjust = 0.7),
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 1/4 (experimental days/replications)',
       x = 'ROI',
       y = 'Rise \u2CA7, s')

rise_boxplot
save_plot('0_boxplot_rise.png', rise_boxplot, base_width = 4, base_height = 4, dpi = 300)
remove(rise_boxplot)

# # fit line demo
# demo_roi <- 'Min up'
# demo_id <- '24_05_30_cell2_ch1'
# 
# inverse_exp <- function(x, A, R0, B){
#   A - R0 * exp(-B * x)
# }
# 
# df.predict.demo <- df.abs.rise.roi %>%
#   filter(roi_name == demo_roi, id == demo_id) %>%
#   select(-id, -roi_name) %>%
#   mutate(predict = inverse_exp(x = t,
#                                A = df.abs.rise.fit$estimate_A[df.abs.rise.fit$roi_name == demo_roi & df.abs.rise.fit$id == demo_id],
#                                R0 = df.abs.rise.fit$estimate_R0[df.abs.rise.fit$roi_name == demo_roi & df.abs.rise.fit$id == demo_id],
#                                B = df.abs.rise.fit$estimate_B[df.abs.rise.fit$roi_name == demo_roi & df.abs.rise.fit$id == demo_id]))
# 
# ggplot() +
#   geom_point(data = df.predict.demo,
#              aes(x = t, y = int), color = 'blue') +
#   geom_line(data = df.predict.demo,
#              aes(x = t, y = predict), color = 'red')
#   
# 
# remove(df.predict.demo, demo_roi, demo_id, inverse_exp)

### DECAY
df.abs.decay.roi <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app %in% seq(18,50), i_app == '100') %>%
  mutate(t = index_app - 16) %>%
  ungroup() %>%
  select(-roi_position, -i_app, -roi, -index_app)

df.abs.decay.fit <- df.abs.decay.roi %>%
  nest_by(roi_name, id) %>%
  mutate(fit = list(nls(int ~ SSasymp(t, A, R0, L),
                        data = data)))  %>%
  reframe(tidy(fit)) %>%
  pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic, p.value)) %>%
  mutate(B = exp(estimate_L),
         tau = 1/B,
         roi_type = as.factor(c('Min','Min','Min','Min',
                                'Mid','Mid','Mid','Mid',
                                'Max','Max','Max','Max',
                                'Mid','Mid','Mid','Mid',
                                'Min','Min','Min','Min')))

df.abs.decay.fit %>%
  group_by(roi_name) %>%
  summarise(median =  median(tau), iqr = IQR(tau))


df.abs.decay.stat <- df.abs.decay.fit %>%
  kruskal_test(tau ~ roi_name) %>%
  add_significance() %>%
  mutate(title = paste('H Test p=', p, sep=''),
         roi_name = factor('Mid up', c('Min up', 'Mid up', 'Max', 'Mid down', 'Min down'),
                           ordered = TRUE),
         tau = 10)

decay_boxplot <- ggplot(data = df.abs.decay.fit, aes(x = roi_name, y = tau)) +
  geom_boxplot(aes(fill = roi_type), alpha = box.alpha) +
  geom_text(data = df.abs.decay.stat,
            aes(label = title)) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        axis.text.x = element_text(angle = 30, vjust = 0.7),
        plot.caption = element_text(size = font.size-4)) +
  labs(caption = 'n = 1/4 (experimental days/replications)',
       x = 'ROI',
       y = 'Decay \u2CA7, s')

decay_boxplot
save_plot('0_boxplot_decay.png', decay_boxplot, base_width = 4, base_height = 4, dpi = 300)
remove(decay_boxplot)

# # fit line demo
# demo_roi <- 'Min down'
# demo_id <- '24_05_30_cell2_ch1'
# 
# ssasymp <- function(x, A, R0, L){
#   A + R0 * exp(-exp(L) * x)
# }
# 
# df.predict.demo <- df.abs.decay.roi %>%
#   filter(roi_name == demo_roi, id == demo_id) %>%
#   select(-id, -roi_name) %>%
#   mutate(predict = ssasymp(x = t,
#                            A = df.abs.decay.fit$estimate_A[df.abs.rise.fit$roi_name == demo_roi & df.abs.rise.fit$id == demo_id],
#                            R0 = df.abs.decay.fit$estimate_R0[df.abs.rise.fit$roi_name == demo_roi & df.abs.rise.fit$id == demo_id],
#                            L = df.abs.decay.fit$estimate_L[df.abs.rise.fit$roi_name == demo_roi & df.abs.rise.fit$id == demo_id]))
# 
# ggplot() +
#   geom_point(data = df.predict.demo,
#              aes(x = t, y = int), color = 'blue') +
#   geom_line(data = df.predict.demo,
#             aes(x = t, y = predict), color = 'red')
# 
# 
# remove(df.predict.demo, demo_roi, demo_id, inverse_exp)


##### I TAU STAT #####
### RISE
df.abs.rise.i.roi <- df.sweep.abs %>%
  mutate(index_app = index_app - 3) %>%
  filter(index_app %in% seq(0,17), roi_type == 'Max') %>%
  mutate(t = index_app) %>%
  ungroup() %>%
  select(-roi_type, -roi_position, -roi_name, -roi_type, -roi, -index_app)

df.abs.rise.i.fit <- df.abs.rise.i.roi %>%
  nest_by(i_app, id) %>%
  mutate(fit = list(nlsLM(int ~ A - R0 * exp(-B * t),
                          start = list(A = 3000, R0 = 1500, B = 0.1),
                          data = data)))  %>%
  reframe(tidy(fit)) %>%
  pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic, p.value)) %>%
  mutate(tau = 1/estimate_B)

df.abs.rise.i.fit %>%
  group_by(i_app) %>%
  summarise(median =  median(tau), iqr = IQR(tau))

df.abs.rise.i.stat <- df.abs.rise.i.fit %>%
  kruskal_test(tau ~ i_app) %>%
  add_significance() %>%
  mutate(title = paste('H Test p=', p, sep=''),
         i_app = factor('50', c('25', '50', '75', '100'),
                        ordered = TRUE),
         tau = 10)

i_rise_boxplot <- ggplot(data = df.abs.rise.i.fit,
       aes(x = i_app, y = tau)) +
  stat_summary(fun = median, size = 0.5, color = 'red2') +
  stat_summary(aes(group = -1),
               fun = median,
               geom = 'line', size = 1, color = 'red2') +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', color = 'red2', width = .3) +
  geom_point(color = 'grey35', size = 1.5) +
  geom_text(data = df.abs.rise.i.stat,
            aes(label = title)) +
  scale_x_discrete(labels=c("25" = "-25", "50" = "-50",
                            "75" = "-75", "100" = "-100")) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam)) +
  labs(caption = 'n = 1/4 (cultures/cells)',
       x = 'I, nA',
       y = 'Rise \u2CA7, s')

save_plot('0_boxplot_rise_i.png', i_rise_boxplot, base_width = 2.5, base_height = 4, dpi = 300)
remove(i_rise_boxplot)

# ggplot(data = df.abs.rise.i.fit, aes(x = i_app, y = tau)) +
#   geom_boxplot(aes(fill = i_app), alpha = box.alpha) +
#   geom_text(data = df.abs.rise.i.stat,
#             aes(label = title)) +
#   scale_fill_manual(values = c('25' = 'red', '50' = 'red2', '75' = 'red3', '100' = 'red4')) +
#   theme_classic() +
#   theme(legend.position = 'none',
#         text=element_text(size = font.size, family = font.fam)) +
#   labs(caption = 'n = 1/4 (cultures/cells)',
#        x = 'I, nA',
#        y = 'Rise \u2CA7, s')

# # fit line demo
# i <- '75'
# demo_id <- '24_05_30_cell2_ch1'
# 
# inverse_exp <- function(x, A, R0, B){
#   A - R0 * exp(-B * x)
# }
# 
# df.predict.demo <- df.abs.rise.i.roi %>%
#   filter(i_app == i, id == demo_id) %>%
#   select(-id, -i_app) %>%
#   mutate(predict = inverse_exp(x = t,
#                                A = df.abs.rise.i.fit$estimate_A[df.abs.rise.fit$i_app == i & df.abs.rise.fit$id == demo_id],
#                                R0 = df.abs.rise.i.fit$estimate_R0[df.abs.rise.fit$i_app == i & df.abs.rise.fit$id == demo_id],
#                                B = df.abs.rise.i.fit$estimate_B[df.abs.rise.fit$i_app == i & df.abs.rise.fit$id == demo_id]))
# 
# ggplot() +
#   geom_point(data = df.predict.demo,
#              aes(x = t, y = int), color = 'blue') +
#   geom_line(data = df.predict.demo,
#              aes(x = t, y = predict), color = 'red')
# 
# remove(df.predict.demo, i, demo_id, inverse_exp)
