# GLT count in hippocampus pyramid layer
# Copyright © 2024 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio/note/projects/PhD/9_2024_GLT_in_TBI/exp/2024_02_14_glt_in_pyramid')


##### DATA PREPROCESSING #####
# df.glt.area <- read.csv('GLTcount_csv.csv', sep = ';') %>%
#   pivot_longer(cols = Cont:TBI_14D_cef,
#                names_to = 'group', values_to = 'area') %>%
#   filter(area != '') %>%
#   mutate(group = as.factor(group), area = as.numeric(as.character(gsub(',','.',area)))) %>%
#   mutate(trauma = as.factor(ifelse(grepl('Cont', group),'Cont','TBI'))) %>%
#   mutate(treat = as.factor(ifelse(grepl('cef', group),'Cef','Placebo'))) %>%
#   mutate(day = case_when((group == 'Cont') | (group == 'Cont_cef') ~ 0,
#                          (group == 'TBI_3D') | (group == 'TBI_3D_cef') ~ 3,
#                          (group == 'TBI_7D') | (group == 'TBI_7D_cef') ~ 7,
#                          (group == 'TBI_14D') | (group == 'TBI_14D_cef') ~ 14))
  
df.glt <- read.csv('GLTcount_compar_all.csv', sep = ';') %>%
  mutate(relative_area = as.numeric(as.character(gsub(',','.',relative_area))),
         relative_intensity = as.numeric(as.character(gsub(',','.',relative_intensity)))) %>%
  mutate(trauma = as.factor(ifelse(grepl('Cont', group),'Cont','TBI'))) %>%
  mutate(treat = as.factor(ifelse(grepl('cef', group),'Cef','Placebo'))) %>%
  mutate(day = case_when((group == 'Cont') | (group == 'Cont_cef') ~ 0,
                         (group == 'TBI_3D') | (group == 'TBI_3D_cef') ~ 3,
                         (group == 'TBI_7D') | (group == 'TBI_7D_cef') ~ 7,
                         (group == 'TBI_14D') | (group == 'TBI_14D_cef') ~ 14))


##### STAT OVERVIEW #####
kruskal.test(relative_area ~ group, data = df.glt)
kruskal.test(relative_intensity ~ group, data = df.glt)

# AREA
# pairwise test between all groups
df.area.stat.tot <- df.glt %>%
  pairwise_wilcox_test(relative_area ~ group, p.adjust.method = 'holm') %>%
  add_significance('p.adj') %>%
  add_y_position(fun = 'mean_sd')

# pairwise test between control and treatment groups only
df.area.stat.ctrl <- df.glt %>%
  pairwise_wilcox_test(area ~ group, p.adjust.method = 'holm', ref.group = 'Cont') %>%
  add_significance('p.adj') %>%
  add_y_position(fun = 'mean_sd')

ggplot() +
  geom_boxplot(data = df.glt,
               aes(x = group, y = relative_area, fill = group), alpha = .5) +
  geom_jitter(data = df.glt,                                                                                          
             aes(x = group, y = relative_area), alpha = .5) +
  stat_pvalue_manual(df.area.stat.tot, label = 'p.adj.signif',
                     hide.ns = TRUE, size = 3) +
  theme(legend.position = 'none')


##### CONT #####
df.ctrl <- df.glt %>%
  filter(trauma == 'Cont') %>%
  droplevels()

# AREA
df.ctrl.area.stat <- df.ctrl %>%
  wilcox_test(relative_area ~ treat) %>%
  add_significance('p') %>%
  add_y_position(fun = 'median_iqr')

ggplot() +
  geom_boxplot(data = df.ctrl,
               aes(y = relative_area, x = treat, fill = treat),
               alpha = .5, outlier.shape = NA) +
  geom_point(data = df.ctrl,
             aes(y = relative_area, x = treat, color = treat),
             position=position_jitterdodge()) +
  stat_pvalue_manual(df.ctrl.area.stat, label = 'p.signif') +
  theme(legend.position = 'none')
  
  
##### TREATMENT VS DAYS #####
# AREA
df.profile <- df.glt %>%
  filter(trauma == 'TBI') %>%
  mutate(day = as.factor(day)) %>%
  droplevels()

df.treat.stat <- df.profile %>%
  group_by(day) %>%
  wilcox_test(relative_area ~ treat) %>%
  add_significance('p') %>%
  add_xy_position(x='day', fun = 'median_iqr')

df.day.stat <- df.profile %>%
  pairwise_wilcox_test(relative_area ~ day, p.adjust.method = 'holm') %>%
  add_significance('p.adj') %>%
  add_xy_position(x='day', fun = 'median')

ggplot() +
  geom_boxplot(data = df.profile,
               aes(y = relative_area, x = as.factor(day), fill = treat),
               alpha = .5, outlier.shape = NA) +
  geom_point(data = df.profile,
              aes(y = relative_area, x = as.factor(day), color = treat),
              position=position_jitterdodge()) +
  stat_summary(data = df.profile,
               fun = median,
               geom = 'line',
               aes(y = relative_area, x = as.factor(day), color = treat, group=treat),
               position = position_dodge(width = 0.75)) +
  stat_summary(data = df.profile,
               fun = median,
               geom = 'point',
               aes(y = relative_area, x = as.factor(day), color = treat, group=treat),
               size=3, position = position_dodge(width = 0.75)) +
  stat_pvalue_manual(df.treat.stat, label = 'p.signif',
                     hide.ns = TRUE, tip.length = 0) +
  stat_pvalue_manual(df.day.stat, label = 'p.adj.signif',
                     hide.ns = TRUE)

# INTENSITY
df.int.treat.stat <- df.profile %>%
  group_by(day) %>%
  wilcox_test(relative_intensity ~ treat) %>%
  add_significance('p') %>%
  add_xy_position(x='day', fun = 'median_iqr')

df.int.day.stat <- df.profile %>%
  pairwise_wilcox_test(relative_intensity ~ day, p.adjust.method = 'holm') %>%
  add_significance('p.adj') %>%
  add_xy_position(x='day', fun = 'median')

ggplot() +
  geom_boxplot(data = df.profile,
               aes(y = relative_intensity, x = as.factor(day), fill = treat),
               alpha = .5, outlier.shape = NA) +
  geom_point(data = df.profile,
             aes(y = relative_intensity, x = as.factor(day), color = treat),
             position=position_jitterdodge()) +
  stat_summary(data = df.profile,
               fun = median,
               geom = 'line',
               aes(y = relative_intensity, x = as.factor(day), color = treat, group=treat),
               position = position_dodge(width = 0.75)) +
  stat_summary(data = df.profile,
               fun = median,
               geom = 'point',
               aes(y = relative_intensity, x = as.factor(day), color = treat, group=treat),
               size=3, position = position_dodge(width = 0.75)) +
  stat_pvalue_manual(df.int.treat.stat, label = 'p.signif',
                     hide.ns = TRUE, tip.length = 0) +
  stat_pvalue_manual(df.int.day.stat, label = 'p.adj.signif',
                     hide.ns = TRUE)


##### CORRELATION STAT #####
model.int.vs.area <- lm(relative_intensity ~ relative_area, data = df.glt)
model.int.vs.area

ggplot(data = df.glt,
       aes(x = relative_area, y = relative_intensity, color = group)) +
  stat_smooth(method = "lm", col = "red") +
  geom_point() 

ggplot(data = df.glt,
       aes(x = layer_area, y = dots_area, color = group)) +
  geom_point()

ggplot(data = df.glt,
       aes(x = layer_intensity, y = dots_intensity, color = group)) +
  geom_point()
