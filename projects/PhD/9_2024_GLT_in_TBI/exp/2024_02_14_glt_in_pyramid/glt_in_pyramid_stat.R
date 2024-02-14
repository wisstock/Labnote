# GLT count in hippocampus pyramid layer
# Copyright © 2023 Borys Olifirov

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
df.glt.area <- read.csv('GLTcount_csv.csv', sep = ';') %>%
  pivot_longer(cols = Cont:TBI_14D_cef,
               names_to = 'group', values_to = 'area') %>%
  filter(area != '') %>%
  mutate(group = as.factor(group), area = as.numeric(as.character(gsub(',','.',area)))) %>%
  mutate(trauma = as.factor(ifelse(grepl('Cont', group),'Cont','TBI'))) %>%
  mutate(treat = as.factor(ifelse(grepl('cef', group),'Cef','Placebo'))) %>%
  mutate(day = case_when((group == 'Cont') | (group == 'Cont_cef') ~ 0,
                         (group == 'TBI_3D') | (group == 'TBI_3D_cef') ~ 3,
                         (group == 'TBI_7D') | (group == 'TBI_7D_cef') ~ 7,
                         (group == 'TBI_14D') | (group == 'TBI_14D_cef') ~ 14))
  

df.glt.int <- read.csv('GLTcount_intens_csv.csv', sep = ';') %>%
  mutate(relative.area = as.numeric(as.character(gsub(',','.',relative.area))),
         relative.int = as.numeric(as.character(gsub(',','.',relative.int))))


##### AREA STAT #####
kruskal.test(area ~ group, data = df.glt.area)

# pairwise test between all groups
df.area.stat.tot <- df.glt.area %>%
  pairwise_wilcox_test(area ~ group, p.adjust.method = 'holm') %>%
  add_significance('p.adj') %>%
  add_y_position(fun = 'mean_sd')

# pairwise test between control and treatment groups only
df.area.stat.ctrl <- df.glt.area %>%
  pairwise_wilcox_test(area ~ group, p.adjust.method = 'holm', ref.group = 'Cont') %>%
  add_significance('p.adj') %>%
  add_y_position(fun = 'mean_sd')

# simple boxplot
ggplot() +
  geom_boxplot(data = df.glt.area,
               aes(x = group, y = area, fill = group), alpha = .5) +
  geom_point(data = df.glt.area,                                                                                          
             aes(x = group, y = area), alpha = .5) +
  stat_pvalue_manual(df.area.stat.ctrl, label = 'p.adj.signif',
                     hide.ns = TRUE, size = 3) +
  theme(legend.position = 'none')

# time profile with boxplots
df.profile <- df.glt.area %>%
  filter(trauma == 'TBI') %>%
  droplevels()

ctrl.median <- median(df.glt.area$area[df.glt.area$trauma == 'Cont' & df.glt.area$treat == 'Placebo'])

ggplot() +
  geom_hline(yintercept = ctrl.median, linetype='dashed') +
  geom_point(data = df.profile,
               aes(y = area, x = day, color = treat),
               alpha = .5) +
  stat_summary(data = df.profile,
               fun = median,
               geom = 'line',
               aes(y = area, x = day, color = treat)) +
  stat_summary(data = df.profile,
               fun = median,
               geom = 'point',
               aes(y = area, x = day, color = treat), size=3)


##### CORRELATION STAT #####
model.int.vs.area <- lm(relative.int ~ relative.area, data = df.glt.int)
model.int.vs.area

ggplot(data = df.glt.int,
       aes(x = relative.area, y = relative.int)) +
  stat_smooth(method = "lm", col = "red") +
  geom_point() 
