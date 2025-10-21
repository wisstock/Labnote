# GLT in astrocytes
#
# COLUMNS FROM PLUGIN
# cell_row = [glt_img.name,                        # id
#             group,                               # group
#             a_region.label,                      # cell_num
#             one_cell_area,                       # cell_area
#             one_dot_num,                         # dot_num
#             one_dot_area,                        # dot_area  це сумарна площа на один астроцит, тобто це не можна використати як параметр окремої точки
#             round(one_dot_rel_area, 3),          # dot_rel_area
#             one_dot_sum_int,                     # dot_sum_int
#             int(one_dot_mean_int),               # dot_mean_int
#             int(one_dot_sum_int / one_dot_num),  # dot_men_int_per_dot
#             int(one_dot_mean_int_dens)]          # dot_mean_int_dens


library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggsci)


setwd('/home/wisstock/bio_note/projects/PhD/9_2024_GLT_in_TBI/exp/2025_08_23_glt_in_astrocytes')

df <- read.csv('astrocyte_count.csv') %>%
  select(-id) %>%
  mutate(dot_rel_area = parse_number(dot_rel_area, locale = locale(decimal_mark = ","))) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(dot_rel_area = dot_area / cell_area)


##### CTRL TEST #####
# RAW
df.ctrl <- df %>%
  filter(group == 'cont') %>%
  droplevels() %>%
  select(-group) %>%
  filter(dot_rel_area < 1,
         dot_sum_int < 1.5e+06) %>%
  droplevels()
  
ggplot(data = df.ctrl,
       aes(x = dot_sum_int, y = dot_rel_area,
           color = treat, fill = treat)) +
  geom_point(alpha = .3) +
  geom_smooth(method = 'loess', se = TRUE, span = 0.85)


ggplot(data = df.ctrl,
       aes(x = dot_rel_area, fill = treat)) +
  geom_density(alpha = .3) +
  scale_x_continuous(limits = c(0,1))




ggplot(data = df.ctrl,
       aes(x = name, y = dot_sum_int, fill = treat)) +
  geom_boxplot() +
  geom_point()

df.ctrl %>%
  group_by(treat) %>%
  kruskal_test(dot_sum_int ~ name) %>%
  add_significance()

df.ctrl %>%
  wilcox_test(dot_sum_int ~ treat) %>%
  add_significance()


ggplot(data = df.ctrl,
       aes(x = name, y = dot_rel_area, fill = treat)) +
  geom_boxplot() +
  geom_point()

df.ctrl %>%
  group_by(treat) %>%
  kruskal_test(dot_rel_area ~ name) %>%
  add_significance()

df.ctrl %>%
  wilcox_test(dot_rel_area ~ treat) %>%
  add_significance()


# dot mean
ggplot(data = df.ctrl,
       aes(x = name, y = dot_mean_int_dens, fill = treat)) +
  geom_boxplot() +
  geom_point()

df.ctrl %>%
  wilcox_test(dot_mean_int_dens ~ treat) %>%
  add_significance()

# MED
df.ctrl.med <- df.ctrl %>%
  group_by(name) %>%
  mutate(med_dot_mean_int = median(dot_mean_int_dens),
         med_dot_sum_int = median(dot_sum_int),
         med_dot_area = median(dot_rel_area),
         med_dot_num = median(dot_num)) %>%
  select(name, treat, med_dot_sum_int, med_dot_area, med_dot_num, med_dot_mean_int) %>%
  distinct() %>%
  droplevels() %>%
  ungroup()

ggplot(data = df.ctrl.med,
       aes(x = treat, y = med_dot_sum_int, fill = treat)) +
  geom_boxplot() +
  geom_point()

df.ctrl.med %>%
  wilcox_test(med_dot_sum_int ~ treat) %>%
  add_significance()


ggplot(data = df.ctrl.med,
       aes(x = treat, y = med_dot_area, fill = treat)) +
  geom_boxplot() +
  geom_point()

df.ctrl.med %>%
  wilcox_test(med_dot_area ~ treat) %>%
  add_significance()



ggplot(data = df.ctrl.med,
       aes(x = treat, y = med_dot_mean_int, fill = treat)) +
  geom_boxplot() +
  geom_point()

df.ctrl.med %>%
  wilcox_test(med_dot_mean_int ~ treat) %>%
  add_significance()


ggplot(data = df.ctrl.med,
       aes(y = med_dot_area, x = med_dot_sum_int, color = treat)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)





##### EXP MED STAT #####
# RAW
df.tbi <-  df %>%
  filter(group != 'cont') %>%
  droplevels()

ggplot(data = df.tbi,
       aes(x = name, y = dot_sum_int, fill = treat)) +
  geom_boxplot() +
  facet_wrap(~group, scales = "free_x")

ggplot(data = df.tbi,
       aes(x = name, y = dot_rel_area, fill = treat)) +
  geom_boxplot() +
  facet_wrap(~group, scales = "free_x")

ggplot(data = df.tbi,
       aes(x = name, y = dot_mean_int_dens, fill = treat)) +
  geom_boxplot() +
  facet_wrap(~group, scales = "free_x")


# ggplot(data = df.tbi,
#        aes(x = dot_sum_int, y = dot_rel_area)) +
#   geom_point(aes(color = treat, fill = treat), alpha = .3) +
#   geom_smooth(aes(group = treat), method = 'loess', se = TRUE, span = 0.85) +
#   facet_wrap(~group, scales = "free_x")




###### MED STAT #####
df.tbi.med <- df.tbi %>%
  group_by(name) %>%
  mutate(med_dot_mean_int = median(dot_mean_int_dens),
         med_dot_sum_int = median(dot_sum_int),
         med_dot_area = median(dot_rel_area),
         med_dot_num = median(dot_num)) %>%
  select(name, treat, group, med_dot_sum_int, med_dot_area, med_dot_num, med_dot_mean_int) %>%
  distinct() %>%
  droplevels() %>%
  ungroup() %>%
  filter(!(med_dot_area < 0.1 & group == 'tbi14')) %>%
  mutate(group = factor(group, levels = c('cont', 'tbi3', 'tbi7', 'tbi14')))

df.ctrl.non.med <- df.ctrl.med %>%
  filter(treat != 'cef') %>%
  mutate(group = 'cont')

df.ctrl.non.med.2nd <- df.ctrl.med %>%
  filter(treat != 'cef') %>%
  mutate(group = 'cont', treat = 'cef')

ctrl.non.val <- df.ctrl.med %>%
  filter(treat != 'cef') %>%
  select(where(is.numeric), -med_dot_num) %>%
  summarise_all(list(median, IQR))

df.tbi.ctrl.med <- rbind(df.tbi.med, df.ctrl.non.med, df.ctrl.non.med.2nd) %>%
  mutate(group = factor(group, levels = c('cont', 'tbi3', 'tbi7', 'tbi14')))

ggplot(data = df.tbi.ctrl.med,
       aes(x = group, y = med_dot_sum_int, fill = treat)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~treat)

ggplot(data = df.tbi.ctrl.med,
       aes(x = group, y = med_dot_mean_int, fill = treat)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~treat)

ggplot(data = df.tbi.med,
       aes(x = group, y = med_dot_area, fill = treat)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~treat)


###### STAT CTRL vs GROUP #####
df.tbi.ctrl.med %>%
  group_by(treat) %>%
  pairwise_wilcox_test(med_dot_sum_int ~ group, ref.group = 'cont')

df.tbi.ctrl.med %>%
  group_by(treat) %>%
  pairwise_wilcox_test(med_dot_mean_int ~ group, ref.group = 'cont')

df.tbi.ctrl.med %>%
  group_by(treat) %>%
  pairwise_wilcox_test(med_dot_area ~ group, ref.group = 'cont')


###### VAL vs TREAT in TBI #####
df.tbi.med %>%
  group_by(group) %>%
  wilcox_test(med_dot_sum_int ~ treat) %>%
  add_significance()

df.tbi.med %>%
  group_by(group) %>%
  wilcox_test(med_dot_area ~ treat) %>%
  add_significance()

df.tbi.med %>%
  group_by(group) %>%
  wilcox_test(med_dot_mean_int ~ treat) %>%
  add_significance()

###### TIME LINES #####
# mean dot int
ggplot(data = df.tbi.med,
       aes(x = group, y = med_dot_mean_int,
           color = treat, group = treat)) +
  geom_hline(yintercept = 354.75, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75)

# mean sum int
ggplot(data = df.tbi.med,
       aes(x = group, y = med_dot_sum_int,
           color = treat, group = treat)) +
  geom_hline(yintercept = 138463.5, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75)

# mean sum int
ggplot(data = df.tbi.med,
       aes(x = group, y = med_dot_area,
           color = treat, group = treat)) +
  geom_hline(yintercept = 0.23, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75)



###### TBI LINE PLOT ######

ggplot(data = df.tbi %>% filter(group == 'tbi3'),
       aes(x=dot_rel_area, y = dot_sum_int, color = treat)) +
  geom_point(alpha = .1) +
  geom_smooth(method = 'lm', se = FALSE)
