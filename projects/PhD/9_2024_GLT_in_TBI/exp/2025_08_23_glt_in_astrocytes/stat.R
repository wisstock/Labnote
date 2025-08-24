# GLT in astrocytes
#
# COLUMNS FROM PLUGIN
# cell_row = [glt_img.name,                        # id
#             group,                               # group
#             a_region.label,                      # cell_num
#             one_cell_area,                       # cell_area
#             one_dot_num,                         # dot_num
#             one_dot_area,                        # dot_area
#             round(one_dot_rel_area, 3),          # dot_rel_area
#             one_dot_sum_int,                     # dot_sum_int
#             int(one_dot_mean_int),               # dot_mean_int
#             int(one_dot_sum_int / one_dot_num),  # dot_men_int_per_dot
#             int(one_dot_mean_int_dens)]          # dot_mean_int_dens


require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
library(gridExtra)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)


setwd('/home/wisstock/bio_note/projects/PhD/9_2024_GLT_in_TBI/exp/2025_08_23_glt_in_astrocytes')

df <- read.csv('astrocyte_count.csv') %>%
  select(-dot_rel_area, -id) %>%
  mutate_if(is.character, as.factor)


##### CTRL TEST #####
# RAW
df.ctrl <- df %>%
  filter(group == 'cont') %>%
  droplevels() %>%
  select(-group) %>%
  filter(dot_area < 3000, dot_sum_int < 1.5e+6) %>%
  droplevels()


ggplot(data = df.ctrl,
       aes(x = dot_sum_int, y = dot_area)) +
  geom_point(aes(color = treat), alpha = .75) +
  geom_smooth(aes(group = treat), method = 'lm', se = FALSE)


ggplot(data = df.ctrl,
       aes(x = name, y = dot_sum_int, fill = treat)) +
  geom_boxplot() +
  geom_point()

df.ctrl %>%
  group_by(treat) %>%
  kruskal_test(dot_sum_int ~ name) %>%
  add_significance()

# MED
df.ctrl.med <- df.ctrl %>%
  group_by(name) %>%
  mutate(med_dot_sum_int = median(dot_sum_int),
         med_dot_area = median(dot_area),
         med_dot_num = median(dot_num),
         med_cell_area = median(cell_area)) %>%
  select(name, treat, med_dot_sum_int, med_dot_area, med_dot_num) %>%
  distinct() %>%
  droplevels() %>%
  ungroup()


ggplot(data = df.ctrl.med,
       aes(x = treat, y = med_dot_sum_int, fill = treat)) +
  geom_boxplot() +
  geom_point()


ggplot(data = df.ctrl.med,
       aes(x = med_dot_area, y = med_dot_sum_int, color = treat)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)


df.ctrl.med %>%
  wilcox_test(med_dot_sum_int ~ treat) %>%
  add_significance()


##### EXP MED STAT #####
# RAW
df.tbi <-  df %>%
  filter(group != 'cont') %>%
  droplevels() %>%
  filter(dot_area < 3000, dot_sum_int < 1.5e+6) %>%
  droplevels()

ggplot(data = df.tbi,
       aes(x = name, y = dot_sum_int, fill = treat)) +
  geom_boxplot() +
  facet_wrap(~group, scales = "free_x")




###### MED STAT #####
df.tbi.med <- df %>%
  group_by(name) %>%
  mutate(tot_dot_num = sum(dot_num),
         tot_dot_int = sum(dot_num),
         tot_glt = tot_dot_num * tot_dot_int,
         med_dot_sum_int = median(dot_sum_int),
         med_dot_area = median(dot_area),
         med_cell_area = median(cell_area),
         group = factor(group, levels = c('cont', 'tbi3', 'tbi7', 'tbi14'))) %>%
  select(name, group, treat,
         med_cell_area, med_dot_sum_int, med_dot_area,
         tot_dot_num, tot_dot_int, tot_glt) %>%
  distinct() %>%
  droplevels() %>%
  ungroup()

ggplot(data = df.tbi.med,
       aes(x = group, y = tot_glt, fill = treat)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~treat)

df.tbi.med %>%
  group_by(treat) %>%
  pairwise_wilcox_test(med_dot_num ~ group, ref.group = 'cont')
