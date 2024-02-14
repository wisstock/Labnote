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
  filter(area != "") %>%
  mutate(group = as.factor(group), area = as.numeric(as.character(gsub(',','.',area))))

# df.glt.int <- read.csv('GLTcount_intens_csv.csv', sep = ';')


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

ggplot() +
  geom_boxplot(data = df.glt.area,
               aes(x = group, y = area, fill = group), alpha = .75) +
  geom_point(data = df.glt.area,
             aes(x = group, y = area), alpha = .75) +
  stat_pvalue_manual(df.area.stat.ctrl, label = 'p.adj.signif',
                     hide.ns = TRUE, size = 3) +
  theme(legend.position = 'none')

