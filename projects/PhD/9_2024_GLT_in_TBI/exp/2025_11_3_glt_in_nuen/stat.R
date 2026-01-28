# GLT in neurons 
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
library(forcats)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggsci)


setwd('/home/wisstock/bio_note/projects/PhD/9_2024_GLT_in_TBI/exp/2025_11_3_glt_in_nuen')

df <- read.csv('df_NeuN.csv') %>%
  mutate_if(is.character, as.factor) %>%
  mutate(group = factor(group, levels = c('Cont', 'TBI_3D', 'TBI_7D', 'TBI_14D')))

df.med <- df %>%
  select(group, treat, relative_area, relative_intensity) %>%
  group_by(group, treat) %>%
  mutate(med_rel_area = median(relative_area),
         med_rel_int = median(relative_intensity)) %>%
  select(-relative_area, -relative_intensity) %>%
  distinct() %>%
  ungroup()

df.summary <- df %>%
  group_by(treat, group) %>%
  summarise(rel_area_median = median(relative_area),
            rel_area_iqr = IQR(relative_area),
            rel_intensity_median = median(relative_intensity),
            rel_intensity_iqr = IQR(relative_intensity),
            slices = n())
write_csv(df.summary, 'summary_neun_tbi_groups.csv')
remove(df.summary)


# plots settings
font.size <- 12  # 19
font.fam <- 'Arial'
box.alpha <- 0.6

cef.color <- 'coral2' 
non.color <- 'deepskyblue3' 

##### CTRL REL AREA #####
ctrl_rel_area_stat <- df %>% filter(group == 'Cont') %>%
  distinct() %>%
  wilcox_test(relative_area ~ treat) %>%
  add_significance() %>%
  add_xy_position()


boxplot_ctrl_rel_area <- ggplot(data = df %>% filter(group == 'Cont'),
       aes(x = fct_relevel(treat, 'none', 'cef'), y = relative_area)) +
  geom_boxplot(aes(fill = treat), , alpha = box.alpha) +
  geom_point(aes(group = treat)) + 
  stat_pvalue_manual(data = ctrl_rel_area_stat, size = font.size-6) +
  scale_fill_manual(name = "Treat",
                    labels = c("Cef.", "Non"),
                    values = c('cef' = cef.color, 'none' = non.color)) +
  scale_x_discrete(labels = c('Non', 'Cef.')) +
  labs(x = 'Control group', y = 'Dots relative area') +
  theme(text=element_text(size = font.size, family = font.fam),
        legend.position = 'None')

boxplot_ctrl_rel_area
save_plot('boxplot_ctrl_rel_area.png', boxplot_ctrl_rel_area,
          base_width = 3.5, base_height = 5, dpi = 300)  # set up plot aspect ratio here
remove(boxplot_ctrl_rel_area, ctrl_rel_area_stat)

##### REL AREA #####
plot_relative_area <- ggplot(data = df %>% filter(group != 'Cont'),
       aes(x = group, y = relative_area,
           color = treat, group = treat)) +
  geom_hline(yintercept = 0.36, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75) +
  scale_color_manual(name = "Treat",
                     labels = c("Cef.", "Non"),
                     values = c('cef' = cef.color, 'none' = non.color)) +
  scale_x_discrete(labels = c('3d afrer TBI', '7d afrer TBI', '14d afrer TBI')) +
  labs(x = 'Group', y = 'Dots relative area') +
  theme(text=element_text(size = font.size, family = font.fam))

plot_relative_area
save_plot('plot_relative_area.png', plot_relative_area,
          base_width = 4.5, base_height = 5, dpi = 300)  # set up plot aspect ratio here


##### CTRL REL INT #####
ctrl_rel_int_stat <- df %>% filter(group == 'Cont') %>%
  distinct() %>%
  wilcox_test(relative_intensity ~ treat) %>%
  add_significance() %>%
  add_xy_position()

boxplot_ctrl_rel_int <- ggplot(data = df %>% filter(group == 'Cont'),
                                aes(x = fct_relevel(treat, 'none', 'cef'), y = relative_intensity)) +
  geom_boxplot(aes(fill = treat), , alpha = box.alpha) +
  geom_point(aes(group = treat)) + 
  stat_pvalue_manual(data = ctrl_rel_int_stat, size = font.size-6) +
  scale_fill_manual(name = "Treat",
                    labels = c("Cef.", "Non"),
                    values = c('cef' = cef.color, 'none' = non.color)) +
  scale_x_discrete(labels = c('Non', 'Cef.')) +
  labs(x = 'Control group', y = 'Dots relative intensity') +
  theme(text=element_text(size = font.size, family = font.fam),
        legend.position = 'None')

boxplot_ctrl_rel_int
save_plot('boxplot_ctrl_rel_int.png', boxplot_ctrl_rel_int,
          base_width = 3.5, base_height = 5, dpi = 300)  # set up plot aspect ratio here
remove(boxplot_ctrl_rel_int, ctrl_rel_int_stat)


##### REL INT #####
plot_relative_int <- ggplot(data = df %>% filter(group != 'Cont'),
       aes(x = group, y = relative_intensity,
           color = treat, group = treat)) +
  geom_hline(yintercept = 0.4, linetype = 'dashed') +
  stat_summary(fun = median,
               geom = 'line', size = 1.25) +
  stat_summary(fun = median,
               geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75) +
  scale_color_manual(name = "Treat",
                     labels = c("Cef.", "Non"),
                     values = c('cef' = cef.color, 'none' = non.color)) +
  scale_x_discrete(labels = c('3d afrer TBI', '7d afrer TBI', '14d afrer TBI')) +
  labs(x = 'Group', y = 'Dots relative intensity') +
  theme(text=element_text(size = font.size, family = font.fam))

plot_relative_int
save_plot('plot_relative_int.png', plot_relative_int,
          base_width = 4.5, base_height = 5, dpi = 300)  # set up plot aspect ratio here
