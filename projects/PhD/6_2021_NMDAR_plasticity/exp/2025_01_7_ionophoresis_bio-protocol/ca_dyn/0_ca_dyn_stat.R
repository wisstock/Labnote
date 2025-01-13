# NMDA ionophoresis, ionophoresis cloud data analysis for bio-protocol
# Copyright © 2025 Borys Olifirov

require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/ca_dyn')

df <- read.csv('0_df_ca_dyn.csv') %>%
  mutate(id = as.factor(id),
         roi = as.factor(roi),
         index = index - 5) %>%
  select(-X)

font.size <- 15
font.fam <- 'Arial'
box.alpha <- 0.6


##### AMP STAT #####
df.max.amp <- df %>%
  filter(index == 20) %>%
  group_by(id, roi) %>%
  mutate(int_mean = mean(int)) %>%
  select(id, roi, int_mean) %>%
  distinct() %>%
  ungroup()

df.amp.stat <- df.max.amp %>%
  kruskal_test(int_mean ~ roi) %>%
  add_significance() %>%
  mutate(title = paste('H Test p=', p, sep=''),
         roi = factor('Max', c('Max', 'Mid', 'Min')),
         int_mean = 2.0)

df.amp.roi.stat <- df %>%
  filter(index == 20, int < 2.6) %>%
  kruskal_test(int ~ roi) %>%
  add_significance() %>%
  mutate(title = paste('H Test p=', p, sep=''),
         roi = factor('Mid', c('Max', 'Mid', 'Min')),
         int = 2.2)

df.amp.roi.comp <- df %>%
  filter(index == 20, int < 2.6) %>%
  pairwise_wilcox_test(int ~ roi, p.adjust.method = 'BH') %>%
  add_significance()


boxplot_max_amp <- ggplot(data = df %>% filter(index == 20, int < 2.6),
                          aes(x = roi, y = int)) +
  geom_boxplot(aes(fill = roi), alpha = box.alpha) +
  # geom_point(color = 'grey35', size = 1.5) +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'line', size = 0.75, linetype = 'dashed', color = 'grey25') +
  stat_summary(aes(group = id),
               fun = median,
               geom = 'point', size = 1.5, color = 'grey25') +
  geom_text(data = df.amp.roi.stat, aes(label = title), size = 1) +
  scale_fill_manual(values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        axis.text.x = element_text(angle = 90, vjust = 0.7),
        plot.caption = element_text(size = font.size-8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
  # labs(caption = 'n = 2/4/34 (cultures/cells/ROIs)')

save_plot('0_boxplot_max_amp.png', boxplot_max_amp, base_width = 2, base_height = 4, dpi = 300)
remove(boxplot_max_amp)


##### REPRESENTATIVE PROFILES #####
prof_rep <- ggplot(data = df %>% filter(id == '24_05_9_cell01_ch0'),
       aes(x = index, y = int, color = roi, fill = roi, group = roi)) +
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
  annotation_custom(ggplotGrob(boxplot_max_amp), xmin = 25, xmax = 45, ymin = 0.5, ymax = 1.2) +
  scale_fill_manual(name = "ROI type",
                    values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  scale_color_manual(name = "ROI type",
                     values = c('Max' = 'red2', 'Mid' = 'green4', 'Min' = 'blue1')) +
  scale_shape_discrete(name = "ROI position") +
  scale_linetype_discrete(name = "ROI position") +
  theme_classic() +
  theme(legend.position = c(0.91, 0.7),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-8)) +
  scale_x_continuous(breaks = c(-5, 0, 10, 20, 30, 40, 50)) +
  labs(x = 'Time, s',  # caption = 'n = 1/1/9 (cultures/cells/ROIs)',
       y = expression(ΔF/F[0]))

save_plot('0_prof_rep.png', prof_rep, base_width = 5, base_height = 4, dpi = 300)
remove(prof_rep)


