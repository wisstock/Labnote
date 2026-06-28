# Sasha's data analysis, HPCA+PSD95

require(dplyr)
require(tidyr)
require(purrr)
require(forcats)
require(rstatix)
library(gridExtra)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)
library(ggmagnify)
library(introdataviz)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_06_19_Sasha_HPCA-PSD')

# PLOT SETTINGS
font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6

psd.color <- 'magenta3'
halo.color <- 'green4'

##### DATA LOADING AND PREPROCESSING #####
df.full <- read.csv('data/HPCA/df.csv') %>%
           select(-X) %>%
           mutate_if(is.character, factor) %>%
           mutate(roi = as.factor(roi),
                  app_factor = as.factor(app),
                  dist_um = dist * 0.16,
                  lab_id = fct_recode(lab_id, 'halo' = 'oreol'),
                  dist_group = as.factor(case_when(dist_um <= 25 ~ 'I',
                                                   (dist_um > 25) & (dist_um <= 50)  ~ 'II',
                                                   dist_um > 50 ~ 'III',
                                                   .default = '0'))) %>%
          filter(app_factor %in% c('0.5','2.5', '5', '10', '20', '30')) %>%
          mutate(rel_time = index - 30)


spines_id <- df.full %>%
  filter(lab_id != 'shaft') %>%
  droplevels() %>%
  mutate(roi_id = interaction(roi, id, sep = '_')) %>%
  group_by(roi_id) %>%
  summarise(has_all = n_distinct(lab_id) == 2) %>%
  filter(has_all) %>%
  droplevels() %>%
  pull(roi_id)

df.spines <- df.full %>%
  filter(lab_id != 'shaft') %>%
  select(-app, -dist) %>%
  droplevels() %>%
  mutate(roi_id = interaction(roi, id, sep = '_')) %>%
  filter(roi_id %in% spines_id) %>%
  group_by(roi_id, index) %>%
  mutate(dist_um = dist_um[lab_id == 'psd'],
         int_diff = abs_int[lab_id == 'halo'] - abs_int[lab_id == 'psd']) %>%
  ungroup() %>%
  mutate(dist_group = as.factor(case_when(dist_um <= 25 ~ 'I',
                                          (dist_um > 25) & (dist_um <= 50)  ~ 'II',
                                          dist_um > 50 ~ 'III',
                                          .default = '0')))

df.shaft <- df.full %>%
  filter(lab_id == 'shaft') %>%
  select(-app, -dist) %>%
  mutate(roi_id = interaction(roi, id, sep = '_shaft_'),
         int_diff = 0)

df.preproc <- rbind(df.spines, df.shaft) %>%
  filter(lab_id != 'shaft') %>%
  droplevels()
remove(df.shaft, df.spines, spines_id)


##### PROFILES FOR 20 ANALYSIS #####
df.20 <-  df.preproc %>%
  filter(app_factor == '20') %>%
  mutate(rel_time = index - 30)

# df
# df.from <- c(xmin = 1, xmax = 4, ymin = -0.001, ymax = 0.075)
# df.to <- c(45, 75, 0.05, 0.15)

plot.int <- ggplot(data = df.20,  #  %>% filter(dist_group != 'I')
       aes(x = rel_time, y = df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_vline(xintercept = 20, linetype = 3) +
  stat_summary(aes(color = lab_id, group = interaction(id,lab_id)),
               fun = median,
               geom = 'line', linewidth = 0.25) +
  stat_summary(aes(color = lab_id, group = lab_id),
               fun = median,
               geom = 'line', linewidth = 0.75) +
  stat_summary(aes(color = lab_id, group = lab_id),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(aes(fill = lab_id, group = lab_id),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .25) +
  scale_y_continuous(breaks = seq(-0.1, 0.4, 0.05)) +
  scale_color_manual(name = "Region",
                     values = c('halo' = halo.color, 'psd' = psd.color, 'shaft' = shaft.color),
                     labels = c('Halo', 'PSD')) +
  scale_fill_manual(name = "Region",
                    values = c('halo' = halo.color, 'psd' = psd.color, 'shaft' = shaft.color),
                    labels = c('Halo', 'PSD')) +
  labs(x = 'Time, s',
       y = 'ΔF/F0',
       caption = 'n = 2/5/372 (cultures/cells/ROIs)') +
  theme_classic() +
  theme(legend.position = c(0.075,0.8),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3))
  # geom_magnify(from = df.from, to = df.to,
  #              axes = "y")

plot.int
save_plot('plot_int.png', plot.int, base_width = 9.2, base_height = 4)



##### LM AGREGATE PLOT #####
# selected frames lm
time.points <- c(0, 1, 3, 5)


# # loess
# ggplot(data = df.preproc %>%
#          filter(rel_time %in% time.points, app_factor == '20') %>%
#          mutate(lab_dist = interaction(lab_id, dist_group, sep = '_')),
#        aes(x = dist_um, y = df,
#            color = lab_id, fill = lab_id)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_point(alpha = .15) +
#   # geom_line(aes(group = roi_id), alpha = .15, color = 'grey') +
#   geom_smooth(method = 'loess', alpha = .15) +
#   scale_y_continuous(limits = c(-0.5, 0.5)) +
#   facet_wrap(~as.factor(rel_time), ncol = 3) +  # ~interaction(as.factor(rel_time), app_factor, sep = '_')
#   labs(title = 'Distance vs ΔF/F0 (app. 20 s)',
#        x = 'Distance (μm)',
#        y = 'ΔF/F0') +
#   scale_color_manual(name = "Spine region",
#                      values = c('halo' = halo.color, 'psd' = psd.color, 'shaft' = shaft.color)) +
#   scale_fill_manual(name = "Spine region",
#                     values = c('halo' = halo.color, 'psd' = psd.color, 'shaft' = shaft.color)) +
#   theme_minimal() +
#   theme(text = element_text(size = 25))


# with lm for intervals
plot.lm <- ggplot(data = df.preproc %>%
         filter(rel_time %in% time.points, app_factor == '20') %>%
         mutate(lab_dist = interaction(lab_id, dist_group, sep = '_')),
       aes(x = dist_um, y = df,
           color = lab_id, fill = lab_id, group = lab_dist)) +
  geom_vline(xintercept = c(25, 50), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = .1) +
  geom_smooth(method = 'lm', linewidth = 1) +
  scale_y_continuous(limits = c(-0.25, 0.3)) +
  scale_color_manual(values = c('psd' = psd.color, 'halo' = halo.color)) +
  scale_fill_manual(values = c('psd' = psd.color, 'halo' = halo.color)) +
  facet_wrap(~rel_time, ncol = 2,
             labeller = as_labeller(c('0' = '0 s after app. start',
                                      '1' = '1 s after app. start',
                                      '3' = '3 s after app. start',
                                      '5' = '5 s after app. start'))) +
  labs(x = 'Distance from electrode tip, μm',
       y = 'ΔF/F0',
       caption = 'n = 2/5/372 (cultures/cells/ROIs)') +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3))

plot.lm
save_plot('plot_lm.png', plot.lm, base_width = 9.2, base_height = 4)


##### ROI BOX #####
df.box.stat <- df.20 %>%
  filter(rel_time %in% time.points) %>%
  group_by(rel_time, dist_group) %>%
  wilcox_test(df ~ lab_id) %>%
  add_significance() %>%
  add_xy_position() %>%
  mutate(groups = dist_group,
         group1 = dist_group,
         group2 = dist_group)
  
plot.box <- ggplot(data = df.20 %>% filter(rel_time %in% time.points),
       aes(x = dist_group, y = df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_split_violin(aes(fill = lab_id), alpha = box.alpha) +
  # geom_boxplot(alpha = box.alpha) +
  # geom_point(alpha = .15) +
  # geom_line(aes(group = roi_id), linewidth = .1) +
  stat_summary(aes(group = lab_id, color = lab_id),
               fun = median,
               geom = 'point', size = 3) +
  stat_pvalue_manual(df.box.stat, x = 'dist_group',
                     remove.bracket = TRUE, hide.ns = TRUE) +
  scale_color_manual(values = c('halo' = halo.color, 'psd' = psd.color, 'shaft' = shaft.color)) +
  scale_fill_manual(values = c('halo' = halo.color, 'psd' = psd.color, 'shaft' = shaft.color)) +
  scale_x_discrete(labels=c("I" = "<25 μm", "II" = "25-50 μm",
                            "III" = ">50 μm")) +
labs(x = 'Distance from electrode tip',
     y = 'ΔF/F0',
     caption = 'n = 2/5 (cultures/cells) <25 μm = 279, 25-50 μm = 402, >50 μm = 191 (ROIs)') +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3)) +
  facet_wrap(~rel_time, ncol = 2,
             labeller = as_labeller(c('0' = '0 s after app. start',
                                      '1' = '1 s after app. start',
                                      '3' = '3 s after app. start',
                                      '5' = '5 s after app. start')))

plot.box
save_plot('plot_box.png', plot.box, base_width = 9.2, base_height = 4)

##### ROI PHASE PORTRAIT STAT #####
df.by.roi <- df.preproc %>%
  filter(app_factor %in% c('20'), dist_group %in% c('I', 'II', 'III'), lab_id != 'shaft') %>%
  select(roi_id, df, lab_id, rel_time, app_factor, dist_group) %>%
  pivot_wider(names_from = lab_id, values_from = df) %>%
  drop_na() %>%
  group_by(roi_id) %>%
  mutate(filter_group = ifelse(((rel_time == 2) & (halo > psd)), TRUE, FALSE),
         rise_group = ifelse(!all(filter_group == FALSE), 'up', 'down')) %>%
  select(-filter_group) %>%
  droplevels() %>%
  ungroup()


track.time <- seq(-10,40)

df.avg.roi <- df.by.roi %>%
  group_by(rel_time) %>%
  summarise(psd_mean = mean(psd, na.rm = TRUE),
            oreol_mean = mean(halo, na.rm = TRUE),
            psd_se = sd(psd, na.rm = TRUE) / sqrt(n()),
            oreol_se = sd(halo, na.rm = TRUE) / sqrt(n()),
            n = n(),
            .groups = 'drop')


plot.phase <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_path(data = df.avg.roi %>% filter(rel_time %in% track.time),
            aes(x = psd_mean, y = oreol_mean, color = rel_time), size = 2) +
  geom_point(data = df.avg.roi %>% filter(rel_time %in% track.time),
             aes(x = psd_mean, y = oreol_mean, color = rel_time), size = 3) +
  geom_linerange(data = df.avg.roi %>% filter(rel_time %in% track.time),
                aes(x = psd_mean, y = oreol_mean,
                    xmin = psd_mean - psd_se*1.96,
                    xmax = psd_mean + psd_se*1.96,
                    color = rel_time)) +
  geom_linerange(data = df.avg.roi %>% filter(rel_time %in% track.time),
                aes(x = psd_mean, y = oreol_mean,
                    ymin = oreol_mean - oreol_se*1.96,
                    ymax = oreol_mean + oreol_se*1.96,
                    color = rel_time)) +
  scale_color_gradient2(low = 'blue2', mid = 'orange', high = "red2",
                        midpoint = 15, limits = c(-10, 40)) +
  labs(x = 'PSD ΔF/F0',
       y = 'Halo ΔF/F0',
       color = 'Time, s',
       caption = 'n = 2/5/372 (cultures/cells/ROIs)') +
  theme_classic() +
  theme(legend.position = c(0.1, 0.5),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3))

plot.phase
save_plot('plot_phase.png', plot.phase, base_width = 9.2, base_height = 4)
