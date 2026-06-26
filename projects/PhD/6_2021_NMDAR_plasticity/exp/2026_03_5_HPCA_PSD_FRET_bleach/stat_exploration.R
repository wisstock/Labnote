# Sasha's data analysis, HPCA+PSD95 FRET control

library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(ggsci)
library(introdataviz)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2026_03_5_HPCA_PSD_FRET_bleach')


##### DF PREPROCESSING #####
df.full <- read.csv('df_combined.csv') %>%
           select(-X) %>%
           mutate_if(is.character, factor) %>%
           mutate(roi = as.factor(roi),
                  roi_id = interaction(id, roi, sep = '_'))
  
# df.full.id.summary <- df.full %>%
#   filter(base == 'simple', index == 10) %>%
#   group_by(id, ch) %>%
#   summarise(n_roi = n_distinct(roi_id),
#             max = max(dF0_int),
#             min = min(dF0_int))

df.ch <- df.full %>%
  filter(base == 'simple', ch %in% c('ch0', 'ch1', 'ch3'),
         data == 'target', base == 'simple') %>%
  droplevels()

df.fret <- df.full %>%
  filter(base == 'simple', ch %in% c('Fc', 'Ea'),
         data == 'target', base == 'simple') %>%
  droplevels()

remove(df.full, df.full.id.summary)

df.bg.raw <- read.csv('df_snr_test.csv')

df.bg.mean <- df.bg.raw %>%
  select(ch0_mean, ch1_mean, ch3_mean) %>%
  mutate(index = seq(1,120)) %>%
  rename_at(.vars = vars(ends_with("_mean")),
            .funs = funs(sub("[_]mean$", "", .))) %>%
  pivot_longer(cols = 1:3, names_to = "ch", values_to = "val") %>%
  mutate(ch = as.factor(ch), wtf = as.factor('mean'))
  
df.bg.sd <- df.bg.raw %>%
  select(ch0_sd, ch1_sd, ch3_sd) %>%
  mutate(index = seq(1,120)) %>%
  rename_at(.vars = vars(ends_with("_sd")),
            .funs = funs(sub("[_]sd$", "", .))) %>%
  pivot_longer(cols = 1:3, names_to = "ch", values_to = "val") %>%
  mutate(ch = as.factor(ch), wtf = as.factor('sd'))

df.bg <- rbind(df.bg.mean, df.bg.sd)

remove(df.bg.raw, df.bg.mean, df.bg.sd)

##### PLOT PARAMETERS #####
start_box <- seq(0, 1)
end_box <- seq(119, 120)

font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6

ch0.color <- 'chartreuse3'
ch1.color <- 'orange2'
ch3.color <- 'firebrick1'
fc.color <- 'skyblue2'
ea.color <- 'coral2'


##### PANNEL 2 PROFILES INT #####


# abs
ggplot(data = df.ch,
       aes(x = index, y = abs_int, colour = ch, fill = ch)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  # annotate("segment", x = start_box[1], xend = tail(start_box, n = 1), y = -15, yend = -15, size = 5,
  #          colour = 'grey30') +
  # annotate("segment", x = end_box[1], xend = tail(end_box, n = 1), y = -15, yend = -15, size = 5,
  #          colour = 'grey30') +
  # annotate("segment", x = 87.5, y = 115, xend = 87.5, yend = 100, 
  #          size = 8, linejoin = "mitre",
  #          arrow = arrow(type = "closed", length = unit(0.01, "npc"))) +
  # annotate("text", x=87.5,y=115, label = "text", color = "white", 
  #          angle = 90, hjust = 1.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c('ch0' = ch0.color, 'ch1' = ch1.color, 'ch3' = ch3.color),
                    labels = c('DD', 'DA', 'AA'),
                    name = 'Channel') +
  scale_colour_manual(values = c('ch0' = ch0.color, 'ch1' = ch1.color, 'ch3' = ch3.color),
                      labels = c('DD', 'DA', 'AA'),
                      name = 'Channel') +
  theme_classic() +
  theme(legend.position = c(0.9,0.8),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3)) +
  labs(x = 'Time, s',
       y = 'Intensity, a.u.',
       caption = 'n = 1/2/155 (cultures/cells/ROIs)')

# SNR
ch0.bg.val <- df.bg$val[df.bg$ch == 'ch0' & df.bg$wtf == 'mean']
ch0.bg.var <- df.bg$val[df.bg$ch == 'ch0' & df.bg$wtf == 'sd']
ch1.bg.val <- df.bg$val[df.bg$ch == 'ch1' & df.bg$wtf == 'mean']
ch1.bg.var <- df.bg$val[df.bg$ch == 'ch1' & df.bg$wtf == 'sd']
ch3.bg.val <- df.bg$val[df.bg$ch == 'ch3' & df.bg$wtf == 'mean']
ch3.bg.var <- df.bg$val[df.bg$ch == 'ch3' & df.bg$wtf == 'sd']

df.snr <- df.ch %>%
  select(ch, index, abs_int, id) %>%
  group_by(ch, index, id) %>%
  mutate(val = mean(abs_int)) %>%
  select(-abs_int) %>%
  distinct() %>%
  ungroup() %>%
  group_by(ch, index) %>%
  mutate(val = mean(val)) %>%
  select(-id) %>%
  distinct() %>%
  ungroup() %>%
  mutate(snr = if_else())

ggplot(data = df.bg %>% filter(wtf == 'mean'),
       aes(x = index, y = val, color = ch)) +
  geom_line()

ggplot(data = df.snr, 
       aes(x = index, y = val, color = ch)) +
  geom_line()

# dF
ggplot(data = df.ch,
       aes(x = index, y = dF0_int, colour = ch, fill = ch)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  scale_fill_manual(values = c('ch0' = ch0.color, 'ch1' = ch1.color, 'ch3' = ch3.color),
                    labels = c('DD', 'DA', 'AA'),
                    name = 'Channel') +
  scale_colour_manual(values = c('ch0' = ch0.color, 'ch1' = ch1.color, 'ch3' = ch3.color),
                      labels = c('DD', 'DA', 'AA'),
                      name = 'Channel') +
  theme_classic() +
  theme(legend.position = c(0.9,0.8),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3)) +
  labs(x = 'Time, s',
       y = 'Intensity, a.u.',
       caption = 'n = 1/2/155 (cultures/cells/ROIs)')


##### PANNEL 2 PROFILES FRET #####
Ea_scale = 25

ggplot(data = df.fret %>% mutate(abs_int = if_else(ch == "Ea", abs_int * Ea_scale, abs_int)),
       aes(x = index, y = abs_int, colour = ch, fill = ch)) +  #  %>% filter(ch == 'Ea')
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  scale_y_continuous(expression(F[c]), 
                     sec.axis = sec_axis(~ . / Ea_scale, name = expression(E[D]))) +
  scale_fill_manual(values = c('Fc' = fc.color, 'Ea' = ea.color),
                    labels = c(expression(F[c]), expression(E[D])),
                    name = 'FRET') +
  scale_colour_manual(values = c('Fc' = fc.color, 'Ea' = ea.color),
                      labels = c(expression(F[c]), expression(E[D])),
                      name = 'FRET') +
  theme_classic() +
  theme(legend.position = c(0.1,0.25),
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3)) +
  labs(caption = 'n = 1/2/155 (cultures/cells/ROIs)',
       x = 'Time, s')  # y = expression(ΔF/F[0])


##### PANNEL 2 BOXPLOT #####
# raw
df.ch.box <- df.ch %>%
  mutate(box = as.factor(case_when(index %in% start_box ~ 'Start',
                         index %in% end_box ~ 'End',
                         .default = 'x'))) %>%
  filter(box != 'x') %>%
  droplevels() %>%
  group_by(ch, box, roi_id) %>%
  mutate(box = factor(box, levels = c('Start', 'End'), ordered = TRUE),
         med_int = median(abs_int),
         med_dF = median(dF_int),
         med_dF0 = median(dF0_int)) %>%
  ungroup() %>%
  select(ch, box, med_int, med_dF, med_dF0, roi_id) %>%
  distinct()
  
length(levels(df.ch.box$roi_id))


df.ch.box.abs.stat <- df.ch.box %>%
  group_by(ch) %>%
  wilcox_test(med_int ~ box) %>%
  add_significance() %>%
  add_xy_position()


ggplot(data = df.ch.box,
       aes(x = box, y = med_int)) +
  geom_boxplot(aes(fill = ch), alpha = box.alpha) +
  stat_summary(aes(group = ch),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(aes(group = ch),
               fun = median,
               geom = 'line', size = 0.5) +
  stat_pvalue_manual(df.ch.box.abs.stat) +
  facet_wrap(~ch, scale = 'free',
             labeller = as_labeller(c("ch0" = "DD", "ch1" = "DA", "ch3" = "AA"))) +
  theme_classic() +
  theme(legend.position = 'none',
        text=element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size-4),
        legend.title = element_text(size = font.size-3)) +
  scale_fill_manual(values = c('ch0' = ch0.color, 'ch1' = ch1.color, 'ch3' = ch3.color)) +
  labs(caption = 'n = 1/2/155 (cultures/cells/ROIs)',
       x = '',
       y = 'Intensity, a.u.')

###### FRET DEVIATION TEST #####
df.fret.box <- df.fret %>%
  mutate(box = as.factor(case_when(index %in% start_box ~ 'Start',
                                   index %in% end_box ~ 'End',
                                   .default = 'x')),
         box = factor(box, levels = c('Start', 'End'), ordered = TRUE)) %>%
  filter(box != 'x') %>%
  droplevels() %>%
  group_by(ch, box, roi_id) %>%
  mutate(med_int = median(abs_int)) %>%
  ungroup() %>%
  select(ch, roi_id, box, med_int) %>%
  distinct() %>%
  group_by(ch, roi_id) %>%
  mutate(abs_diff = med_int - med_int[box == 'Start']) %>%
  select(ch, box, roi_id, abs_diff) %>%
  ungroup() %>%
  filter(box == 'End') %>%
  droplevels() %>%
  distinct()
  # pivot_wider(names_from = ch, values_from = abs_diff)

length(levels(df.fret.box$roi_id))

df.fret.box.abs.stat <- df.fret.box %>%
  group_by(ch) %>%
  sign_test(abs_diff ~ 1, mu = 0) %>%
  add_significance()


Ea_diff_scale = 100

ggplot(data = df.fret.box %>% mutate(abs_int = if_else(ch == "Ea", abs_diff * Ea_diff_scale, abs_diff)),
         aes(x = ch, y = abs_diff)) +
    geom_boxplot(aes(fill = ch), alpha = box.alpha) +
  scale_y_continuous(expression(F[c]), 
                     sec.axis = sec_axis(~ . / Ea_diff_scale, name = expression(E[D])))
  

##### Fc NOISE ANALYSIS #####
df.fret <- df.full %>%
  filter(base == 'simple', ch %in% c('Fc', 'Ea')) %>%
  droplevels()


ggplot(data = df.fret %>% filter(ch == 'Fc'),
       aes(x = time, y = abs_int, color = data, fill = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  ylim(c(-0.25,5))


df.fret.box <- df.fret %>%
  mutate(box = as.factor(case_when(index %in% start_box ~ 'start',
                                   index %in% end_box ~ 'end',
                                   .default = 'x'))) %>%
  filter(box != 'x') %>%
  droplevels() %>%
  group_by(ch, box, roi_id, data) %>%
  mutate(box = factor(box, levels = c('start', 'end'), ordered = TRUE),
         med_int = median(abs_int),
         med_dF = median(dF_int + 0.000001),
         med_dF0 = median(dF0_int + 0.000001)) %>%
  ungroup() %>%
  select(ch, box, med_int, med_dF, med_dF0, roi_id, data) %>%
  distinct()

ggplot(data = df.fret.box %>% filter(ch == 'Fc'),
       aes(x = box, y = med_int)) +
  geom_boxplot(aes(fill = data)) +
  stat_summary(aes(group = data),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(aes(group = data),
               fun = median,
               geom = 'line', size = 0.5) +
  facet_wrap(~data, scale = 'free')


ggplot(data = df.fret.box %>% filter(ch == 'Ea'),
       aes(x = box, y = med_int)) +
  geom_boxplot(aes(fill = data)) +
  stat_summary(aes(group = data),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(aes(group = data),
               fun = median,
               geom = 'line', size = 0.5) +
  facet_wrap(~data, scale = 'free')


##### 2D #####
df.wide.box <- df.full %>%
  mutate(box = as.factor(case_when(index %in% start_box ~ 'start',
                                   index %in% end_box ~ 'end',
                                   .default = 'x'))) %>%
  filter(box != 'x') %>%
  droplevels() %>%
  filter(base == 'simple', ch %in% c('ch3', 'ch0', 'Fc', 'Ea')) %>%
  select(ch, box, roi_id, data, abs_int) %>%
  group_by(ch, box, roi_id, data) %>%
  mutate(box = factor(box, levels = c('start', 'end'), ordered = TRUE),
         med_int = median(abs_int)) %>%
  ungroup() %>%
  select(ch, box, data, med_int, roi_id) %>%
  distinct() %>%
  pivot_wider(names_from = ch, values_from = med_int)

ggplot(data = df.wide.box,
       aes(x = Ea, y = ch3, colour = data)) +
  geom_point(alpha = .5)
scale_y_continuous(trans = "pseudo_log") +
  scale_x_continuous(trans = "pseudo_log")



df.wide <- df.full %>%
  filter(base == 'simple', ch %in% c('ch3', 'ch0', 'Fc', 'Ea')) %>%
  select(id, roi, index, abs_int, ch, data) %>%
  pivot_wider(names_from = ch, values_from = abs_int) %>%
  mutate(rel = ch0/ch3)


ggplot(data = df.wide %>% filter(Fc != 0, ch3 != 0),
       aes(x = Fc, y = ch3, colour = data)) +
  geom_point(alpha = .5)
  scale_y_continuous(trans = "pseudo_log") +
  scale_x_continuous(trans = "pseudo_log")
  
ggplot(data = df.wide %>% filter(Fc != 0, ch3 != 0),
       aes(x = Fc, y = ch0, colour = data)) +
  geom_point(alpha = .5)
  scale_y_continuous(trans = "pseudo_log") +
  scale_x_continuous(trans = "pseudo_log")
  
  
  ##### EXPLORATION #####
  ggplot(data = df.full %>% filter(base == 'simple'),
         aes(x = time, y = abs_int, color = ch, fill = ch,
             linetype = data, shape = data)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun = median,
                 geom = 'line', linewidth = 0.3) +
    stat_summary(fun = median,
                 geom = 'point', size = 1) +
    stat_summary(fun.min = function(z) {quantile(z,0.25)},
                 fun.max = function(z) {quantile(z,0.75)},
                 fun = median,
                 geom = 'ribbon', linewidth = 0, alpha = .15) +
    theme_minimal() +
    facet_wrap(~id)
  
  
  ggplot(data = df.full %>% filter(base == 'simple', ch == 'ch0'),
         aes(x = time, y = abs_int, color = data, fill = data,
             group = roi_id)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun = median,
                 geom = 'line', linewidth = 0.3) +
    # stat_summary(fun = median,
    #              geom = 'point', size = 1) +
    # stat_summary(fun.min = function(z) {quantile(z,0.25)},
    #              fun.max = function(z) {quantile(z,0.75)},
    #              fun = median,
    #              geom = 'ribbon', linewidth = 0, alpha = .15) +
    theme_minimal() +
    facet_wrap(~id)
