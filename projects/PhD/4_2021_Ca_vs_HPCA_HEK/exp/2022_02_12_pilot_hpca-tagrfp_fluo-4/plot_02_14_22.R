# HPCA-TagRFP + Fluo-4 data analysis
# Copyright © 2022 Borys Olifirov
#
# Profile data frame structure:
# 'ID' - recording ID
# 'power' - 405 nm stimulation power (%)
# 'ch' - channel (FP or Ca dye)
# 'frame' - frame num
# 'time' - frame time (s)
# 'mask' - mask type (master, up, down)
# 'mask_region' - mask region (1 for master or down)
# 'mean' - mask mean intensity
# 'delta' - mask ΔF/F
# 'rel' - mask mean / master mask mean
#
# Area data frame structure:
# 'ID' - recording ID
# 'stim_frame' - stimulation frame number
# 'mask' - mask type (up or down)
# 'area' - mask region area (in px)
# 'rel_area' - relative mask area in master mask region only (mask / master mask)

require(dplyr)
require(tidyr)
require(purrr)

require(rstatix)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/astria/bio/note/projects/PhD/6_2021_PIP2_HPCA_HEK/exp/2022_02_12_pilot_hpca-tagrfp_fluo-4')

##### GLOBAL PLOT OPTIONS #####
font_family <- "ubuntu mono"
font_size <- 15

line_size <- 1
line_size_light <- 1.5


##### DATA PREPROCESSING #####
# profile section
df.total <- read.csv('profile_02_2_2022.csv') %>%
  mutate(ID = as.factor(ID),
         mask = as.factor(mask),
         mask_region = as.factor(mask_region))
df.ca <- df.total %>%
         filter(ch == 'ca')
df.fp <- df.total %>%
         filter(ch == 'fp')
remove(df.total)

# up mask region selection
membrane.region <- list('cell1_02_2_2022' = c(1),
                        'cell2_02_2_2022' = c(1, 2),
                        'cell4_02_2_2022' = c(1, 2, 4, 5),
                        'cell5_02_2_2022' = c(1),
                        'cell6_02_2_2022' = c(2),
                        'cell7_02_2_2022' = c(1, 3, 5))

df.fp.up.region <- df.fp %>%
  data.frame() %>%
  filter(mask == 'up')
df.fp.cyto.region <- df.fp.up.region %>%
  data.frame()

for (cell.id in names(membrane.region)) {
  df.fp.up.region <- df.fp.up.region %>%
    filter(!(ID == cell.id & !(mask_region %in% membrane.region[[cell.id]]))) %>%
    select(ID, frame, ch, time, mask_region, mean, delta, rel) %>%
    mutate(loc = 'm')
  
  df.fp.cyto.region <- df.fp.cyto.region %>%
    filter(!(ID == cell.id & mask_region %in% membrane.region[[cell.id]])) %>%
    select(ID, frame, ch, time, mask_region, mean, delta, rel) %>%
    mutate(loc = 'c')
}

df.fp.up.region <- rbind(df.fp.up.region, df.fp.cyto.region) %>%
  mutate(loc = as.factor(loc))
remove(cell.id, df.fp.cyto.region)

# px-wise up mask intensity section
df.px <- read.csv('px_02_2_2022.csv') %>%
  mutate(ID = as.factor(ID),
         mask_region = as.factor(mask_region),
         stim = as.factor(stim))

df.px.cyto <- df.px %>% data.frame()
for (cell.id in names(membrane.region)) {
  df.px <- df.px %>%
    filter(!(ID == cell.id & !(mask_region %in% membrane.region[[cell.id]]))) %>%
    mutate(loc = 'm')
  
  df.px.cyto <- df.px.cyto %>%
    filter(!(ID == cell.id & mask_region %in% membrane.region[[cell.id]])) %>%
    mutate(loc = 'c')
}
df.px <- rbind(df.px, df.px.cyto) %>%
  mutate(loc = as.factor(loc))
remove(df.px.cyto, membrane.region)

# area section
df.area <- read.csv('area_02_2_2022.csv') %>%
  mutate(ID = as.factor(ID),
         stim_num = as.factor(stim_num),
         stim_frame = as.factor(stim_frame),
         mask = as.factor(mask))

# translocation section
# in progress

##### AREA ANALYSIS #####
df.area.crop <- df.area %>%
  filter(rel_area != 0, mask != 'master') %>%
  select(ID, mask, stim_num, rel_area) %>%
  droplevels()

df.area.stat <- df.area.crop %>%
  group_by(mask) %>%
  kruskal_test(rel_area ~ stim_num)

area.plot <- ggplot(df.area.crop,
       aes(x = stim_num, y = rel_area*100,
           color = ID, group = ID)) +
  geom_point() +
  facet_wrap(~mask, scales='free_y', nrow = 1) +
  labs(x = 'Stimulation number', y = 'Relative mask area (%)') +
  scale_color_startrek() +
  scale_fill_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(1, .38),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('mask_area.png', area.plot, 
       width = 21, height = 15, units = 'cm')

remove(df.area.crop, df.area.stat, area.plot)


##### UP REGIONS PX-WISE #####
# https://stats.stackexchange.com/questions/539011/suitable-approach-to-cluster-histogram-like-dataset-using-hdbscan-implementation

# 1st vs. 3d stims
df.px.int <- df.px %>%
  select(-delta) %>%
  group_by(ID, stim) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = stim, values_from = int) %>%
  select(-row) %>%
  rename(stim_1 = '1', stim_2 = '2', stim_3 = '3') %>%
  ungroup()

df.px.delta <- df.px %>%
  select(-int) %>%
  group_by(ID, stim) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = stim, values_from = delta) %>%
  select(-row) %>%
  rename(stim_1 = '1', stim_2 = '2', stim_3 = '3') %>%
  ungroup()

# total delta vs. delta    
total.delta.scatter <- ggplot(df.px.delta, aes(x = stim_1,
                                               y = stim_3,
                                               color = loc,
                                               shape = ID)) +
  geom_point(alpha = .75, size = 2) +
  labs(x = '1st stim. ΔF/F0', y = '3d stim. ΔF/F0') +
  scale_color_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family, face="bold", size=font_size))
ggsave('delta_vs_delta_total.png',
       total.delta.scatter, 
       width = 22, height = 15, units = 'cm')
remove(total.delta.scatter)

# total int vs. int
total.int.scatter <- ggplot(df.px.int, aes(x = stim_1,
                                           y = stim_3,
                                           color = loc,
                                           shape = ID)) +
  geom_point(alpha = .5, size = 2) +
  labs(x = '1st stim. int, AU', y = '3d stim. int, AU') +
  scale_color_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family, face="bold", size=font_size))
ggsave('int_vs_int_total.png',
       total.int.scatter, 
       width = 22, height = 15, units = 'cm')
remove(total.int.scatter)

delta.cell <- ggplot(df.px.delta %>% filter(ID == 'cell5_02_2_2022'),
       aes(x = stim_1, y = stim_3, color = loc, shape = mask_region)) +
  geom_point(alpha = .5, size = 2) +
  labs(x = '1st stim. ΔF/F0', y = '3d stim. ΔF/F0') +
  scale_color_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family, face="bold", size=font_size))
ggsave('delta_vs_delta_cell5.png',
       delta.cell, 
       width = 22, height = 15, units = 'cm')
remove(delta.cell)


ggplot(demo.px, aes(x = stim_1, y = stim_3, color = loc, shape = ID)) +
  geom_point(alpha = .5, size = 2)
  

# density plot
for (rec.id in levels(df.px$ID)) {
  region.density <- ggplot(df.px %>% filter(ID == rec.id, stim == 3), aes(x = delta,
                                             fill = mask_region)) +
    geom_density(alpha = .5, color = 'black') +
    scale_x_continuous(name = 'ΔF/F0',
                       limits = c(0, 1.05),
                       breaks = seq(0, 100, .1),
                       expand = c(0,0)) +
    scale_y_continuous(name = 'Density',
                       breaks = NULL,
                       expand = c(0,0)) +
    scale_fill_simpsons() +
    theme_classic() +
    theme(panel.grid.major = element_line(linetype = 'dotted',
                                          size = 0.1,
                                          colour = "black"),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          text = element_text(family=font_family, face="bold", size=font_size))

  ggsave(sprintf('density_stim3_%s.png', rec.id),
         region.density, 
         width = 12, height = 15, units = 'cm')  
}
remove(rec.id, region.density)

##### Ca PROFILE #####
ca.profile.ctrl <- ggplot(df.ca, aes(x = time, y = delta, color = ID)) +
  geom_line(size = line_size) +
  geom_point() +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(0, 60),
                     breaks = seq(0, 100, 5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.1, 3.5),
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  scale_color_startrek() +
  scale_fill_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.175, 1),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('ca_profile_ctrl.png', ca.profile.ctrl, 
       width = 28, height = 12, units = 'cm')
remove(ca.profile.ctrl)
 

##### BLEACHING CTRL #####
master.ctrl <- ggplot() +
  geom_line(data = df.fp %>% filter(mask == 'master'),
            aes(x = time,y = delta*10, color = ID),
            linetype = "dashed", size = line_size) +
  geom_point(data = df.fp %>% filter(mask == 'master'),
             aes(x = time,y = delta*10, color = ID)) +
  geom_line(data = df.ca,
            aes(x = time, y = delta, color = ID),
            size = line_size) +
  geom_point(data = df.ca,
             aes(x = time, y = delta, color = ID)) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(0, 60),
                     breaks = seq(0, 100, 5),
                     expand = c(0,0)) +
  scale_y_continuous(name = "Fluo-4 ΔF/F0",
                     limits = c(-1.75, 3.5),
                     breaks = seq(-100, 100, 1),
                     sec.axis = sec_axis(trans=~./10, name="HPCA-TagRFP ΔF/F0")) +
  scale_color_startrek() +
  scale_fill_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.18, 1),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))
ggsave('master_mask_ctrl.png', master.ctrl, 
       width = 28, height = 12, units = 'cm')
remove(master.ctrl)


##### UP PROFILE PLOT #####
# # all up regions
# df.fp.up.total <- df.fp %>%
#   filter(mask == 'up') %>%
#   select(c(ID, time, delta)) %>%
#   group_by(ID, time) %>%
#   summarise_all(list(mean, sd), na.rm = TRUE) %>%
#   rename(mean = fn1, sd = fn2) %>%
#   ungroup()

df.fp.up.total <- df.fp.up.region %>%
  select(c(ID, time, delta)) %>%
  group_by(ID, time) %>%
  summarise_all(list(mean, sd), na.rm = TRUE) %>%
  rename(mean = fn1, sd = fn2) %>%
  ungroup()
  
up.profile.ctrl <- ggplot(df.fp.up.total, aes(x = time, y = mean, color = ID)) +
  geom_line(size = line_size) +
  geom_point() +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(0, 60),
                     breaks = seq(0, 100, 5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.1, 0.61),
                     breaks = seq(-100, 100, .1),
                     expand = c(0,0)) +
  scale_color_startrek() +
  scale_fill_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.175, 1),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('up_profile_ctrl.png', up.profile.ctrl, 
       width = 28, height = 12, units = 'cm')
remove(up.profile.ctrl)


##### NORM Ca + UP #####
# density plot
for (rec.id in levels(df.px$ID)) {
  norm.int <- ggplot() +
    geom_line(data = df.fp.up.total %>% filter(ID == rec.id) %>% mutate(mean = mean/max(mean)),
              aes(x = time, y = mean),
              size = line_size, linetype = "dashed",) +
    geom_line(data = df.ca %>% filter(ID == rec.id) %>% mutate(delta = delta/max(delta)),
              aes(x = time, y = delta),
              size = line_size) +
    scale_x_continuous(name = 'Time (s)',
                       limits = c(0, 60),
                       breaks = seq(0, 100, 5),
                       expand = c(0,0)) +
    scale_y_continuous(name = 'Norm. ΔF/F0',
                       limits = c(-0.1, 1),
                       breaks = seq(-100, 100, .2),
                       expand = c(0,0)) +
    theme_classic() +
    theme(panel.grid.major = element_line(linetype = 'dotted',
                                          size = 0.1,
                                          colour = "black"),
          text = element_text(family=font_family, face="bold", size=font_size))
  
  ggsave(sprintf('norm_int_%s.png', rec.id),
         norm.int, 
         width = 28, height = 8, units = 'cm')  
}
remove(rec.id, norm.int)


##### UP INSERTION MAX BOXPLOT #####
df.fp.max <- df.fp.up.region %>%
  group_by(ID, mask_region, loc) %>%
  summarise(max_region = max(delta),
            Int_0 = mean(delta[1:8])) %>%
  group_by(ID, loc) %>%
  summarise(int_0 = mean(Int_0),
            int_r = mean(max_region)) %>%
  pivot_longer(cols = 3:4, names_to = 'val', values_to = 'delta') %>%
  unite(loc, val, col = 'sub', sep = '_') %>%
  mutate(sub = as.factor(sub)) %>%
  filter(ID != 'cell5_02_2_2022') %>%
  ungroup() %>%
  mutate(sub = ordered(sub, c('c_int_0', 'm_int_0', 'c_int_r', 'm_int_r')))

df.fp.max.stat <- df.fp.max %>%
  pairwise_wilcox_test(delta ~ sub, p.adjust.method = 'BH') %>%
  add_xy_position(fun = 'mean_sd', scales = 'free')

insertion.max <- ggplot(df.fp.max, aes(x = sub, y = delta)) +
  geom_boxplot(fill = 'grey', alpha = .5) +
  geom_point(aes(color = ID, group = ID), size = 3) +
  geom_line(aes(color = ID, group = ID), size = line_size) +
  stat_pvalue_manual(df.fp.max.stat, label = 'p.adj.signif',
                     hide.ns = TRUE, size = 5, family=font_family) +
  scale_x_discrete(name = NULL,
                   labels = c('TGN F0','PM F0','TGN Fmax', 'PM Fmax')) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.01, 0.75),
                     breaks = seq(-100, 100, .1),
                     expand = c(0,0)) +
  scale_color_startrek() +
  scale_fill_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('fp_region_max.png', insertion.max, 
       width = 10, height = 15, units = 'cm')
remove(df.fp.max, df.fp.max.stat, insertion.max)


##### UP TRANSLOCATION MAX BOXPLOT #####
# NOT WORKONG!!!!!!
df.trans.max <- df.fp.up.region %>%
  group_by(ID, mask_region, loc) %>%
  summarise(max_region = max(rel),
            int_0 = mean(rel[1:8])) %>%
  group_by(ID, loc) %>%
  summarise(int_0 = mean(int_0),
            int_r = mean(max_region)) %>%
  pivot_longer(cols = 3:4, names_to = 'val', values_to = 'rel') %>%
  unite(loc, val, col = 'sub', sep = '_') %>%
  mutate(sub = as.factor(sub)) %>%
  ungroup() %>%
  mutate(sub = ordered(sub, c('c_int_0', 'm_int_0', 'c_int_r', 'm_int_r')))

ggplot(df.trans.max, aes(x = sub, y = rel)) +
  geom_boxplot(fill = 'grey', alpha = .5) +
  geom_point(aes(color = ID, group = ID), size = 3) +
  geom_line(aes(color = ID, group = ID), size = line_size)

##### FP vs .Ca #####
tail.frame <- 22
bad.cells <- c('cell5_02_2_2022', 'cell4_02_2_2022')

df.tail <- df.fp.up.region %>%
  filter(frame >= tail.frame, loc == 'c') %>%
  select(ID, ch, frame, delta) %>%
  group_by(ID, ch, frame) %>%
  summarise(delta = mean(delta)) %>%
  ungroup() %>%
  rbind(., df.ca %>%
          select(ID, ch, frame, delta) %>%
          filter(frame >= tail.frame)) %>%
  ungroup()

df.dd <- df.tail %>%
  filter(delta >= 0.05) %>%
  pivot_wider(names_from = ch, values_from = delta) %>%
  drop_na() %>%
  select(ID, fp, ca)

# data log transform
df.log <- df.dd %>%
  group_by(ID) %>%
  mutate(log_y = log10(fp/(1-fp))) %>%
  mutate(log_c = log10(ca))  %>%
  ungroup()

# linear model for ungrouped data
general.lm <- lm(log_y ~ log_c, data = df.log %>%
                   filter(!(ID %in% bad.cells)))
general.lm <- summary(general.lm)
general.lm

# multiple linear models
df.id.lm <- df.log %>%
  group_by(ID) %>%
  do(tidy(lm(log_y ~ log_c, .))) %>%
  select(ID, term, estimate, std.error) %>%
  ungroup()

# multiple lm predictions
df.hill <- df.id.lm %>%
  select(ID, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(K = '(Intercept)',
         n = 'log_c') %>%
  group_by(ID) %>%
  mutate(K = round(10^K, digits = 3),
         n = round(n, digits = 3))

# # multiple lm coefficient z-test
# compare.coeff <- function(b1,se1,b2,se2){
#   return((b1-b2)/sqrt(se1^2+se2^2))
# }
# p_value = 2*pnorm(-abs(compare.coeff(b.TBI,se.TBI,b.cef,se.cef)))
# 
# df.lm.comparison <- df.id.lm %>%
#   filter(!(ID %in% bad.cells), term == 'log_c') %>%
#   rename(n = estimate, se = std.error) %>%
#   select(ID, n, se) %>%
#   droplevels()

# df.lm.comparison
# df.lm.comparison$ID 
# interaction(df.lm.comparison$ID, df.lm.comparison$ID)

# line model plot
hill.table <- ggtexttable(df.hill,
            rows = NULL,
            theme = ttheme("light"))

lm.mod.label <- sprintf('R^2=%.3f \n K=%.3f \n n=%.3f+/-%.3f',
                        general.lm$adj.r.squared,       # R^2
                        10^general.lm$coefficients[1],  # Kd
                        general.lm$coefficients[2],     # n
                        general.lm$coefficients[4])     # n SD

lm.dd <- ggplot(df.dd) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_smooth(data = df.log %>% filter(!(ID %in% bad.cells)),
              mapping = aes(x = log_c, y = log_y),
              method = 'lm', alpha = .2, color = 'grey') +
  geom_point(mapping = aes(y = log10(fp/(1-fp)),
                           x = log10(ca),
                           color = ID,
                           fill = ID)) +
  geom_smooth(mapping = aes(y = log10(fp/(1-fp)),
                            x = log10(ca),
                            color = ID,
                            fill = ID),
              method = 'lm', se = FALSE) +
  geom_text(x = -0.43, y = -0.15,
            family = font_family, face = "bold",
            label = lm.mod.label) +
  labs(x = 'log(L)',
       y = 'log(Y/(1-Y))') +
  scale_color_startrek() +
  scale_fill_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

hill.fit <- ggarrange(lm.dd, hill.table,
          labels = c("A", "B"),
          ncol = 2, nrow = 1, widths = c(3, 1))

ggsave('dd_lm.png', hill.fit, 
       width = 30, height = 10, units = 'cm')

# native dd plot
native.dd <- ggplot() +
  geom_smooth(data = df.dd %>% filter(fp>0), mapping = aes(x = ca, y = fp),
              colour = 'grey', size = line_size + 0.5, alpha = .2,
              method='glm', method.args = list(family=binomial)) +
  geom_line(data = df.dd, mapping = aes(x = ca, y = fp, color = ID),
            size = line_size) +
  geom_point(data = df.dd, mapping = aes(x = ca, y = fp, color = ID)) +
  xlab('Fluo-4 ΔF/F0') +
  ylab('HPCA ΔF/F0') +
  scale_color_startrek() +
  scale_fill_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.25, 1),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('dd_native.png', native.dd, 
       width = 20, height = 15, units = 'cm')

remove(native.dd, lm.dd, df.tail, hill.table, df.hill, hill.fit, df.id.lm, df.log, general.lm)
