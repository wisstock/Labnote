# HPCA-TagRFP + Fluo-4 profiles data analysis, img for FENS forum 2022
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

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/astria/bio/institute/my_pub/conf/2022_FENS/HPCA/poster/scripts')

##### GLOBAL PLOT OPTIONS #####
font_family <- "ubuntu mono"
font_size <- 30

line_size <- 1
line_size_light <- 1.5

ca.blue <- c('#0000FF', '#6767FF', '#00AEFF', '#366BFF')
fp.red <- c('#FF0000', '#B30018', '#FF2323', '#DE3800')


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
  filter(mask == 'up')
df.fp.cyto.region <- df.fp.up.region %>%
  data.frame()

for (cell.id in names(membrane.region)) {
  df.fp.up.region <- df.fp.up.region %>%
    filter(!(ID == cell.id & !(mask_region %in% membrane.region[[cell.id]]))) %>%
    mutate(loc = 'm')
  
  df.fp.cyto.region <- df.fp.cyto.region %>%
    filter(!(ID == cell.id & mask_region %in% membrane.region[[cell.id]])) %>%
    mutate(loc = 'c')
}

df.fp.up <- rbind(df.fp.up.region, df.fp.cyto.region) %>%
  mutate(loc = as.factor(loc))
remove(cell.id, df.fp.cyto.region, df.fp.up.region)


##### UP INSERTION MAX BOXPLOT #####
df.fp.max <- df.fp.up %>%
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
  mutate(sub = ordered(sub, c('m_int_0', 'm_int_r', 'c_int_0', 'c_int_r'))) %>%
  mutate(loc = recode(sub, 'c_int_0' = 'c',
                           'c_int_r' = 'c',
                            'm_int_0' = 'm',
                            'm_int_r' = 'm')) %>%
  unite(ID, loc, col = 'loc_id', sep = '_', remove = FALSE)

df.fp.max.stat <- df.fp.max %>%
  pairwise_wilcox_test(delta ~ sub, p.adjust.method = 'BH') %>%
  add_xy_position(fun = 'mean_sd', scales = 'free')

insertion.max <- ggplot(df.fp.max, aes(x = sub, y = delta)) +
  geom_point(aes(group = sub, color = loc),
             size = 3) +
  geom_line(aes(group = loc_id, color = loc, linetype = loc),
            size = line_size) +
  geom_boxplot(aes(group = sub, fill = loc),
               alpha = .5) +
  stat_pvalue_manual(df.fp.max.stat, label = 'p.adj.signif',
                     hide.ns = TRUE, size = 5, family=font_family) +
  scale_fill_manual(name = NULL,
                    values = fp.red) +
  scale_color_manual(name = NULL,
                     values = fp.red) +
  scale_x_discrete(name = NULL,
                   labels = c('PM Fmax', 'PM Fmax', 'IMs F0','IMs Fmax')) +
  scale_y_continuous(name = 'HPCA-tagRFP ΔF/F0',
                     limits = c(-0.01, 0.75),
                     breaks = seq(-100, 100, .1),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = 'none',
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('ins_max.png', insertion.max, 
       width = 29.8, height = 11.1, units = 'cm', dpi = 300)
# remove(df.fp.max, df.fp.max.stat, insertion.max)




##### REPRESENTATIVE CELL FILTERING #####
rep.cell.ID <- 'cell7_02_2_2022'
df.fp <- df.fp %>% filter(ID == rep.cell.ID)
df.fp.up <- df.fp.up %>% filter(ID == rep.cell.ID)

df.ca <- df.ca %>%
  filter(ID == rep.cell.ID)
ca.delta.der <- append(df.ca$delta[-1] - rev(rev(df.ca$delta)[-1]), 0, after = 0)
df.ca <- df.ca %>%
  mutate(der = ca.delta.der)
remove(ca.delta.der)



##### CTRL CA & FP PROFILES #####
stim.i <- c(12, 23, 33)
frame.i <- c(6, 20, 46)

ca.prof <- ggplot() +
  geom_line(data = df.ca,
            aes(x = time, y = delta, color='Ca2+ dye'),
            size = line_size,
            color = 'blue') +
  geom_point(data = df.ca,
             aes(x = time, y = delta),
             show.legend=FALSE,
             color = 'blue') +
  geom_line(data = df.fp %>% filter(mask == 'master'),
            aes(x = time, y = delta*5),
            size = line_size,
            color = 'red') +
  geom_point(data = df.fp %>% filter(mask == 'master'),
             aes(x = time, y = delta*5),
             show.legend=FALSE,
             color = 'red') +
  geom_segment(aes(x = stim.i, xend = stim.i,
                   y = c(1.9, 1.9, 1.9), yend = c(1.5, 1.5, 1.5)),
               arrow = arrow(length = unit(0.2, "cm")),
               size = line_size) +
  geom_segment(aes(x = frame.i, xend = frame.i,
                   y = c(-0.35, -0.35, -0.35), yend = c(-0.2, -0.2, -0.2)),
               arrow = arrow(length = unit(0.3, "cm")),
               size = line_size,
               color = '#808080') +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(0, 60),
                     breaks = seq(0, 100, 5),
                     expand = c(0,0)) +
  scale_y_continuous(name = "Fluo-4 ΔF/F0",
                     limits = c(-0.5, 2),
                     breaks = seq(-100, 100, 0.5),
                     sec.axis = sec_axis(trans=~./5, name="HPCA ΔF/F0")) +
  scale_color_discrete(name="") +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family, face="bold", size=font_size),
        axis.line.y = element_line(color = "blue"),
        axis.ticks.y = element_line(color = "blue"),
        axis.line.y.right = element_line(color = "red"),
        axis.ticks.y.right = element_line(color = "red"))

ggsave('ca_profile.png', ca.prof, 
       width = 32.5, height = 10, units = 'cm', dpi = 300)
remove(ca.prof, stim.i, frame.i)


##### NORM PROFILES #####
df.fp.up.norm <- df.fp.up %>%
  select(ID, time, delta, loc) %>% 
  group_by(ID, time, loc) %>%
  summarise(loc_mean = mean(delta), loc_sd = IQR(delta)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(loc_mean = loc_mean / max(loc_mean)) %>%
  mutate(loc = ordered(loc, c('m', 'c'))) %>%
  mutate(loc = recode(loc, 'm' = 'PM', 'c' = 'IMs'))
  # mutate(loc = revalue(cyl, c('PM' = 'm', 'IMs' = 'c')))

df.ca.norm <- df.ca %>%
  select(ID, time, delta) %>%
  group_by(ID) %>%
  mutate(delta_norm = delta / max(delta))

norm.prof <- ggplot() +
  geom_line(data = df.ca.norm,
            aes(x = time, y = delta_norm),
            color = 'blue',
            size = line_size) +
  geom_point(data = df.ca.norm,
             aes(x = time, y = delta_norm),
             color = 'blue') +
  geom_line(data = df.fp.up.norm,
            aes(x = time, y = loc_mean, linetype = loc, color = loc),
            size = line_size) +
  geom_point(data = df.fp.up.norm,
             aes(x = time, y = loc_mean, color = loc)) +
  geom_ribbon(data = df.fp.up.norm,
              aes(x = time,
                  ymin = loc_mean-loc_sd,
                  ymax = loc_mean+loc_sd,
                  linetype = loc,
                  color = loc,
                  fill = loc),
              alpha=.15,
              size=0.5) +
  scale_fill_manual(name = NULL,
                    values = fp.red,
                    labels = c('Plasma membrane (PM)',
                               'Intracellular membranes (IMs)')) +
  scale_color_manual(name = NULL,
                     values = fp.red,
                     labels = c('Plasma membrane (PM)',
                                'Intracellular membranes (IMs)')) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(0, 60),
                     breaks = seq(0, 100, 5),
                     expand = c(0,0)) +
  scale_y_continuous(name = "Norm. ΔF/F0",
                     limits = c(-0.2, 1.1),
                     breaks = seq(-100, 100, 0.25)) +
  guides(linetype = FALSE) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.43, 1),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('norm_fp_profile.png', norm.prof, 
       width = 36.3, height = 10.5, units = 'cm', dpi = 300)
# remove(df.fp.up.norm, df.ca.norm, norm.prof)
