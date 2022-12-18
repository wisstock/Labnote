# HPCA-TagRFP + Fluo-4 data analysis
# Copyright © 2022 Borys Olifirov
#
# px-wise intensity and dist data frame structure:
# 'ID' - recording ID
# 'stim' - stimulation number
# 'mask_element' - up mask element number
# 'd' - px distance from nucleus region border
# 'int' - px intensity
# 'delta' - px delta F

require(dplyr)
require(tidyr)
require(purrr)

require(rstatix)  # -

require(ggplot2)
require(ggpubr)  # -
require(cowplot)
require(ggsci)

setwd('/home/wisstock/Bio/note/projects/PhD/6_2021_PIP2_HPCA_HEK/exp/2022_12_15_dist_and_FENS_pic/scripts')

##### GLOBAL PLOT OPTIONS #####
font_family <- "ubuntu mono"
font_size <- 30

line_size <- 1
line_size_light <- 1.5


##### DATA PREPROCESSING #####
# up mask region selection
membrane.region <- list('cell4_02_2_2022' = c(1, 2, 4, 5),
                        'cell5_02_2_2022' = c(1),
                        'cell6_02_2_2022' = c(2),
                        'cell7_02_2_2022' = c(1, 3, 5))


# px-wise up mask intensity section
df.px <- read.csv('px_05_19_2022.csv') %>%
  filter(ID !='cell4_02_2_2022') %>%
  mutate(ID = as.factor(ID),
         mask_element = as.factor(mask_element),
         stim = as.factor(stim)) %>%
  mutate(d = d*0.138)

df.px.cyto <- df.px %>% data.frame()
for (cell.id in names(membrane.region)) {
  df.px <- df.px %>%
    filter(!(ID == cell.id & !(mask_element %in% membrane.region[[cell.id]]))) %>%
    mutate(loc = 'm')
  
  df.px.cyto <- df.px.cyto %>%
    filter(!(ID == cell.id & mask_element %in% membrane.region[[cell.id]])) %>%
    mutate(loc = 'c')
}
df.px <- rbind(df.px, df.px.cyto) %>%
  mutate(loc = as.factor(loc))
remove(df.px.cyto, membrane.region, cell.id)

df.stim.delta <- df.px %>%
  select(-int) %>%
  group_by(ID, stim) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = stim, values_from = delta) %>%
  select(-row) %>%
  rename(stim_1 = '1', stim_2 = '2', stim_3 = '3') %>%
  mutate(ddelta = stim_3 - stim_1,
         d = round(d, digits = 2),
         d_fac = as.factor(d)) %>%
  ungroup()

##### STIM VS D PX-WISE#####
# delta.vs.d <-
ggplot(df.px %>% filter(loc == 'm', stim != '2'),
       aes(x = d, y = delta, color = stim)) +
  geom_point(alpha = .2, size = 1) +
  geom_smooth(aes(linetype = stim)) +
  facet_wrap(~ ID) +
  scale_x_continuous(name = 'Distance (um)',
                     limits = c(0, 3),
                     breaks = seq(0, 100, 0.5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(-0.5, 2),
                     breaks = seq(0, 100, 0.25),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family,
                            size=font_size))


##### REPRESENTATIVE CELL #####
rep.cell.ID <- 'cell7_02_2_2022'
df.px.rep <- df.px %>% filter(ID == rep.cell.ID)

fp.red <- c('#FF0000', '#FF8400', '#FF2323', '#DE3800')

delta.d <- ggplot(df.px.rep %>% filter(loc == 'm', stim != 2),
       aes(x = d, y = delta, color = stim, fill = stim)) +
  geom_point(aes(shape = stim), alpha = .25, size = 2) +
  geom_smooth() +
  scale_fill_manual(name = NULL,
                    values = fp.red,
                    labels = c('After 1st stimul','After 3d stimul')) +
  scale_color_manual(name = NULL,
                     values = fp.red,
                     labels = c('After 1st stimul','After 3d stimul')) +
  scale_x_continuous(name = 'Distance (um)',
                     limits = c(0, 5),
                     breaks = seq(0, 100, 0.5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'HPCA-tagRFP ΔF/F0',
                     limits = c(-0.1, 1.25),
                     breaks = seq(0, 100, 0.25),
                     expand = c(0,0)) +
  guides(shape = FALSE) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.7, 1),
        legend.justification = c("left", "top"),
        text = element_text(family=font_family,
                            face="bold",
                            size=font_size))

ggsave('delta_vs_d.png',
       delta.d,
       width = 37.3, height = 22.8, units = 'cm', dpi = 300)
# remove(delta.vs.d)


##### STIM VS. DELTA MEAN #####
df.delta.mean <- df.px %>%
  filter(loc == 'm') %>%
  select(ID, stim, d, delta) %>%
  mutate(d = round(d, digits = 2),
         d_f = as.factor(d)) %>%
  group_by(ID, stim, d_f) %>%
  mutate(stim_mean = mean(delta),
         stim_sd = sd(delta)) %>%
  ungroup()


ggplot(df.delta.mean) +
  geom_line(aes(x = delta, y = d, color = stim),
            size = line_size) +
  facet_wrap(~ ID, scales = "free")
geom_point(data = df.fp.up.norm,
           aes(x = time, y = loc_mean),
           color = 'red') +
  geom_ribbon(data = df.fp.up.norm,
              aes(x = time,
                  ymin = loc_mean-loc_sd,
                  ymax = loc_mean+loc_sd,
                  linetype = loc),
              color = 'red',
              fill = 'red',
              alpha=.15,
              size=0.5) +
  
  



##### DDELTA VS D PX-WISE#####
d.delta.vs.d <- ggplot(df.stim.delta %>% filter(loc == 'm'),
                       aes(x = d, y = ddelta,
                           color = mask_element,
                           shape = mask_element)) +
  geom_point(alpha = .5, size = 2) +
  facet_wrap(~ ID, scales = "free") +
  labs(x = 'Distance, um', y = 'Δ(ΔF/F0)') +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family, face="bold",
                            size=font_size),
        plot.background = element_rect(fill = "black"),
        axis.title.x = element_text(colour = "white"),
        axis.title.y = element_text(colour = "white"),
        axis.text.x = element_text(face="bold", color='white'),
        axis.text.y = element_text(face="bold", color="white"))

  
ggsave('d_delta_vs_d.png',
       d.delta.vs.d,
       width = 28, height = 18, units = 'cm')
remove(d.delta.vs.d)


##### 1ST VS 3D STIM PX-WISE#####
# https://stats.stackexchange.com/questions/539011/suitable-approach-to-cluster-histogram-like-dataset-using-hdbscan-implementation
# total delta vs. delta    
ggplot(df.stim.delta %>% filter(loc == 'm'),
       aes(x = stim_1, y = stim_3,
           color = d)) +
  geom_point(alpha = .6, size = 2) +
  facet_wrap(~ ID, scales = "free") +
  scale_colour_gradient(low = 'blue', high = 'yellow')

  labs(x = '1st stim. ΔF/F0', y = '3d stim. ΔF/F0') +
  scale_color_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family, face="bold", size=font_size))
# ggsave('delta_vs_delta_total.png',
#        total.delta.scatter, 
#        width = 22, height = 15, units = 'cm')
# remove(total.delta.scatter)

ggplot(df.px.delta %>% filter(ID == 'cell7_02_2_2022'),
       aes(x = stim_1, y = stim_3, color = d, shape = mask_element)) +
  geom_point(alpha = .5, size = 2) +
  labs(x = '1st stim. ΔF/F0', y = '3d stim. ΔF/F0') +
  scale_color_startrek() +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        text = element_text(family=font_family, face="bold", size=font_size))
#        delta.cell, 
#        width = 22, height = 15, units = 'cm')
# remove(delta.cell)


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
