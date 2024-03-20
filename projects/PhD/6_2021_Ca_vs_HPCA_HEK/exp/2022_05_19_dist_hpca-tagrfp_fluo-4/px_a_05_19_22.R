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

# require(rstatix)  # -

require(ggplot2)
# require(ggpubr)  # -
require(cowplot)
require(ggsci)

setwd('/home/astria/bio/note/projects/PhD/6_2021_PIP2_HPCA_HEK/exp/2022_05_19_dist_hpca-tagrfp_fluo-4')

##### GLOBAL PLOT OPTIONS #####
font_family <- "ubuntu mono"
font_size <- 15

line_size <- 1
line_size_light <- 1.5


##### DATA PREPROCESSING #####
# up mask region selection
membrane.region <- list('cell2_02_2_2022' = c(1, 2),
                        'cell4_02_2_2022' = c(1, 2, 4, 5),
                        'cell5_02_2_2022' = c(1),
                        'cell6_02_2_2022' = c(2),
                        'cell7_02_2_2022' = c(1, 3, 5))


# px-wise up mask intensity section
df.px <- read.csv('px_05_19_2022.csv') %>%
  filter(ID !='cell4_02_2_2022') %>%
  mutate(ID = as.factor(ID),
         mask_element = as.factor(mask_element),
         stim = as.factor(stim))

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
  mutate(ddelta = stim_3 - stim_1) %>%
  ungroup()


##### STIM VS D #####
delta.vs.d <- ggplot(df.px %>% filter(loc == 'm', stim != '2'),
       aes(x = d, y = delta, color = stim, shape = mask_element)) +
  geom_point(alpha = .5, size = 2) +
  facet_wrap(~ ID, scales = "free")
ggsave('delta_vs_d.png',
       delta.vs.d,
       width = 28, height = 18, units = 'cm')
remove(delta.vs.d)


##### DDELTA VS D #####
d.delta.vs.d <- ggplot(df.stim.delta %>% filter(loc == 'm'),
                       aes(x = d, y = ddelta,
                           color = mask_element,
                           shape = mask_element)) +
  geom_point(alpha = .5, size = 2) +
  facet_wrap(~ ID, scales = "free")
ggsave('d_delta_vs_d.png',
       d.delta.vs.d,
       width = 28, height = 18, units = 'cm')
remove(d.delta.vs.d)


##### 1ST VS 3D STIM #####
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
