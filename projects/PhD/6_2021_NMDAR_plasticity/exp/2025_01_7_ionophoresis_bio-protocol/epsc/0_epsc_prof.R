# NMDA ionophoresis, EPSCs record plot for bio-protocol
# Copyright © 2025 Borys Olifirov

require(dplyr)
require(tidyr)

require(ggplot2)


setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_01_7_ionophoresis_bio-protocol/epsc')

font.size <- 17
font.fam <- 'Arial'
box.alpha <- 0.6

base.01 <- read.csv('traces/base_01.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))
base.02 <- read.csv('traces/base_02.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))

wo.01 <- read.csv('traces/wo_01.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))
wo.02 <- read.csv('traces/wo_01.csv', dec = ",") %>%
  slice(which(row_number() %% 2 == 1))


df.rep <- rbind(base.01 %>% slice(1:31000) %>% mutate(time = 'base'),
                wo.02 %>% slice(1:31000) %>% mutate(time = 'wo')) %>%
  mutate(time = as.factor(time)) %>%
  group_by(time) %>%
  mutate(index = seq(1, 31000),
         pa = if_else(time == 'base', pa + 80, pa - 80))
  

##### REPRESENTATIVE EVENTS #####
ggplot(data = df.rep %>% filter(time == 'wo', index %in% seq(16550, 16610)),
       aes(x = index, y = pa, group = time)) +
  geom_line() +
  theme_classic()

# base
df.base.rep <- df.rep %>%
  filter(time == 'base') %>%
  mutate(event = as.factor(case_when(between(index, 89, 149) ~ '1',
                                     between(index, 4282, 4342) ~ '2',
                                     between(index, 16085, 16145) ~ '3',
                                     between(index, 11372, 11432) ~ '4',
                                     between(index, 25920, 25980) ~ '5',
                                     .default = '-'))) %>%
  filter(event != '-') %>%
  droplevels() %>%
  group_by(event) %>%
  mutate(e_index = seq(1,61)) %>%
  ungroup()

event_base <- ggplot(data = df.base.rep %>% mutate(pa = pa-160) %>% filter(e_index %in% seq(15, 53)),
       aes(x = e_index, y = pa)) +
  geom_line(aes(group = event), size = 0.15) +
  stat_summary(fun = median,
               geom = 'line', size = 2, color = 'coral') +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25, fill = 'coral') +
  theme_void() +
  annotate(geom = 'rect',  # 10 pA scale bar
           xmin = 15, xmax = 16,
           ymin = -35, ymax = -25,
           fill = 'black') +
  annotate(geom = 'rect',  # 5 ms scale bar
           xmin = 15, xmax = 23.3,
           ymin = -35, ymax = -34.25,
           fill = 'black') +
  scale_y_continuous(limits = c(-35, 15))
  

event_base
save_plot('0_pic_base_epsc.png', event_base, base_width = 4, base_height = 7, dpi = 300)


# wo
df.wo.rep <- df.rep %>%
  filter(time == 'wo') %>%
  mutate(event = as.factor(case_when(between(index, 1543, 1603) ~ '1',
                                     between(index, 2195, 2255) ~ '2',
                                     between(index, 5140, 5200) ~ '3',
                                     between(index, 10293, 10353) ~ '4',
                                     between(index, 16550, 16610) ~ '5',
                                     .default = '-'))) %>%
  filter(event != '-') %>%
  droplevels() %>%
  group_by(event) %>%
  mutate(e_index = seq(1,61)) %>%
  ungroup()

event_wo <- ggplot(data = df.wo.rep %>% mutate(pa = pa-22) %>% filter(e_index %in% seq(17, 55)),
       aes(x = e_index, y = pa)) +
  geom_line(aes(group = event), size = 0.15) +
  stat_summary(fun = median,
               geom = 'line', size = 2, color = 'deeppink2') +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25, fill = 'deeppink2') +
  theme_void() +
  scale_y_continuous(limits = c(-35, 15))

event_wo
save_plot('0_pic_wo_epsc.png', event_wo, base_width = 4, base_height = 7, dpi = 300)

# df.base.rep <- base.01 %>%
#   mutate(index = seq(1,31706)) %>%
#   mutate(event = as.factor(case_when(between(index, 90*2, 150*2) ~ '1',
#                                      between(index, 4280*2, 4340*2) ~ '2',
#                                      between(index, 16080*2, 16140*2) ~ '3',
#                                      between(index, 11370*2, 11430*2) ~ '4',
#                                      between(index, 25920*2, 25980*2) ~ '5',
#                                      .default = '-'))) %>%
#   filter(event != '-') %>%
#   droplevels() %>%
#   group_by(event) %>%
#   mutate(e_index = seq(1,121)) %>%
#   ungroup()


##### INTERVAL PROF #####
prof_epsc <- ggplot(data = df.rep,
       aes(x = index, y = pa, group = time)) +
  geom_line() +
  annotate(geom = 'rect',  # 20 pA scale bar
           xmin = 0, xmax = 280,
           ymin = -44, ymax = -24,
           fill = 'black') +
  annotate(geom = 'rect',  # 10 s scale bar
           xmin = 0, xmax = 3300,
           ymin = -44, ymax = -40,
           fill = 'black') +
  theme_void() +
  theme(legend.position = 'none')

prof_epsc
save_plot('0_pic_prof_epsc.png', prof_epsc, base_width = 7, base_height = 4.4, dpi = 300)
remove(prof_epsc)