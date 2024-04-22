# mEPSC analysis
# Copyright © 2024 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)
require(dunn.test)

setwd('/home/wisstock/bio/note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_04_19_mEPSC/-40_Ca')

df.full <- bind_rows(read.csv('220512_s13.csv'),
                     read.csv('220512_s23.csv')) %>%
  mutate(Sample = as.factor(Sample),
         Group = as.factor(Group)) %>%
  mutate(Group = factor(Group, levels = c("-70ctrl", "-40ctrl", "-40post", "-70post")))


##### AMPLITUDE STAT #####
df.amp <- df.full %>%
  select(Amp, Sample, Group, TimPeak)


df.amp.norm <- df.amp %>%
  group_by(Sample) %>%
  mutate(Amed = Amp / median(Amp[Group == '-70ctrl'])) %>%
  ungroup()

ggplot() +
  geom_boxplot(data = df.amp.norm,
               aes(x = Group, y = Amed, fill = Sample,
                   group = interaction(Group, Sample))) +
  stat_summary(data = df.amp.norm,
               fun = median,
               geom = 'line',
               aes(x = Group, y = Amed, color = Sample,
                   group = Sample),
               position = position_dodge(width = 0.75)) +
  stat_summary(data = df.amp.norm,
               fun = median,
               geom = 'point',
               aes(x = Group, y = Amed, color = Sample,
                   group = Sample),
               alpha = 0.75,
               position = position_dodge(width = 0.75))


df.amp.70 <- df.amp %>%
  select(-TimPeak)
  filter((Group == '-70ctrl') | (Group == '-70post')) %>%
  droplevels()
df.amp.40 <- df.amp %>%
  select(-TimPeak) %>%
  filter((Group == '-70ctrl') | (Group == '-70post')) %>%
  select(Amp, Sample, Group) %>%
  droplevels()

ggplot() +
  geom_point(data = df.amp,
             aes(x = TimPeak, y = Amp, color = Group, alpha = .75)) +
  facet_wrap(facets = vars(Sample), nrow = nlevels(df.full$Sample), strip.position = 'right')


# box stat
df.amp.ctrl <- df.amp %>%
  group_by(Group) %>%
  wilcox_test(Amp~Sample)

df.amp.stat <- df.amp %>%
  group_by(Sample) %>%
  wilcox_test(Amp~Group) %>%
  add_significance() %>%
  add_y_position(fun = 'median_iqr', scales = 'free')


ggplot() +
  geom_boxplot(data = df.full,
               aes(x = Group, y = Amp, fill = Group))

##### IEI STAT #####
df.iei <- df.full %>%
  select(IEI, Sample, Group, TimPeak)


# prof
ggplot() +
  geom_point(data = df.iei,
             aes(x = TimPeak, y = IEI, color = Group, alpha = .75))


ggplot() +
  geom_boxplot(data = df.full,
               aes(x = Group, y = IEI, fill = Group))

##### AREA STAT #####