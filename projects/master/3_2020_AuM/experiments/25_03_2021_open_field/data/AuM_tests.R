# AuM behavior tests results analysis.
# Copyright © 2021 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(ggplot2)
require(ggsci)
require(ggpubr)

setwd('/home/astria/bio/note/projects/3_2020_AuM/experiments/25_03_2021_open_field/data')

#### OPEN FIELD ####
df.of <- read.csv('of.csv') %>%
         rename(Distance = Distance..m) %>%
         subset(select = c(Group, Day, Distance)) %>%
         mutate(Group = as.factor(Group)) %>%
         mutate(Group = recode_factor(Group, '1' = 'Sham+AuM', '2' = 'SNI+AuM'))

df.of.stat <- df.of %>%
           group_by(Day, Group) %>%
           summarise_all(list(mean, sd), na.rm = TRUE) %>%
           rename(mean = fn1, sd = fn2)

ggplot(df.pf.stat) +
  geom_point(aes(x = Day,
                 y = mean,
                 colour = Group)) +
  geom_line(aes(x = Day,
                y = mean,
                colour = Group),
            size = 1) +
  geom_errorbar(aes(x = Day,
                    ymin = mean-sd, ymax = mean+sd,
                    colour = Group),
                size = 0.5,
                width = 0.7) +
  labs(title = 'Open field test',
       y = 'Distance (m)') +
  scale_color_jco() +
  theme_minimal()

#### VON FREY ####
df.vf <- read.csv('von_frey.csv') %>%
         subset(select = c(Group, Day, Paw, Threshold)) %>%
         mutate(Group = recode_factor(Group, '1' = 'Sham+AuM', '2' = 'SNI+AuM')) %>%
         mutate(Threshold = as.numeric(Threshold)) %>%
         drop_na()

df.vf.stat <- df.vf %>%
              group_by(Group, Day, Paw) %>%
              summarise_all(list(mean, sd)) %>%
              rename(mean = fn1, sd = fn2)

sham.paw <- ggplot(filter(df.vf.stat, Group == 'Sham+AuM')) +
  geom_vline(xintercept = 5) +
  geom_point(aes(x = Day,
                 y = mean,
                 colour = Paw)) +
  geom_line(aes(x = Day,
                y = mean,
                colour = Paw),
            size = 1) +
  geom_errorbar(aes(x = Day,
                    ymin = mean-sd, ymax = mean+sd,
                    colour = Paw),
                size = 0.5,
                width = 0.7) +
  scale_x_continuous(breaks = seq(-100, 100, 2)) +
  labs(title = 'Sham',
       y = 'Threshold (g)') +
  scale_color_jco() +
  theme_minimal()

sni.paw <- ggplot(filter(df.vf.stat, Group == 'SNI+AuM')) +
  geom_vline(xintercept = 5) +
  geom_point(aes(x = Day,
                 y = mean,
                 colour = Paw)) +
  geom_line(aes(x = Day,
                y = mean,
                colour = Paw),
            size = 1) +
  geom_errorbar(aes(x = Day,
                    ymin = mean-sd, ymax = mean+sd,
                    colour = Paw),
                size = 0.5,
                width = 0.7) +
  scale_x_continuous(breaks = seq(-100, 100, 2)) +
  labs(title = 'SNI',
       y = 'Threshold (g)') +
  scale_color_jco() +
  theme_minimal()

ggarrange(sham.paw, sni.paw, ncol = 1, nrow = 2)


ggplot(filter(df.vf.stat, Paw == 'R' & Day != -7)) +
  geom_vline(xintercept = 5) +
  geom_point(aes(x = Day,
                 y = mean,
                 colour = Group)) +
  geom_line(aes(x = Day,
                y = mean,
                colour = Group),
            size = 1) +
  geom_errorbar(aes(x = Day,
                    ymin = mean-sd, ymax = mean+sd,
                    colour = Group),
                size = .5,
                width = .7) +
  scale_x_continuous(breaks = seq(-100, 100, 2)) +
  labs(title = 'Operated paw',
       y = 'Threshold (g)') +
  scale_color_jco() +
  theme_minimal(text=element_text(family="Arial"))
