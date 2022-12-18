# Flyo-4 data analysis, uncaging with different stimulation power and stimulation area
# Copyright © 2020 Borys Olifirov

require(ggplot2)
# require(ggpubr)
require(dplyr)
require(magrittr)

setwd('/home/astria/Bio/Note/projects/HPCA_transient/experiments/test_uncaging_areas_old_HEK/scripts')

df.80 <- read.csv('results_80.csv')
df.1_2 <- read.csv('results_1-2.csv')
df.1_5 <- read.csv('results_1-5.csv')


##### 405 nm 80% #####
df.80 <- subset(df.80, select = c(area, time, int))
df.80.sum <- df.80 %>% group_by(time, area) %>% summarise_all(list(mean, sd), na.rm = TRUE)

ggplot(df.80.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(area),
                  colour=as.factor(area)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(area)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(area)),
             size=1.5) +
  geom_vline(xintercept = 6, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  labs(y ='Intensity', x = 'Time (s)', title = '405 nm 80%') +
  guides(fill=guide_legend(title='Stimulation area'),
         colour=guide_legend(title='Stimulation area')) +
  theme_minimal(base_size = 18,
                base_family = 'oswald')


##### AREA 1/2 ####
df.1_2 <- subset(df.1_2, select = c(area, time, int))
df.1_2.sum <- df.1_2 %>% group_by(time, area) %>% summarise_all(list(mean, sd))

ggplot(df.1_2.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(area),
                  colour=as.factor(area)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(area)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(area)),
             size=1.5) +
  geom_vline(xintercept = 6, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  labs(y ='Intensity', x = 'Time (s)', title = 'Area 1/2') +
  guides(fill=guide_legend(title='Laser power'),
         colour=guide_legend(title='Laser power')) +
  theme_minimal(base_size = 18,
                base_family = 'oswald')


##### AREA 1/5 #####
df.1_5 <- subset(df.1_5, select = c(area, time, int))
df.1_5.sum <- df.1_5 %>% group_by(time, area) %>% summarise_all(list(mean, sd))

ggplot(df.1_5.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(area),
                  colour=as.factor(area)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(area)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(area)),
             size=1.5) +
  geom_vline(xintercept = 6, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  labs(y ='Intensity', x = 'Time (s)', title = 'Area 1/5') +
  guides(fill=guide_legend(title='Laser power'),
         colour=guide_legend(title='Laser power')) +
  theme_minimal(base_size = 18,
                base_family = 'oswald')
