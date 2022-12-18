# Flyo-4 data analysis
# Copyright © 2020 Borys Olifirov

require(ggplot2)
require(ggpubr)
require(dplyr)
require(magrittr)

setwd('/home/astria/Bio/Note/Projects/HEPCA/07.2020/07_26_2020')

df.100 <- read.csv('100_us.csv')
df.exp <- read.csv('exposures.csv')
df.blc <- read.csv('bleaching.csv')


##### EXPOSURE #####
df.exp <- subset(df.exp, select = c(exp, time, int))
df.exp.sum <- df.exp %>% group_by(time, exp) %>% summarise_all(list(mean, sd))

ggplot(df.exp.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(exp),
                  colour=as.factor(exp)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(exp)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(exp)),
             size=1.5) +
  geom_vline(xintercept = 25, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 180),
                     breaks = seq(0, 180, 20)) +
  labs(y ='Intensity', x = 'Time (s)') +
  guides(fill=guide_legend(title='Exposure (us/px)'),
         colour=guide_legend(title='Exposure (us/px)')) +
  theme_minimal(base_size = 18,
                base_family = 'oswald') +
  ggpubr::color_palette("jco") +
  ggpubr::fill_palette("jco")

# scale_y_continuous(limits = c(200, 900),
#                    labels = function(x) paste0(round(x/max(df.exp.sum$fn1)*100,0), "%"))

##### CYCLES ####
df.100 <- subset(df.100, select = c(cycl, time, int))
df.100.sum <- df.100 %>% group_by(time, cycl) %>% summarise_all(list(mean, sd))

ggplot(df.100.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(cycl),
                  colour=as.factor(cycl)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(cycl)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(cycl)),
             size=1.5) +
  geom_vline(xintercept = 25, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 180),
                     breaks = seq(0, 180, 20)) +
  labs(y ='Intensity', x = 'Time (s)') +
  guides(fill=guide_legend(title='Exposure cycles'),
         colour=guide_legend(title='Exposure cycles')) +
  theme_minimal(base_size = 18,
                base_family = 'oswald') +
  ggpubr::color_palette("jco") +
  ggpubr::fill_palette("jco")

##### BLEACHING #####
df.blc <- subset(df.blc, select = c(cycl, time, int))
df.blc.sum <- df.blc %>% group_by(time, cycl) %>% summarise_all(list(mean, sd))

ggplot(df.blc.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(cycl),
                  colour=as.factor(cycl)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(cycl)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(cycl)),
             size=1.5) +
  geom_vline(xintercept = 25, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 180),
                     breaks = seq(0, 180, 20)) +
  labs(y ='Intensity', x = 'Time (s)') +
  guides(fill=guide_legend(title='Iteration'),
         colour=guide_legend(title='Iteration')) +
  theme_minimal(base_size = 18,
                base_family = 'oswald') +
  ggpubr::color_palette("jco") +
  ggpubr::fill_palette("jco")
