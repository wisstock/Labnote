# Flyo-4 data analysis, AM form loading protocols without (1) and with (2) PA
# Copyright © 2020 Borys Olifirov

require(ggplot2)
require(ggpubr)
require(dplyr)
require(magrittr)


setwd('/home/astria/Bio/Note/Projects/HPCA/09.2020/09_21_2020')

df <- read.csv('results.csv')

ggplot(df, aes(x=time, y=int, colour = as.factor(cell))) +
  geom_line(na.rm = TRUE, size = 1)


df <- subset(df, cell != 6, select = c('PA', 'time', 'int'))
df.sum <- df %>% group_by(time, PA) %>% summarise_all(list(mean, sd))

ggplot(df.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(PA),
                  colour=as.factor(PA)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(PA)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(PA)),
             size=1.5) +
  geom_vline(xintercept = 23, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 180),
                     breaks = seq(0, 180, 20)) +
  labs(y ='Intensity', x = 'Time (s)') +
  scale_fill_discrete(name = 'Protocol', labels = c('no PA', 'PA')) +
  scale_colour_discrete(name = 'Protocol', labels = c('no PA', 'PA')) +
  theme_minimal(base_size = 18,
                base_family = 'oswald')


  ggpubr::color_palette("jco") +
  ggpubr::fill_palette("jco")
