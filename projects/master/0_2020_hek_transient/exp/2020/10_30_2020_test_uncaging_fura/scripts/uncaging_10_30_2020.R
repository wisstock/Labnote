# Flyo-4 and Fura Red data analysis, uncaging with different stimulation power
# Copyright © 2020 Borys Olifirov

require(ggplot2)
require(dplyr)
require(magrittr)

setwd('/home/astria/bio/note/projects/hek_transient/experiments/test_uncaging_fura/scripts')

df.raw <- read.csv('results.csv')



df <- subset(df.raw, select = c(power, time, int))
df.sum <- df %>% group_by(time, power) %>% summarise_all(list(mean, sd), na.rm = TRUE)
df.sum[df.sum < 0] <- 0

ggplot(df.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(power),
                  colour=as.factor(power)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(power)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(power)),
             size=1.5) +
  geom_vline(xintercept = 4, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  labs(y ='ΔF/F0', x = 'Time (s)', title = 'Fluo-4 transient') +
  guides(fill=guide_legend(title='405 nm power (%)'),
         colour=guide_legend(title='405 nm power (%)')) +
  theme_minimal(base_size = 18)


df_dd <- subset(df.sum, time == 4, select = c(power, fn1, fn2))

ggplot(df_dd, aes(x=power)) +
  geom_ribbon(aes(ymin=fn1-(fn2/fn1),
                  ymax=fn1+(fn2/fn1)),
              alpha=.5) +
  geom_line(aes(y=fn1)) +
  geom_point(aes(y=fn1)) +
  labs(y ='ΔF/F0', x = '405 power (%)',
       title = 'Fluo-4 exposure-effect',
       subtitle = 'Gray area - corrected sd (sd / F/F0)') +
  theme_minimal(base_size = 18)

