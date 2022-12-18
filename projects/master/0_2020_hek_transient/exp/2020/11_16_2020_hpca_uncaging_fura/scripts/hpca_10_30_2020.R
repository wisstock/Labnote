# Flyo-4 and Fura Red data analysis, uncaging with different stimulation power
# Copyright © 2020 Borys Olifirov

require(ggplot2)
require(dplyr)
require(magrittr)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/11_16_2020_hpca_uncaging_fura/scripts')

df.fluo <- read.csv('results_fluo.csv')
df.fluo <- subset(df.fluo, feature == 5)
df.hpca <- read.csv('results_hpca.csv')
df.hpca <- subset(df.hpca,feature == 'reg')
df <- rbind(df.fluo,df.hpca)

df.feature <- subset(df, select = c(feature, time, int))
df.sum <- df.feature %>% group_by(time, feature) %>% summarise_all(list(mean, sd), na.rm = TRUE)

##### ALL CURVES #####
ggplot(df.sum, aes(x=time)) +
  geom_line(aes(y=fn1, colour=as.factor(feature))) +
  geom_point(aes(y=fn1, colour=as.factor(feature))) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  colour=as.factor(feature),
                  fill=as.factor(feature)),
              alpha=.3,
              size=.2) +
  geom_vline(xintercept = 3, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  labs(y ='ΔF/F0', x = 'Time (s)') +
  scale_colour_discrete(name  = 'Registration',
                        breaks = c('5', 'reg'),
                        labels=c('Fluo-4 (405 nm 2%)', '456 nm only')) +
  scale_fill_discrete(name  = 'Registration',
                      breaks = c('5', 'reg'),
                      labels=c('Fluo-4 (405 nm 2%)', '456 nm only')) +
  theme_minimal(base_size = 18)



# df_dd <- subset(df.sum, time == 4, select = c(power, fn1, fn2))
# 
# ggplot(df_dd, aes(x=power)) +
#   geom_ribbon(aes(ymin=fn1-(fn2/fn1),
#                   ymax=fn1+(fn2/fn1)),
#               alpha=.5) +
#   geom_line(aes(y=fn1)) +
#   geom_point(aes(y=fn1)) +
#   labs(y ='ΔF/F0', x = '405 power (%)',
#        title = 'Fluo-4 exposure-effect',
#        subtitle = 'Gray area - corrected sd (sd / F/F0)') +
#   theme_minimal(base_size = 18)

