# HPCA-EYFP data analysis, uncaging with different stimulation power
# Copyright © 2021 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(ggplot2)
require(ggsci)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/1_03_2021_hpca-yfp_analysis/data')

df.fluo <- read.csv('results_fluo.csv') %>%  # 20.10.2020, results of test uncaging, stimulation after 5th frame
           subset(select = c(power, time, int)) %>%
           group_by(time, power) %>%
           summarise_all(list(mean, sd), na.rm = TRUE) %>%
           rename(mean = fn1, sd = fn2, frame = time)
df.fluo$time <- as.numeric(df.fluo$frame - 3)
df.fluo$power <- as.factor(df.fluo$power)

# rbind(read.csv('results_22_01_1-0.csv'),  # results 22.01.2021, 1.0 s/frame
#       read.csv('results_27_01_1-0.csv'),  # results 27.01.2021, 1.0 s/frame
#       read.csv('results_17_03_1-0.csv')) %>%  # results 17.03.2021, 1.0 s/frame

df.yfp <- read.csv('results_17_03.csv') %>%
          filter(cell == 'cell7_01_17_03' |
                 cell == 'cell9_01_17_03' |
                 cell == 'cell10_01_17_03' |
                 cell == 'cell12_01_17_03' |
                 cell == 'cell13_01_17_03' |
                 cell == 'cell20_01_17_03')
df.yfp$power <- as.factor(df.yfp$power)


##### MEAN #####
# TIME SERIES
selected.mask <- 'up'
coeff <- 7

df.mean.yfp <- subset(df.yfp, select = c(power, time, mask, rel)) %>%
               group_by(time, power, mask) %>%
               summarise_all(list(mean, sd), na.rm = TRUE) %>%
               rename(mean = fn1, sd = fn2)

  ggplot() +
  geom_line(data = filter(df.mean.yfp, mask == selected.mask),
            aes(x = time, y = mean, colour = power, fill = power),
            size = 1) +
  geom_ribbon(data = filter(df.mean.yfp, mask == selected.mask),
              aes(x = time,
                  ymin = mean-sd,
                  ymax = mean+sd,
                  fill = as.factor(power),
                  colour = as.factor(power)),
              alpha=.2,
              size=0) +
  geom_line(data = subset(df.fluo, power == c('20', '50')),
            aes(x = time, y = mean/coeff, colour = power, fill = power),
            size = 1,
            linetype = 5,
            alpha = 0.6) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(-20, 90),
                     breaks = seq(-100, 100, 5)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0',
                     breaks = seq(-100, 100, 0.05),
                     sec.axis = sec_axis(trans = ~.*coeff,
                                       name = 'Fluo-4 ΔF/F0',
                                       breaks = seq(-100, 100, 0.5))) +
  labs(title = sprintf("Mask %s mean", selected.mask),
       colour = '405 nm power (%)',
       fill = '405 nm power (%)') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('mean_mask_%s.png', selected.mask, selected.power),
#        width = 30, height = 15, units = 'cm')


# DOSE DEP
start.time <- 4
end.time <- 100

time.slice.yfp <- df.mean.yfp %>%
                  filter(mask == selected.mask,
                         time>=start.time,
                         time<=end.time) %>%
                  subset(select = c(power, mean, time)) %>%
                  rename(mean_yfp = mean)
time.slice.fluo <- df.fluo %>%
                   filter(time>=start.time & time<=end.time) %>%
                   subset(select = c(power, mean, time)) %>%
                   rename(mean_fluo = mean)
dose.dep.mean <- left_join(time.slice.yfp, time.slice.fluo,
                            by = c('power', 'time')) %>%
                 subset(select = c(power, mean_yfp, mean_fluo))
dose.dep.rel <- dose.dep.mean

write.csv(file = 'DR_50_17_03.csv', dose.dep.rel)

ggplot(dose.dep.rel, aes(x = mean_fluo, y = mean_yfp, colour = as.factor(power))) +
  geom_line(size = 1) +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     limits = c(0.2, 2.1),
                     breaks = seq(-100, 100, .1)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0',
                     limits = c(0, 0.35),
                     breaks = seq(-100, 100, 0.05)) +
  labs(title = "Dose dep mean",
       color = '405 nm power (%)') + 
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('dose_dep_mean.png', selected.power),
#        width = 30, height = 15, units = 'cm')


##### RAW #####
# TIME SERIES
selected.power <- 50
selected.mask <- 'up'

df.yfp.norm <- df.yfp %>%
               filter(power == selected.power) %>%
               subset(select = c(cell, time, rel, mask)) %>%
               group_by(cell, mask) %>%
               mutate(norm = rel/max(rel))

ggplot() +
  geom_line(data = filter(df.yfp.norm, mask == selected.mask),
            aes(x=time, y= norm, colour=as.factor(cell)),
            size=1) +
  geom_line(data=subset(df.fluo, (power == selected.power)),
            aes(x=time, y=mean/max(mean)),
            size = 1.2,
            linetype = 5,
            color = 'black',
            alpha = .6) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(-20, 90),
                     breaks = seq(-100, 100, 5)) +
  scale_y_continuous(name = 'Normalized ΔF/F0',
                     limits = c(-0.1, 1),
                     breaks = seq(-1, 1, .1)) +
  labs(title = sprintf("Mask %s, 405 nm power %s", selected.mask, selected.power),
       colour = 'Cell') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('raw_mask_%s_power_%s.png', selected.mask, selected.power),
#        width = 30, height = 15, units = 'cm')

# DOSE DEP
start.time <- 5
end.time <- 95

time.slice.yfp <- df.yfp.norm %>%
                  filter(mask == selected.mask,
                         time>=start.time,
                         time<=end.time) %>%
                  subset(select = c(cell, norm, time)) %>%
                  rename(int_yfp = norm)
time.slice.fluo <- df.fluo %>%
                   filter(power == selected.power,
                          time>=start.time,
                          time<=end.time) %>%
                   subset(select = c(mean, time)) %>%
                   rename(int_fluo = mean)
dose.dep.norm <- left_join(time.slice.yfp, time.slice.fluo,
                           by = c('time')) %>%
                 subset(select = c(cell, int_yfp, int_fluo))

ggplot() +
  geom_line(data = dose.dep.norm, 
            aes(x = int_fluo, y = int_yfp, colour = cell),
            size = 1) +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     limits = c(0.3, 2),
                     breaks = seq(0, 10, .1)) +
  scale_y_continuous(name = 'Normalized HPCA-YFP ΔF/F0',
                     limits = c(-0.1, 1),
                     breaks = seq(-1, 1, .1)) +
  labs(title = sprintf("Dose dep, 405 nm power %s", selected.power),
       colour = 'Cell') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('dose_dep_raw_power_%s.png', selected.power),
#        width = 30, height = 15, units = 'cm')
