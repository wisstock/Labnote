# HPCA-EYFP data analysis, uncaging with different stimulation power
# & Hill equation fitting
# Copyright © 2021 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(ggplot2)
require(ggsci)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/03_26_2021_hill_analysis/data')

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
           cell == 'cell20_01_17_03')  # regect bad cells
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

dose.dep.50 <- subset(dose.dep.rel, power == '50',
                      select = c(mean_yfp, mean_fluo)) %>%
  mutate(norm_yfp = mean_yfp/max(mean_yfp)) %>%
  mutate(log_yfp = log10(norm_yfp/(1-norm_yfp))) %>%
  mutate(log_fluo = log10(mean_fluo))
dose.dep.50[dose.dep.50 == Inf] <- 0.0

hill.mod <- lm(log_yfp ~ log_fluo, data = dose.dep.50)
print(hill.mod)

ggplot(dose.dep.50, aes(x = log_fluo, y = log_yfp)) +
  geom_point() +
  stat_smooth(method = 'lm')

# ggsave(sprintf('dose_dep_mean.png', selected.power),
#        width = 30, height = 15, units = 'cm')

##### HILL #####
Ka <- 0.5
n <- 2
obs <- 20
hill.sd <- 0.01
hill <- function(x){((x^n)/(Ka^n+x^n))}

model.hill <- data.frame(c = seq(0, 1, 1/obs)) %>%
              mutate(y = hill(c)) %>%
              mutate(c_noise = c + rnorm(obs+1, mean = 0, sd = hill.sd)) %>%
              mutate(y_noise = y + rnorm(obs+1, mean = 0, sd = hill.sd))

ggplot(model.hill) +
  geom_point(aes(x = c_noise, y = y_noise)) +
  geom_line(aes(x = c, y = y), color = 'red') +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  labs(title = sprintf('Dose-responce plot, Ka=%s, n=%s', Ka, n),
     subtitle = sprintf('Red - model data without noise, dots - model data + gaussian noise (sd=%s)', hill.sd), 
     x = 'C',
     y = 'Y') +
  theme_minimal()

# linear model building
hill.mod <- lm(log10(y_noise/(1-y_noise)) ~ log10(c_noise), data = model.hill)
summary(hill.mod)

ggplot(model.hill) + 
  geom_point(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise)))) +
  geom_line(aes(x = log10(c), y = log10(y/(1-y))), color = 'red') +
  stat_smooth(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise))),
              method = 'lm') +
  scale_x_continuous(limits = c(-1, 0),
                     breaks = seq(-10, 10, 0.1)) +
  scale_y_continuous(limits = c(-3.5, 1),
                     breaks = seq(-10, 10, 0.5)) +
  labs(title = sprintf('Hill plot, Ka = %s, n = %s', Ka, n),
       subtitle = 'Red - model data without noise, blue - linear model', 
       x = 'log(C)',
       y = 'log(Y/1-Y)') +
  theme_minimal()

ggplot() +
  stat_function(data = data.frame(x = c(0, c)), aes(x),
                 fun = function(x){(x^n)/(Ka^n+x^n)}) +
  scale_x_continuous(limits = c(0, c),
                     breaks = seq(0, c, c/10)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, c, c/10)) +
  labs(title = sprintf('Ka = %s, n = %s', Ka, n),
       x = 'Ligand concentration',
       y = 'Binded target') +
  theme_minimal()



