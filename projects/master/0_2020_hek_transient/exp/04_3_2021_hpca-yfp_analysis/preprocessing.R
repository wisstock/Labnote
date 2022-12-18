# HPCA-EYFP data analysis, uncaging with different stimulation power
# Copyright © 2021 Borys Olifirov
#
# Data frame structure:
# 'ID'     full unice recording ID
# 'date'   experiment date
# 'cell'   cell number with date suffix
# 'rep'    number of stimulation for this cell
# 'power'  405 nm laser power
# 'stimul' number of frame directly before stimulation
# 'frame'  number of frame
# 'time'   absolute time of frame recording, stimulation frame used as 0
# 'mask'   mask type
# 'int'    mean intensity of mask
# 'delta'  ΔF/F0 for int data
# 'rel'    ΔF/F0 for mask/cell mask, cell photo bleaching compensation

require(dplyr)
require(tidyr)
require(purrr)

require(TTR)
require(rstatix)

require(ggplot2)
require(ggpubr)
require(ggsci)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/04_3_2021_hpca-yfp_analysis')

#### DATA PREPROCESSING ####
df.fluo <- read.csv('results_fluo.csv') %>%
           subset(select = c(power, time, int)) %>%
           group_by(time, power) %>%
           summarise_all(list(mean, sd), na.rm = TRUE) %>%
           rename(mean = fn1, sd = fn2, frame = time) %>%
           mutate(time = frame - 3) %>%
           mutate(power = as.factor(power))

df.yfp.04 <- rbind(read.csv('results_3_04.csv')) %>%  # ,read.csv('results_17_03.csv')
          mutate(power = as.factor(power)) %>%
          mutate(cell = as.factor(cell)) %>%
          mutate(date = as.factor(date)) %>%
          mutate(rep = as.factor(rep)) %>%
          subset(select = c(ID, cell, date, rep, power, time, mask, delta, rel, int))

df.yfp.03 <- read.csv('results_17_03.csv') %>%
             mutate(power = as.factor(power)) %>%
             mutate(cell = as.factor(cell)) %>%
             mutate(date = as.factor(date)) %>%
             mutate(rep = as.factor(rep)) %>%
             subset(select = c(ID, cell, date, rep, power, time, mask, delta, rel, int))

# DROP BAD CELLS
df.yfp.04.drop <- df.yfp.04 %>%
  filter(ID == 'cell1_01_3_04' |  # 405 nm 50%
           ID == 'cell1_02_3_04' |
           ID == 'cell1_03_3_04' |
           ID == 'cell4_01_3_04' |
           ID == 'cell7_01_3_04' |
           ID == 'cell2_01_3_04' |  # 405 nm 75%
           ID == 'cell2_02_3_04' |
           ID == 'cell2_03_3_04' |
           ID == 'cell8_01_3_04' |
           ID == 'cell8_02_3_04' |
           ID == 'cell3_01_3_04' |  # 405 nm 100%
           ID == 'cell3_02_3_04' |
           ID == 'cell6_01_3_04' |
           ID == 'cell9_01_3_04')

df.yfp.03.drop <- df.yfp.03 %>%
  filter(ID == 'cell7_01_17_03'  |  # 405 nm 20%
           ID == 'cell9_01_17_03'  |
           ID == 'cell10_01_17_03')  # |  # 405 nm 50%
           # ID == 'cell12_01_17_03' |
           # ID == 'cell13_01_17_03') 
# ALL GOOD CELLS               
df.yfp <- rbind(df.yfp.03.drop, df.yfp.04.drop)        


##### RAW #####
# REPEATED STIMULATIONS
selected.power <- '100'
selected.mask <- 'up'

ggplot(filter(df.yfp, power == selected.power, mask == selected.mask)) +
  geom_vline(xintercept = 0) +
  geom_line(aes(x=time, y= rel, colour=ID),
            size=1.2) +
  scale_x_continuous(name = 'Time (s)',
                    limits = c(-20, 90),
                   breaks = seq(-100, 100, 5)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0',
                     breaks = seq(-100, 100, 0.1)) +
  labs(title = sprintf("Mask %s, 405 nm power %s", selected.mask, selected.power),
       colour = 'Recording ID') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_light() +
  theme(legend.position = 'bottom',
        text=element_text(family="Times New Roman", face="bold", size=12))

#  facet_grid(rows = vars(cell))

df.yfp.norm <- df.yfp %>%
  filter(power == selected.power,
         mask == selected.mask) %>%
  subset(select = c(ID, date, cell, mask, rep, time, rel)) %>%
  group_by(cell, rep, mask) %>%
  mutate(delta_norm = rel / max(rel))

# MAXIMUM INCREASE BOXPLOT
df.yfp.max <- df.yfp %>%
              subset(select = c(ID, power, time, mask, rel)) %>%
              group_by(power, mask, ID) %>%
              filter(rel == max(rel) & mask != 'cell')

kruskal.test(rel ~ power, data = filter(df.yfp.max, mask == 'up'))

ggplot(df.yfp.max) +
  geom_boxplot(aes(x = power, y = rel, fill = power),
               alpha = .6) +
  geom_point(aes(x = power, y = rel, color = power)) +
  facet_grid(cols = vars(mask)) +
  scale_color_jco() +
  scale_fill_jco() +
  theme_light()

# TRANS
df.trans <- df.yfp %>%
            filter(mask == 'up') %>%
            subset(select = c(ID, power, time, delta))

selected.power <- '75'

ggplot(filter(df.trans, power == selected.power)) +
  geom_line(aes(x = time, y = delta, colour = ID), size = 1.2) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(-20, 90),
                     breaks = seq(-100, 100, 5)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0',
                     breaks = seq(-100, 100, 0.1)) +
  labs(title = sprintf("Up/Cell, 405 nm power %s", selected.power),
       colour = 'Recording ID') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_light() +
  theme(legend.position = 'top',
        legend.justification = 'left')

df.trans.max <- df.trans %>%
                group_by(ID, power) %>%
                mutate(max_trans = max(delta)) %>%
                subset(power != '20', select = c(power, max_trans)) %>%
                unique()

kruskal.test(trans_delta ~ power, data = df.trans.delta)

ggplot(df.trans.max, aes(x = power, y = max_trans, fill = power)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c('50', '75'),
                                        c('75', '100'),
                                        c('50', '100'))) +
  scale_y_continuous(limits = c(0.15, .36),
                     breaks = seq(-1, 1, 0.05)) +
  labs(title = sprintf("Translocations", selected.mask),
       x = 'Power (%)', 
       y = 'ΔF',
       colour = '405 nm power (%)',
       fill = '405 nm power (%)') +
  scale_fill_jco() +
  theme_light()

# DOSE DEP
start.time <- 5
end.time <- 95
decim.num <- 10
selected.power <- '100'

# w.fun <- function(x, a, b){(exp(-x/a)+b)/max(exp(-x/a)+b)}
# a <- 30
# b <- 5
#
# ggscatter(data.frame(y = w.fun(seq(0, 100, 1), a, b), x = seq(0, 100, 1)),
#          x = 'x', y = 'y')

time.slice.yfp <- df.yfp.good %>%
                  filter(power == selected.power,
                         mask == selected.mask,
                         time>=start.time,
                         time<=end.time) %>%
                  subset(select = c(ID, time, delta)) %>%
                  group_by(ID) %>%
                  mutate(norm = rel/max(rel)) %>%
                  mutate(dec = HMA(norm, n = 10)) %>%  # w = w.fun(seq(0, length(norm)-1), a, b))
                  rename(yfp = delta)
time.slice.fluo <- df.fluo %>%
                   filter(power == selected.power,
                   time>=start.time,
                   time<=end.time) %>%
                   subset(select = c(mean, time)) %>%
                   rename(fluo = mean)

dose.dep.norm <- left_join(time.slice.yfp, time.slice.fluo,
                           by = c('time')) %>%
                 subset(select = c(ID, yfp, norm, dec, fluo))

ggplot() +
  geom_line(data = dose.dep.norm, 
            aes(x = fluo, y = dec, colour = ID),
            size = 1) +
  geom_line(data = dose.dep.norm, 
            aes(x = fluo, y = norm, colour = ID),
            size = 0.2,
            linetype = 5) +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     limits = c(0.5, 2),
                     breaks = seq(0, 10, .1)) +
  scale_y_continuous(name = 'Normalized HPCA-YFP ΔF/F0',
                     limits = c(-0.5, 1),
                     breaks = seq(-1, 1, .1)) +
  labs(title = sprintf("Dose dep, 405 nm power %s", selected.power),
       colour = 'Cell') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('dose_dep_raw_power_%s.png', selected.power),
#        width = 30, height = 15, units = 'cm')

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
            aes(x = time, y = mean, colour = power),
            size = 1) +
  geom_ribbon(data = filter(df.mean.yfp, mask == selected.mask),
              aes(x = time,
                  ymin = mean-sd,
                  ymax = mean+sd,
                  fill = as.factor(power),
                  colour = as.factor(power)),
              alpha=.2,
              size=0) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(-20, 90),
                     breaks = seq(-100, 100, 5)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0') +
  labs(title = sprintf("Mask %s mean", selected.mask),
       colour = '405 nm power (%)',
       fill = '405 nm power (%)') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_light() +
  theme(legend.position = 'top',
        legend.justification = 'left')

# MAXIMUM INCREASE BARPLOT
df.mean.max <- df.mean.yfp %>%
               group_by(power, mask) %>%
               filter(mean == max(mean) & mask != 'cell')

selected.mask <- 'up'
ggplot(filter(df.mean.max, mask == selected.mask, power != '20'),
       aes(x = power, y = mean, fill = power)) + 
  geom_bar(stat = "identity", alpha = .5) +
  geom_errorbar(aes(x = power, ymin = mean - sd, ymax = mean + sd), width = 0.25) +
  geom_point(data = filter(df.yfp.max, mask == selected.mask, power != '20'),
             aes(x = power, y = rel), size = 1.5) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
  labs(title = sprintf("Mask %s maximum insertion", selected.mask),
       x = 'Power (%)', 
       y = 'HPCA-YFP ΔF/F0',
       colour = '405 nm power (%)',
       fill = '405 nm power (%)') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_light()
  
               

#geom_line(data = subset(df.fluo, power == c('20', '50')),
#          aes(x = time, y = mean/coeff, colour = power, fill = power),
#          size = 1,
#          linetype = 5,
#          alpha = 0.6) +

#sec.axis = sec_axis(trans = ~.*coeff,
#                    name = 'Fluo-4 ΔF/F0',
#                    breaks = seq(-100, 100, 0.5))

# ggsave(sprintf('mean_mask_%s.png', selected.mask, selected.power),
#        width = 30, height = 15, units = 'cm')


# DOSE DEP
start.time <- 7
end.time <- 75

time.slice.yfp <- df.mean.yfp %>%
                  filter(mask == selected.mask,
                         time>=start.time,
                         time<=end.time) %>%
                  subset(select = c(power, mean, sd, time)) %>%
                  rename(mean_yfp = mean, sd_yfp = sd)
time.slice.fluo <- df.fluo %>%
                   filter(time>=start.time & time<=end.time) %>%
                   subset(select = c(power, mean, time)) %>%
                   rename(mean_fluo = mean)
dose.dep.mean <- left_join(time.slice.yfp, time.slice.fluo,
                            by = c('power', 'time')) %>%
                 subset(select = c(power, mean_yfp, sd_yfp, mean_fluo))
                 # group_by(power) %>%
                 # mutate(mean_yfp = (mean_yfp - min(mean_yfp))/(max(mean_yfp)-min(mean_yfp)))

ggplot(filter(dose.dep.mean, power != '75')) +
  geom_line(aes(x = mean_fluo, y = mean_yfp,
                colour = power),
            size = 1) +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     limits = c(0.2, 4),
                     breaks = seq(-100, 100, .5)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0',
                     limits = c(0, 0.45),
                     breaks = seq(-100, 100, 0.05)) +
  labs(title = "Dose dep mean",
       color = '405 nm power (%)') + 
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# HILL FIT
hill.mod <- lm(log10(mean_yfp/(1-mean_yfp)) ~ log10(mean_fluo),
               data = filter(dose.dep.mean, power == '50'))

model.n <- hill.mod$coefficients[[2]]
intercept <- hill.mod$coefficients[[1]]
model.ka <- intercept/model.n

df.hill <- dose.dep.mean %>%
           filter(power != '75') %>%
           group_by(power) %>%
           mutate(log_y = log10(mean_yfp/(1-mean_yfp))) %>%
           mutate(log_c = log10(mean_fluo)) %>%
           do(tidy(lm(log_y ~ log_c, .)))

ggplot(filter(dose.dep.mean, power != '75'),
       aes(y = log10(mean_yfp/(1-mean_yfp)),
           x = log10(mean_fluo),
           color = power,
           fill = power))+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = 'log(Fluo-4 ΔF/F0)',
       y = 'log(HPCA-YFP ΔF/F0 / (1-HPCA-YFP ΔF/F0))',
       title = 'Hill log form',
       color = '405 nm power (%)',
       fill = '405 nm power (%)') + 
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()
