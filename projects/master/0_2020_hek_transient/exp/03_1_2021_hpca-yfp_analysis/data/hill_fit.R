# Experiments with Hill equation.
# Copyright © 2021 Borys Olifirov
#
#       [C]^n
# Y = -----------                   natural form
#    Ka^n + [C]^n
#
#      y
# log(---) = log([C]) + n*log(Ka)   linear form
#     1-Y
#
# Y - fraction of occupied target protein
# [C] - ligand (Ca2+) concentration
# n - Hill coeficient
# Ka - ligand concentration producing half occupation
#
# Kd = Ka^n, Kd - dissociation constant
setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/1_03_2021_hpca-yfp_analysis/data')

df.experimental <- read.csv('DR_50_17_03.csv') %>%
  filter(power == 50) %>%
  subset(select = c(mean_yfp, mean_fluo)) %>%
  rename(yfp = mean_yfp, fluo = mean_fluo) %>%
  mutate(norm_yfp = yfp/max(yfp))

ggscatter(df.experimental, x = 'fluo', y = 'yfp')
max.fluo <- round(max(df.experimental$fluo), digit = 1)
min.fluo <- round(min(df.experimental$fluo), digit =1)
max.yfp <- round(max(df.experimental$yfp), digit =2)
min.yfp <- round(min(df.experimental$yfp), digit =2)


Ka <- 0.5
n <- 4

start.val <- 0
end.val <- 1
obs <- 20  # number of observations
hill.sd <- 0.01  # additive gaussian noise sd
hill <- function(x){((x^n)/(Ka^n+x^n))}  # Hill equation function

model.hill <- data.frame(c = seq(start.val, end.val, (end.val-start.val)/obs)) %>%
  mutate(y = hill(c)) %>%
  mutate(c_noise = c + rnorm(obs+1, mean = 0, sd = hill.sd)) %>%
  mutate(y_noise = y + rnorm(obs+1, mean = 0, sd = hill.sd))


#### Dose-response plot building ####
scale <- 1
y.shift <- 0 
x.shift <- 0

ggplot() +
  geom_point(data = df.experimental,
             aes(x = fluo, y = yfp), # fluo+x.shift, y = (yfp/scale)+y.shift),
             color = 'blue')

DRP <- ggplot(model.hill) +
  geom_point(aes(x = c_noise, y = y_noise), alpha = .6) +
  stat_function(data = data.frame(x = c(start.val, end.val)),
                aes(x),
                fun = function(x){(x^n)/(Ka^n+x^n)},
                color = 'red',
                size = 1) +
  scale_x_continuous(name = 'C',
                     limits = c(0, 1)) +
  scale_y_continuous(name = 'Y',
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  labs(title = sprintf('Dose-response plot, Ka=%s, n=%s', Ka, n),
       subtitle = sprintf('red - theoretical curve,\ndots - generated data + noise (sd=%s)', hill.sd), 
       x = 'Fluo-4 ΔF/F0') +
  theme_minimal()

ggplot() +
  stat_function(data = data.frame(x = c(start.val, end.val)),
                aes(x),
                fun = function(x){(x^n)/(Ka^n+x^n)},
                color = 'red',
                size = 1) +
  geom_point(data = df.experimental,
             aes(x = fluo, y = norm_yfp), # fluo+x.shift, y = (yfp/scale)+y.shift),
             color = 'blue') +
  scale_x_continuous(name = 'C',
                     limits = c(0, 1),
                     sec.axis = sec_axis(trans = ~.*1,
                                         breaks = seq(-(min.fluo*5), max.fluo*5, 0.1),
                                         name = 'Fluo-4 ΔF/F0')) +
  scale_y_continuous(name = 'Y',
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     sec.axis = sec_axis(trans = ~.*scale,
                                         breaks = seq(-1*scale, 1*scale, 0.05),
                                         name = 'HPCA-YFP ΔF/F0')) +
  labs(title = sprintf('Dose-response plot, Ka=%s, n=%s', Ka, n),
       subtitle = sprintf('red - theoretical curve,\ndots - generated data + noise (sd=%s)', hill.sd), 
       x = 'Fluo-4 ΔF/F0') +
  theme_minimal()

#### Linear model building ####
hill.mod <- lm(log10(y_noise/(1-y_noise)) ~ log10(c_noise), data = model.hill)
model.n <- hill.mod$coefficients[[2]]
intercept <- hill.mod$coefficients[[1]]
model.ka <- 10^(intercept/model.n)
# df.hill.mod <- data.frame('Parameter' = c('Ka', 'n'), 'Val.' = c(model.ka, model.n))
# model.res <- ggtexttable(df.hill.mod, rows = NULL, theme = ttheme("blank"))

#### Hill plot building ####
HP <- ggplot(model.hill) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = log10(c), y = log10(y/(1-y))),
            color = 'red',
            size = 1) +
  stat_smooth(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise))),
              method = 'lm',
              color = 'blue',
              size = 1.2,
              alpha = 0.6) +
  geom_point(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise))),
             alpha = .6) +
  scale_x_continuous(breaks = seq(-10, 10, 0.25)) +  # limits = c(-2, 0),
  scale_y_continuous(breaks = seq(-10, 10, 0.25)) +  # limits = c(-3.5, 1),
  labs(title = 'Hill plot',
       subtitle = sprintf('blue - linear model \nmodel Ka=%.2f, model n=%.2f', model.ka, model.n), 
       x = 'log(C)',
       y = 'log(Y/1-Y)') +
  theme_minimal()


ggarrange(DRP, HP, ncol = 2, nrow = 1)
