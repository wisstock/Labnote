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

require(ggplot2)
require(ggpubr)
require(dplyr)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/1_03_2021_hpca-yfp_analysis/data')

Ka <- 0.35
n <- 2

start.val <- 0
end.val <- 1
exp.mean <- 0.1  # mean value for exponential distributed C data
obs <- 60  # number of observations
hill.sd <- 0.01  # additive gaussian noise sd
hill <- function(x){((x^n)/(Ka^n+x^n))}  # Hill equation function

#### Model data ####
model.hill <- data.frame(c = seq(start.val, end.val, (end.val-start.val)/obs)) %>%
  mutate(y = hill(c)) %>%
  mutate(c_noise = c + rnorm(obs+1, mean = 0, sd = hill.sd)) %>%  # add noise to C
  mutate(y_noise = y + rnorm(obs+1, mean = 0, sd = hill.sd)) %>%  # add noise to Y
  mutate(exp_c = rexp(obs+1, rate = 1/exp.mean)) %>%
  mutate(exp_y = hill(exp_c) + rnorm(obs+1, mean = 0, sd = hill.sd))

#### Dose-response plot building ####
DRP <- ggplot(model.hill) +
  stat_function(data = data.frame(x = c(start.val, end.val)),
                aes(x),
                fun = function(x){(x^n)/(Ka^n+x^n)},
                color = 'red',
                size = 1) +
  geom_point(aes(x = c_noise, y = y_noise),
             colour = 'black',
             alpha = .6,
             size = 2) +
  geom_point(aes(x = exp_c, y = exp_y),
             colour = 'blue',
             alpha = .6,
             size = 2.5) +
  scale_x_continuous(name = 'C',
                     limits = c(0, 1)) +
  scale_y_continuous(name = 'Y',
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  labs(title = sprintf('Dose-response plot, Ka=%s, n=%s, SD=%s', Ka, n, hill.sd),
       subtitle = sprintf('red curve - theoretical curve,\nblack - uniform distributed C,\nblue - exp. distributed C (mean=%s)', exp.mean)) +
  theme_minimal()

#### Linear model building ####
# uniform distributed data
hill.uni.lm <- glm(log10(y_noise/(1-y_noise)) ~ log10(c_noise), data = model.hill)
uni.model.n <- hill.uni.lm$coefficients[[2]]
uni.intercept <- hill.uni.lm$coefficients[[1]]
uni.model.ka <- 10^(uni.intercept/uni.model.n)  # calculate model Ka from intercep (n*log(Ka))

# exponential distributed data
hill.exp.lm <- glm(log10(exp_y/(1-exp_y)) ~ log10(exp_c), data = model.hill)
exp.model.n <- hill.exp.lm$coefficients[[2]]
exp.intercept <- hill.exp.lm$coefficients[[1]]
exp.model.ka <- 10^(exp.intercept/exp.model.n)

#### Hill plot building ####
HP <- ggplot(model.hill) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = log10(c), y = log10(y/(1-y))),
            color = 'red',
            size = 1) +
  stat_smooth(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise))),
              method = 'lm',
              color = 'black',
              size = 1.2,
              alpha = 0.6) +
  stat_smooth(aes(x = log10(exp_c), y = log10(exp_y/(1-exp_y))),
              method = 'lm',
              color = 'blue',
              size = 1.2,
              alpha = 0.6,
              linetype = 2) +
  geom_point(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise))),
             alpha = .6) +
  geom_point(aes(x = log10(exp_c), y = log10(exp_y/(1-exp_y))),
             alpha = .6,
             colour = 'blue') +
  scale_x_continuous(breaks = seq(-10, 10, 0.25)) +  # limits = c(-2, 0),
  scale_y_continuous(breaks = seq(-10, 10, 0.5)) +  # limits = c(-3.5, 1),
  labs(title = 'Hill plot',
       subtitle = sprintf('uniform C model Ka=%.2f, n=%.2f \nexponential C model Ka=%.2f, n=%.2f',
                          uni.model.ka, uni.model.n, exp.model.ka, exp.model.n), 
       x = 'log(C)',
       y = 'log(Y/1-Y)') +
  theme_minimal()

#### Final plot ####
ggarrange(DRP, HP, ncol = 2, nrow = 1)  # print two panel plot
summary(hill.uni.lm)  # print linear model data in console
summary(hill.exp.lm)
