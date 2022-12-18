# Experimental data decimation and fitting.
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

set.seed(100)
setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/03_26_2021_hill_analysis/data')

df.experimental <- read.csv('DR_50_17_03.csv') %>%
  filter(power == 50) %>%
  subset(select = c(mean_yfp, mean_fluo)) %>%
  rename(yfp = mean_yfp, fluo = mean_fluo) %>%
  mutate(yfp_norm = yfp/max(yfp)) %>%
  mutate(fluo_norm = fluo/max(fluo))

#### Exponential Fluo-4 transient generating ####
shapiro.test(log(df.experimental$fluo))
trans.mean <- mean(log(df.experimental$fluo))
trans.sd <- sd(log(df.experimental$fluo))

df.exp <- data.frame(obs = df.experimental$fluo_norm,
                     exp = rexp(length(df.experimental$fluo_norm),
                                       rate = 1/mean(df.experimental$fluo_norm)))


ggplot(df.exp) +
  geom_density(aes(x = obs),
                 fill = 'blue',
                 colour = 'blue',
                 alpha = .5) +
  geom_density(aes(x = exp),
                 fill = 'red',
                 colour = 'red',
                 alpha = .5)

#### Transient lover limit ####
# all extimations based om normilized fluo data
low.limits <- seq(0., 0.8, 0.05)
df.crop.residuals <- data.frame(a = as.numeric(),
                                intercept = as.numeric(),
                                residuals_sd = as.numeric(),
                                R_adj = as.numeric(),
                                lower_limit = as.numeric())

for (lim in low.limits){
  df.selected.lim <- df.experimental %>%
    subset(select = c(yfp_norm, fluo_norm)) %>%
    filter(fluo_norm >= lim) %>%
    mutate(fluo_norm = log10(fluo_norm)) %>%
    mutate(yfp_norm = log10(yfp_norm/(1-yfp_norm))) %>%
    filter(yfp_norm != Inf)
  lim.lm <- lm(yfp_norm ~ fluo_norm,
               data = df.selected.lim)
  
  df.crop.residuals <- rbind(df.crop.residuals,
                             data.frame(a = lim.lm$coefficients[2],
                                        intercept = lim.lm$coefficients[1],
                                        residuals_sd = sigma(lim.lm),
                                        R_adj = summary(lim.lm)$adj.r.squared,
                                        lower_limit = lim))
}
a <- ggplot(df.crop.residuals) +
  geom_line(aes(x = lower_limit, y = a),
            colour = 'blue') +
  geom_point(aes(x = lower_limit, y = a),
             colour = 'blue') +
  labs(y = 'n',
       x = 'Lover limit') +
  theme_minimal()

inter <- ggplot(df.crop.residuals) +
  geom_line(aes(x = lower_limit, y = intercept),
            colour = 'blue') +
  geom_point(aes(x = lower_limit, y = intercept),
             colour = 'blue') +
  labs(y = 'Intercept',
       x = 'Lover limit') +
  theme_minimal()

R <- ggplot(df.crop.residuals) +
  geom_line(aes(x = lower_limit, y = R_adj),
            colour = 'blue') +
  geom_point(aes(x = lower_limit, y = R_adj),
             colour = 'blue') +
  labs(y = 'Adjasted R',
       x = 'Lover limit') +
  theme_minimal()

resSD <- ggplot(df.crop.residuals) +
  geom_line(aes(x = lower_limit, y = residuals_sd),
            colour = 'blue') +
  geom_point(aes(x = lower_limit, y = residuals_sd),
             colour = 'blue') +
  labs(y = 'Residuals SD',
       x = 'Lover limit') +
  theme_minimal()
#ggarrange(ggarrange(R, resSD, ncol = 2, nrow = 1),
#          ggscatter(data = df.experimental, x = 'fluo_norm', y = 'yfp_norm')+
#            theme_minimal(),
#          ncol = 1, nrow = 2)
ggarrange(a, inter, R, resSD, ncol = 1, nrow = 4)


#### Transient decimation ####
#group.num <- 15
df.dec <- df.experimental %>%
          subset(select = c(yfp_norm, fluo_norm)) %>%
          filter(fluo_norm >= .6) %>%
          mutate(fluo_norm = log10(fluo_norm)) %>%
          mutate(yfp_norm = log10(yfp_norm/(1-yfp_norm))) %>%

          #mutate(group = ntile(fluo_norm, group.num)) %>%
          #group_by(group) %>%
          #summarise_all(list(mean, sd), na.rm = TRUE) %>%
          #rename(yfp_mean = yfp_norm_fn1,
          #       yfp_sd = yfp_norm_fn2,
          #       fluo_mean = fluo_norm_fn1,
          #       fluo_sd = fluo_norm_fn2)
          # filter(abs(fluo - mean(fluo)) == min(abs(fluo - mean(fluo))))
  
ggplot() +
  stat_smooth(data = df.experimental,
              aes(x = log10(fluo_norm), y = log10(yfp_norm/(1-yfp_norm))),
              method = 'lm',
              colour = 'blue',
              alpha = .3) +
  geom_point(data = df.experimental,
             aes(x = log10(fluo_norm), y = log10(yfp_norm/(1-yfp_norm))),
             colour = 'blue',
             alpha = .3) +
  stat_smooth(data = df.dec,
              aes(x = log10(fluo_norm), y = log10(yfp_norm/(1-yfp_norm))),
              method = 'lm',
              colour = 'red',
              alpha = .3) +
  geom_point(data = df.dec,
             aes(x = log10(fluo_norm), y = log10(yfp_norm/(1-yfp_norm))),
             color = 'red',
             alpha = .3,
             size = 4) +
  theme_minimal()




