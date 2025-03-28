
require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)

setwd('/home/wisstock/bio_note/projects/PhD/10_SPIM_UC2/exp/25_03_28_T-SPIM_1st_beads')

df <- rbind(read.csv('bead_stat/bead_01.csv') %>% mutate(bead = '01'),
            read.csv('bead_stat/bead_02.csv') %>% mutate(bead = '02'),
            read.csv('bead_stat/bead_03.csv') %>% mutate(bead = '03'),
            read.csv('bead_stat/bead_04.csv') %>% mutate(bead = '04'),
            read.csv('bead_stat/bead_05.csv') %>% mutate(bead = '05'),
            read.csv('bead_stat/bead_06.csv') %>% mutate(bead = '06')) %>%
      mutate(bead = as.factor(bead)) %>%
      group_by(bead) %>%
      mutate(norm = (int-min(int))/(max(int)-min(int)),
             px_adj = case_when(bead == '01' ~ px-1,
                                bead == '02' ~ px+1,
                                bead == '03' ~ px-1,
                                bead == '04' ~ px-1,
                                bead == '05' ~ px-2,
                                bead == '06' ~ px,)) %>%
      ungroup()

px.size <- 5.56 / 4


##### GAUSSIAN FIT #####
df.peak <- df %>% filter(px_adj > -2, px < 14) %>%
  mutate(px_norm = px_adj - 6,
         px_size = px_norm * px.size)

px.median <- df.peak %>%
  group_by(px_adj) %>%
  mutate(norm_median = median(norm)) %>%
  ungroup() %>%
  select(px_size, norm_median) %>%
  distinct()

fitG = function(x,y,mu,sig,scale){
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    optim(c(mu,sig,scale),f)
  }

px.fit <- fitG(x = px.median$px_size, y = px.median$norm_median,
               mu = 0, sig = 2, scale = 1)
px.fit$par

fit_x <- seq(-10,10, 0.1)
fit_val <- dnorm(fit_x,
                 mean = px.fit$par[1],
                 sd = px.fit$par[2]) * px.fit$par[3]
plot(fit_x, fit_val)

df.g.fit <- data.frame(x = fit_x, y = fit_val)

##### PLOT #####
ggplot(df.peak, aes(x = px_size, y = norm)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_vline(xintercept = 2, linetype = 2) +
  geom_vline(xintercept = -2.1, linetype = 2) +
  geom_line(data = df.g.fit, aes(x = x, y = y), color = 'blue', size = 1.5) +
  geom_line(aes(group = bead), size = 0.25, color = 'grey30') +
  geom_point(aes(group = bead), size = 0.5, color = 'grey30') +
  # stat_summary(fun = median,
  #              geom = 'line', size = 1, color = 'grey25') +
  # stat_summary(fun = median,
  #              geom = 'point', size = 1.5, color = 'grey25') +
  # stat_summary(fun.min = function(z) { quantile(z,0.25) },
  #              fun.max = function(z) { quantile(z,0.75) },
  #              fun = median,
  #              geom = 'ribbon', size = 0, alpha = .25, color = 'grey25') +
  scale_x_continuous(breaks = seq(-10, 10, 2.5),
                     limits = c(-10, 10)) +
  theme_minimal() +
  theme(legend.position = 'none') +
  labs(y = 'I norm.', x = 'Distance, um',
       title = 'Yellow-green 100 nm beads (n=6)',
       subtitle = '4x 0.1 NA (expected lateral resolution ~3 um at 500 nm) \nEstimated FWHM ~4.1 um')
