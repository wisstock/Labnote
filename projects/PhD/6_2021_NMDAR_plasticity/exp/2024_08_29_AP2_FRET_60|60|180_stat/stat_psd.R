# NMDA ionophoresis, HPCA traslocation in spines
# Copyright © 2024 Borys Olifirov, Oleksandra Fedchenko

require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_08_29_AP2_FRET_60|60|180_stat')


##### DATA PREPROCESSING #####
# psd
df.60 <- bind_rows(read.csv('./psd_sasha/cell1_HPCA_in_PSD95.csv'),
                   read.csv('./psd_sasha/cell2_HPCA_in_PSD95.csv')) %>%
  mutate(protocol = '60s')
df.5 <- read.csv('./psd_sasha/cell3_HPCA_in_PSD95.csv') %>%
  mutate(protocol = '5s')
df.10 <- read.csv('./psd_sasha/cell4_HPCA_in_PSD95.csv') %>%
  mutate(protocol = '10s')
df.15 <- read.csv('./psd_sasha/cell5_HPCA_in_PSD95.csv') %>%
  mutate(protocol = '15s')

df.02.Hz <- bind_rows(df.60,
                      df.5,
                      df.10,
                      df.15) %>%
  mutate(time = time * 5,
         rel_time = (index-5) * 5,
         site_type = as.factor('PSD'),
         cell_id = id) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id)) %>%
  select(-X)
remove(df.60, df.5, df.10, df.15)

df.15 <- read.csv('./psd_sasha/cell6_HPCA_in_PSD95.csv') %>%
  mutate(protocol = '15s',
         rel_time = (time * 2) - 25)
df.60 <- read.csv('./psd_sasha/cell7_HPCA_in_PSD95.csv') %>%
  mutate(protocol = '60s',
         rel_time = (index-5) * 2)

df.05.Hz <- bind_rows(df.60,
                      df.15) %>%
  mutate(time = time * 2,
         site_type = as.factor('PSD'),
         cell_id = id) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id)) %>%
  select(-X)
remove(df.60, df.15)

df.psd <- bind_rows(df.02.Hz,
                    df.05.Hz)
remove(df.02.Hz, df.05.Hz)

# oreol
df.60 <- bind_rows(read.csv('./psd_sasha/cell1_HPCA_in_PSD_oreol.csv'),
                   read.csv('./psd_sasha/cell2_HPCA_in_PSD_oreol.csv')) %>%
  mutate(protocol = '60s')
df.5 <- read.csv('./psd_sasha/cell3_HPCA_in_PSD_oreol.csv') %>%
  mutate(protocol = '5s')
df.10 <- read.csv('./psd_sasha/cell4_HPCA_in_PSD_oreol.csv') %>%
  mutate(protocol = '10s')
df.15 <- read.csv('./psd_sasha/cell5_HPCA_in_PSD_oreol.csv') %>%
  mutate(protocol = '15s')

df.02.Hz <- bind_rows(df.60,
                      df.5,
                      df.10,
                      df.15) %>%
  mutate(time = time * 5,
         rel_time = (index-5) * 5,
         site_type = as.factor('Oreol'),
         cell_id = id) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id)) %>%
  select(-X)
remove(df.60, df.5, df.10, df.15)

df.15 <- read.csv('./psd_sasha/cell6_HPCA_in_PSD_oreol.csv') %>%
  mutate(protocol = '15s',
         rel_time = (time * 2) - 25)
df.60 <- read.csv('./psd_sasha/cell7_HPCA_in_PSD_oreol.csv') %>%
  mutate(protocol = '60s',
         rel_time = (index-5) * 2)

df.05.Hz <- bind_rows(df.60,
                      df.15) %>%
  mutate(time = time * 2,
         site_type = as.factor('Oreol'),
         cell_id = id) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id)) %>%
  select(-X)
remove(df.60, df.15)

df.oreol <- bind_rows(df.02.Hz,
                      df.05.Hz)
remove(df.02.Hz, df.05.Hz)

# shaft
df.60 <- bind_rows(read.csv('./psd_sasha/cell1_HPCA_in_shaft.csv'),
                   read.csv('./psd_sasha/cell2_HPCA_in_shaft.csv')) %>%
  mutate(protocol = '60s')
df.5 <- read.csv('./psd_sasha/cell3_HPCA_in_shaft.csv') %>%
  mutate(protocol = '5s')
df.10 <- read.csv('./psd_sasha/cell4_HPCA_in_shaft.csv') %>%
  mutate(protocol = '10s')
df.15 <- read.csv('./psd_sasha/cell5_HPCA_in_shaft.csv') %>%
  mutate(protocol = '15s')

df.02.Hz <- bind_rows(df.60,
                      df.5,
                      df.10,
                      df.15) %>%
  mutate(time = time * 5,
         rel_time = (index-5) * 5,
         site_type = as.factor('Shaft'),
         cell_id = id) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id)) %>%
  select(-X)
remove(df.60, df.5, df.10, df.15)

df.15 <- read.csv('./psd_sasha/cell6_HPCA_in_shaft.csv') %>%
  mutate(protocol = '15s',
         rel_time = (time * 2) - 25)
df.60 <- read.csv('./psd_sasha/cell7_HPCA_in_shaft.csv') %>%
  mutate(protocol = '60s',
         rel_time = (index-5) * 2)

df.05.Hz <- bind_rows(df.60,
                      df.15) %>%
  mutate(time = time * 2,
         site_type = as.factor('Shaft'),
         cell_id = id) %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id)) %>%
  select(-X)
remove(df.60, df.15)

df.shaft <- bind_rows(df.02.Hz,
                      df.05.Hz)
remove(df.02.Hz, df.05.Hz)

df.full <- bind_rows(df.psd, df.oreol, df.shaft) %>%
  mutate(protocol = as.factor(protocol),
         cell_id = as.factor(cell_id),
         protocol = factor(protocol, c('5s', '10s', '15s', '60s'), ordered = TRUE))
remove(df.psd, df.oreol, df.shaft)
write.csv(df.full, file = 'psd_sasha.csv')

###### CTRL PLOTS #####
df.sites.plot <- df.full %>%
  group_by(cell_id, index, site_type) %>%
  mutate(int_med = median(int)) %>%
  ungroup() %>%
  select(site_type, int_med, index, rel_time, cell_id)

ggplot(data = df.full,
       aes(x = rel_time, y = int,
           color = cell_id, fill = cell_id,
           group = cell_id)) +
  stat_summary(fun = median,
               geom = 'line', size = 0.5) +
  stat_summary(fun = median,
               geom = 'point', size = 0.5) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .15) +
  annotate('rect',
           xmin = 2, xmax = 8,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, color = 'black', size = 0) +
  geom_vline(xintercept = 0, lty =2) +
  theme_classic() +
  labs(color = 'Cell ID',
       fill = 'Cell ID',
       x = 'Time, s',
       y = expression(ΔF/F[0])) +
  facet_wrap(~site_type+protocol, nrow = 3)


##### HIGH AND LOW ROIS #####
calc.mixmdl <- function(input_vector) {
  mixmdl <- normalmixEM(input_vector)
  input_dens <- density(input_vector)
  comp1 <- dnorm(x = input_dens$x,
                 mean = mixmdl$mu[1],
                 sd = mixmdl$sigma[1]) * mixmdl$lambda[1]
  comp2 <- dnorm(x = input_dens$x,
                 mean = mixmdl$mu[2],
                 sd = mixmdl$sigma[2]) * mixmdl$lambda[2]
  return(list(data.frame(val = input_dens$x,
                         raw = input_dens$y,
                         comp1 = comp1,
                         comp2 = comp2,
                         comp_comb = comp1+comp2),
              mixmdl))
}

df.zero.gmm.optim <- bind_rows(data.frame('AIC' = c(-356.22427190298845,
                                                    -393.28511084735896,
                                                    -390.9785059885147,
                                                    -388.02984719016854,
                                                    -381.97731567833415,
                                                    -379.0864807891088,
                                                    -378.55286650449705,
                                                    -382.86313089972316,
                                                    -378.03845948137337,
                                                    -375.6694049228803),
                                          'BIC' = c(-349.34811328514206,
                                                    -376.094714302743,
                                                    -363.47387151712917,
                                                    -350.2109747920134,
                                                    -333.8442053534094,
                                                    -320.63913253741447,
                                                    -309.79128032603313,
                                                    -303.78730679448967,
                                                    -288.6483974493703,
                                                    -275.9651049641076),
                                          'site_type' = as.factor('Shaft')),
                               data.frame('AIC' = c(-968.2260725044782,
                                                    -1100.2451231661007,
                                                    -1092.6582612456436,
                                                    -1093.2391543894623,
                                                    -1087.9474910969343,
                                                    -1092.7763098260775,
                                                    -1086.830643423723,
                                                    -1075.2268182779221,
                                                    -1075.116724776967,
                                                    -1081.1473613680782),
                                          'BIC' = c(-959.5453538970228,
                                                    -1078.543326647462,
                                                    -1057.9353868158216,
                                                    -1045.495202048457,
                                                    -1027.1824608447457,
                                                    -1018.9902016627058,
                                                    -1000.0234573491679,
                                                    -975.3985542921838,
                                                    -962.2673828800455,
                                                    -955.2769415599735),
                                          'site_type' = as.factor('PSD')),
                               data.frame('AIC' = c(-479.6944097286916,
                                                    -559.1299534419587,
                                                    -556.0633332836582,
                                                    -562.8328278529899,
                                                    -557.5268463906223,
                                                    -562.8841556606095,
                                                    -556.893817381008,
                                                    -551.5276175892523,
                                                    -556.7372960313539,
                                                    -554.4587623515083),
                                          'BIC' = c(-472.13910508224626,
                                                    -540.2416918258455,
                                                    -525.8421146978769,
                                                    -521.2786522975407,
                                                    -504.6397138655051,
                                                    -498.66406616582435,
                                                    -481.3407709165549,
                                                    -464.64161415513115,
                                                    -458.51833562756485,
                                                    -444.90684497805125),
                                          'site_type' = as.factor('Oreol')))

df.zero <- df.full %>%
  filter(rel_time >= 2 & rel_time <= 7) %>%
  group_by(index, roi_id, site_type) %>%
  mutate(int = median(int)) %>%
  ungroup() %>%
  select(roi_id, int, cell_id, protocol, site_type)
  
ggplot() +
  geom_density(data = df.zero,
               aes(x = int, fill = site_type), alpha = .5)


# plot(calc.mixmdl.optim(df.zero$int[df.zero$site_type == 'Shaft']), type = 'line')

pm <- calc.mixmdl(df.zero$int[df.zero$site_type == 'PSD']) 
df.psd.mixmdl <- pm[[1]]
psd.mixmgl <- pm[[2]]
om <- calc.mixmdl(df.zero$int[df.zero$site_type == 'Oreol']) 
df.oreol.mixmdl <- om[[1]]
oreol.mixmgl <- om[[2]]
sm <- calc.mixmdl(df.zero$int[df.zero$site_type == 'Shaft']) 
df.shaft.mixmdl <- sm[[1]]
shaft.mixmgl <- sm[[2]]
remove(pm, om, sm)


# PSD HIST
psd.hist <- ggplot() +
  geom_histogram(data = df.zero %>% filter(site_type == 'PSD'),
                 aes(x = int),
                 color = 'black',
                 fill = 'magenta',
                 alpha = .5) +
    geom_line(data = df.psd.mixmdl,
              aes(x = val, y = comp1 * 30),
              lty = 2, size = 0.75) +
    geom_line(data = df.psd.mixmdl,
              aes(x = val, y = comp2 * 30),
              lty = 2, size = 0.75) +
    geom_line(data = df.psd.mixmdl,
              aes(x = val, y = comp_comb * 30),
              size = 0.75) +
  scale_x_continuous(limits = c(-0.25, 0.6)) +
  theme_classic() +
  labs(title = 'PSD',
       x = expression(ΔF/F[0]),
       y = '# ROIs')

psd.bic <- ggplot() +
  geom_line(data = df.zero.gmm.optim %>% filter(site_type == 'PSD'),
            aes(x = seq(1,10), y = BIC), color = 'magenta') +
  geom_point(data = df.zero.gmm.optim %>% filter(site_type == 'PSD'),
             aes(x = seq(1,10), y = BIC), color = 'magenta') +
  scale_x_continuous(breaks = seq(1,10)) + 
  theme_classic() +
  labs(x = '# components',
       y = 'BIC')

psd.gmm.plot <- ggdraw(psd.hist) +
  draw_plot(psd.bic, x = 0.6, y = 0.5, width = 0.3, height = 0.4) +
  draw_plot_label(c("A", "Aa"),
                  c(0, 0.55),
                  c(1, 0.92),
                  size = 14)

# PSD HIST
oreol.hist <- ggplot() +
  geom_histogram(data = df.zero %>% filter(site_type == 'Oreol'),
                 aes(x = int),
                 color = 'black',
                 fill = 'orange',
                 alpha = .5) +
  geom_line(data = df.oreol.mixmdl,
            aes(x = val, y = comp1 * 13),
            lty = 2, size = 0.75) +
  geom_line(data = df.oreol.mixmdl,
            aes(x = val, y = comp2 * 13),
            lty = 2, size = 0.75) +
  geom_line(data = df.oreol.mixmdl,
            aes(x = val, y = comp_comb * 13),
            size = 0.75) +
  scale_x_continuous(limits = c(-0.25, 0.6)) +
  theme_classic() +
  labs(title = 'Oreol',
       x = expression(ΔF/F[0]),
       y = '# ROIs')

oreol.bic <- ggplot() +
  geom_line(data = df.zero.gmm.optim %>% filter(site_type == 'Oreol'),
            aes(x = seq(1,10), y = BIC), color = 'orange') +
  geom_point(data = df.zero.gmm.optim %>% filter(site_type == 'Oreol'),
             aes(x = seq(1,10), y = BIC), color = 'orange') +
  scale_x_continuous(breaks = seq(1,10)) + 
  theme_classic() +
  labs(x = '# components',
       y = 'BIC')

orl.gmm.plot <- ggdraw(oreol.hist) +
  draw_plot(oreol.bic, x = 0.6, y = 0.5, width = 0.3, height = 0.4) +
  draw_plot_label(c("B", "Ba"),
                  c(0, 0.55),
                  c(1, 0.92),
                  size = 14)

# SHAFT HIST
shaft.hist <- ggplot() +
  geom_histogram(data = df.zero %>% filter(site_type == 'Shaft'),
                 aes(x = int),
                 color = 'black',
                 fill = 'cyan',
                 alpha = .5) +
  geom_line(data = df.shaft.mixmdl,
            aes(x = val, y = comp1 * 13),
            lty = 2, size = 0.75) +
  geom_line(data = df.shaft.mixmdl,
            aes(x = val, y = comp2 * 13),
            lty = 2, size = 0.75) +
  geom_line(data = df.shaft.mixmdl,
            aes(x = val, y = comp_comb * 13),
            size = 0.75) +
  scale_x_continuous(limits = c(-0.25, 0.6)) +
  theme_classic() +
  labs(title = 'Shaft',
       x = expression(ΔF/F[0]),
       y = '# ROIs')

shaft.bic <- ggplot() +
  geom_line(data = df.zero.gmm.optim %>% filter(site_type == 'Shaft'),
            aes(x = seq(1,10), y = BIC), color = 'cyan') +
  geom_point(data = df.zero.gmm.optim %>% filter(site_type == 'Shaft'),
             aes(x = seq(1,10), y = BIC), color = 'cyan') +
  scale_x_continuous(breaks = seq(1,10)) + 
  theme_classic() +
  labs(x = '# components',
       y = 'BIC')

shf.gmm.plot <- ggdraw(shaft.hist) +
  draw_plot(shaft.bic, x = 0.6, y = 0.5, width = 0.3, height = 0.4) +
  draw_plot_label(c("C", "Ca"),
                  c(0, 0.55),
                  c(1, 0.92),
                  size = 14)

plot_grid(psd.gmm.plot, orl.gmm.plot, shf.gmm.plot, nrow = 3)


psd.threshold <- 0.18
oreol.threshold <- 0.3
shaft.threshold <- 0.32

df.full <- df.full %>%
  group_by(roi_id, site_type) %>%
  mutate(LH = case_when((site_type == 'PSD') & (int >= psd.threshold)  ~ 'H',
                        (site_type == 'PSD') & (int < psd.threshold)  ~ 'L',
                        (site_type == 'Oreol') & (int >= oreol.threshold)  ~ 'H',
                        (site_type == 'Oreol') & (int < oreol.threshold)  ~ 'L',
                        (site_type == 'Shaft') & (int >= shaft.threshold)  ~ 'H',
                        (site_type == 'Shaft') & (int < shaft.threshold)  ~ 'L',
                        .default = 'out')) %>%
  filter(LH != 'out') %>%
  droplevels()

