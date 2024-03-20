# HPCA-TagRFP + Fluo-4 rim profile analysis
# Copyright © 2022 Borys Olifirov
#
# Profile data frame structure:
# 'ID' - recording ID
# 'time' - frame time (s)
# 'rim_d' - rim point mask ΔF/F
# 'dist' - distance from nucleus border

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/Bio/note/projects/PhD/6_2021_PIP2_HPCA_HEK/exp/2022_12_15_dist_and_FENS_pic/scripts')

##### GLOBAL PLOT OPTIONS #####
font_family <- "ubuntu mono"
font_size <- 30

line_size <- 1
line_size_light <- 1.5


##### DATA PREPROCESSING #####
df.fp <- read.csv('rim_fp_1px_02_2_2022.csv') %>%
  mutate(ch = 'fp')
df.ca <- read.csv('rim_ca_5px_02_2_2022.csv') %>%
  mutate(ch = 'ca')
df.total <- rbind(df.fp, df.ca) %>%
  mutate(ch = as.factor(ch),
         ID = as.factor(ID)) %>%
  mutate(i = as.numeric(i))
remove(df.fp, df.ca)

# d.line <- df.total %>%
#   select(ID, i, d) %>%
#   group_by(ID) %>%
#   unique() %>%
#   mutate(line = 1)

# 'cell2_02_2_2022'
# 'cell5_02_2_2022'
# 'cell6_02_2_2022'
# 'cell7_02_2_2022'

##### ONE CELL PLOT #####
sel.cell <- 'cell7_02_2_2022'
lab.line <- c('I', 'II', 'III')
lab.pos <- c(163, 195, 215)
sel.line <- c(c(159, 166), c(192, 198), c(212, 218))
stim.i <- c(14, 24, 34)

# distantion bar
dist.colors <- colorRampPalette(c('#FFFFFF', '#0000FF'))

dist.bar <- ggplot(d.line %>% filter(ID == sel.cell)) +
  geom_raster(aes(y = line, x = i, fill = d)) +
  scale_fill_gradientn(colours = dist.colors(2),
                       name="Distance (um)") +
  scale_x_continuous(name = 'Distance bar',
                     breaks = seq(0, 1000, 25),
                     expand = c(0,0)) +
  scale_y_reverse(name = '',
                  breaks = c(),
                  expand = c(0,0)) +
  theme(legend.box="horizontal",
        legend.position="bottom",
        panel.border = element_rect(fill = NA),
        plot.margin = unit(c(0, 0, -5, 132), "pt"),
        text = element_text(family=font_family,
                            face="bold",
                            size=font_size)) +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5))

# FP rim profile
gr.colors <- colorRampPalette(c('#08FF00', '#08F100', '#000000', '#EC0000', '#FF0000'))


df.cell.fp <- df.total %>%
  filter(ch == 'fp', ID == sel.cell)
fp.lim <- max(abs(df.cell.fp$delta))

fp.prof <- ggplot(df.total %>%
                    filter(ch == 'fp', ID == sel.cell)) +  #   facet_wrap(~ ID, scales = "free") +
  geom_raster(aes(x = i, y = time, fill = delta)) +
  geom_hline(yintercept = stim.i, color = 'white', size = line_size) +
  geom_vline(xintercept = sel.line,
             color = 'white', linetype = 'dashed', size = line_size) +
  theme(legend.position = "none") +
  scale_fill_gradientn(colours = gr.colors(5), 
                       limits=c(-fp.lim,fp.lim),
                       name="ΔF/F0") +
  scale_x_continuous(name = '',
                     position = "top",
                     breaks = c(),
                     expand = c(0,0)) +
  scale_y_reverse(name = 'Time (s)',
                  breaks = seq(0, 100, 4),
                  expand = c(0,0)) +
  theme(legend.position="left",
        panel.border = element_rect(fill = NA),
        plot.margin = unit(c(-13, 0, -5, 5.5), "pt"),
        text = element_text(family=font_family,
                            face="bold",
                            size=font_size)) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 15))

# Ca rim profile
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

ca.prof <- ggplot(df.total %>% filter(ch == 'ca', ID == sel.cell)) +
  geom_raster(aes(x = i, y = time, fill = delta)) +
  geom_hline(yintercept = stim.i, color = 'white', size = line_size) +
  geom_vline(xintercept = sel.line,
             color = 'white', linetype = 'dashed', size = line_size) +
  annotate("label", x = lab.pos,
           y = 2,
           label = lab.line,
           size = font_size-18) +
  theme(legend.position = "none") +
  scale_fill_gradientn(colours = jet.colors(7),
                       name="ΔF/F0") +
  scale_x_continuous(breaks = c(),
                     expand = c(0,0)) +
  scale_y_reverse(name = 'Time (s)',
                  breaks = seq(0, 100, 4),
                  expand = c(0,0)) +
  theme(legend.position="left",
        panel.border = element_rect(fill = NA),
        plot.margin = unit(c(0, 0, -5, 5.5), "pt"),
        text = element_text(family=font_family,
                            face="bold",
                            size=font_size)) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 15))

all.prof <- plot_grid(ca.prof, fp.prof, dist.bar,
                      ncol = 1,
                      rel_heights = c(0.4, 0.4, 0.15))
ggsave('all_rim.png', all.prof, 
       width = 36.3, height = 28, units = 'cm', dpi = 300)
# remove(all.prof, fp.prof, ca.prof, dist.bar, d.line)


##### STIMUL PDF #####
fp.red <- c('#FF0000', '#FF5117', '#FF6666')
ca.blue <- c('#0000FF', '#6767FF', '#00AEFF', '#366BFF')

stim.t <- c(6, 20, 46)

df.time <- df.total %>%
  filter(time %in% stim.t,
         ID == 'cell7_02_2_2022') %>%
  select(time, delta, ch)

# FP
# fp.rim <-
ggplot(df.time %>% filter(ch == 'fp'), aes(x = delta, fill = as.factor(time))) +
  geom_density(alpha = .5) +
  scale_fill_manual(name = NULL,
                    values = fp.red,
                    labels = c('Native','After 1st stimul','After 3d stimul')) +
  scale_x_continuous(name = 'HPCA-tagRFP ΔF/F0',
                     limits = c(-0.5, 2.5),
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'Density',
                     breaks = c(),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.7, 1),
        legend.justification = c("left", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('fp_rim.png', fp.rim, 
       width = 37.3, height = 10.6, units = 'cm', dpi = 300)
  
# Ca
ca.rim <- ggplot(df.time %>% filter(ch == 'ca'), aes(x = delta, fill = as.factor(time))) +
  geom_density(alpha = .5) +
  scale_fill_manual(name = NULL,
                    values = ca.blue,
                    labels = c('Native','After 1st stimul','After 3d stimul')) +
  scale_x_continuous(name = 'ΔF/F0',
                     limits = c(-0.5, 2.5),
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'Density',
                     breaks = c(),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.7, 1),
        legend.justification = c("left", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('ca_rim.png', ca.rim, 
       width = 37.3, height = 10.6, units = 'cm', dpi = 300)


##### SELECTED LINE #####
line.1 <- sel.line[1:2]
line.2 <- sel.line[3:4]
line.3 <- sel.line[5:6]
red.blue <- c('#0000FF', '#FF0000')

df.line <- df.total %>%
  filter(ID == sel.cell, i %in% sel.line) %>%
  mutate(band = case_when(i %in% line.1 ~ 'I',
                          i %in% line.2 ~ 'II',
                          i %in% line.3 ~ 'III')) %>%
  select(time, delta, band, ch) %>%
  group_by(time, band, ch) %>%
  mutate(band_mean = mean(delta),
         band_sd = sd(delta))

# band.grid <-
ggplot(df.line,
       aes(x = time, y = band_mean, color = ch, fill = ch)) +
  geom_vline(xintercept = stim.i) +
  geom_line(size = line_size) +
  geom_point() +
  geom_ribbon(aes(ymin = band_mean-band_sd,
                  ymax = band_mean+band_sd),
              alpha=.15,
              size=0.5) +
  facet_wrap(~band, nrow = 3) +
  scale_fill_manual(name = NULL,
                    values = red.blue,
                    labels = c('Fluo-4','HPCA-tagRFP')) +
  scale_color_manual(name = NULL,
                    values = red.blue,
                    labels = c('Fluo-4','HPCA-tagRFP')) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(0, 60),
                     breaks = seq(0, 100, 5),
                     expand = c(0,0)) +
  scale_y_continuous(name = "ΔF/F0",
                     breaks = seq(-100, 100, 1)) +
  scale_linetype_discrete(name="") +
  theme_classic() +
  theme(legend.box="horizontal",
        panel.grid.major = element_line(linetype = 'dotted',
                                        size = 0.1,
                                        colour = "black"),
        legend.position = c(.23, 1),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('band.png', band.grid, 
       width = 36.3, height = 17.2, units = 'cm', dpi = 300)
# remove(all.prof, fp.prof, ca.prof, dist.bar, d.line)


##### DD #####
dd.total <- df.total %>%
  select(-time, -d) %>%
  group_by(ID, i) %>%
  pivot_wider(names_from = ch, values_from = delta, values_fn = median)


# rep.dd <-
  ggplot(dd.total %>% filter(ID == 'cell7_02_2_2022'),
                 aes(x = ca, y = fp)) +
  geom_point(color = 'grey', alpha = .25) +
  geom_smooth(color = 'red', fill = 'red',
              method = 'loess') +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'HPCA ΔF/F0',
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        linewidth = 0.1,
                                        colour = "black"),
        legend.position = "none",
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('rep_dd.png', rep.dd, 
       width = 18.7, height = 8, units = 'cm')

##### DD WITH DIST #####
df.d.dd <- df.total %>%
  filter(ID == 'cell7_02_2_2022')

bin_num <- 3  # round(1 + log(length(df.total$d)), digits = 0)
d_interval <- round((max(df.d.dd$d) - min(df.d.dd$d)) / bin_num, digits = 3)
bin_list <- seq(0, bin_num, 1) * d_interval
last_bin <- sprintf("(%.2f,%.2f]", tail(bin_list, n=1), max(df.d.dd$d))

d.dd.total <- df.d.dd %>%
  mutate(d.interval = cut(d, breaks = bin_list)) %>%
  select(-time, -d) %>%
  group_by(ID, i, d.interval) %>%
  pivot_wider(names_from = ch, values_from = delta, values_fn = median) %>%
  ungroup() %>%
  select(-i) %>%
  mutate(d.interval = as.character(d.interval)) %>%
  mutate(d.interval = replace_na(d.interval, last_bin)) %>%
  mutate(d.interval = as.factor(d.interval)) %>%
  filter(ca <= 2) %>%
  drop_na()

### RAW PLOT
dd_inteval_raw <- ggplot(d.dd.total,
       aes(x = ca, y = fp, fill = d.interval, color = d.interval)) +
  geom_point(alpha = .25) +
  geom_smooth(method = 'loess') +  # facet_wrap(~ ID, scales = 'free') +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'HPCA ΔF/F0',
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        linewidth = 0.1,
                                        colour = "black"),
        legend.position = c(.15, .75),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('rep_dd_interval.png', dd_inteval_raw, 
       width = 28, height = 15, units = 'cm')


### LM FORM
d.dd.linear <- d.dd.total %>%
  group_by(ID, d.interval) %>%
  mutate(log_y = log10(fp/(1-fp))) %>%
  mutate(log_c = log10(ca))

d.dd.lm <- d.dd.linear %>%
  drop_na() %>%
  do(tidy(lm(log_y ~ log_c, .)))
  
# model.n <- hill.mod$coefficients[[2]]
# intercept <- hill.mod$coefficients[[1]]
# model.ka <- intercept/model.n

### LM PLOT
dd_inteval_lm <- ggplot(d.dd.linear,
       aes(x = log_c, y = log_y, fill = d.interval, color = d.interval)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = .3) +
  geom_smooth(method = 'lm') +
  scale_x_continuous(limits = c(-0.5, 0.3)) +
  scale_y_continuous(limits = c(-3, 1)) +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = 'dotted',
                                        linewidth = 0.1,
                                        colour = "black"),
        legend.position = c(-1, .75),
        text = element_text(family=font_family, face="bold", size=font_size))

ggsave('rep_dd_interval_lm.png', dd_inteval_lm, 
       width = 28, height = 15, units = 'cm')
  
### n PLOT
ggplot(d.dd.lm %>% filter(term == 'log_c' & estimate > 0),
       aes(x = ID, y = estimate, color = d.interval)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymax = estimate + std.error,
                    ymin = estimate - std.error),
                width = 0.1)
