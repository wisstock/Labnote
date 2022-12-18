# Fluo-4 data analysis
# Copyright © 2021 Borys Olifirov
#
# Require data frame structure:
# 'cell'   cell number with date suffix
# 'rep'    number of stimulation for this cell
# 'power'  405 nm laser power
# 'time'   absolute time of frame recording, stimulation frame used as 0
# 'delta'  ΔF/F0 for int data

require(dplyr)
require(tidyr)
require(purrr)

require(rstatix)

require(ggplot2)
require(ggpubr)
require(cowplot)

setwd('/home/astria/bio/note/projects/PhD/6_2021_PIP2_HPCA_HEK/exp/11_22_2021_test_uncaging')


##### DATA PROCESSING #####
df.not.40 <- read.csv('results_10_4.csv') %>%
             filter(power != '40')

df.full <- read.csv('results_09_29.csv') %>%
           rbind(read.csv('results_11_13.csv')) %>%
           rbind(read.csv('results_11_21.csv')) %>%  # bad cellsrbind(read.csv('results_11_22.csv')) %>%
           rbind(df.not.40) %>%
           mutate(power = as.factor(power))

df.fluo <- df.full %>%
           filter(power != '90'&
           power != '95' &
           power != '100') %>%
           droplevels()

df.fluo.abs.baseline <- df.fluo %>%
                        select(power, time, int) %>%
                        filter(time < 0) %>%
                        group_by(time, power) %>%
                        summarise_all(list(mean, sd), na.rm = TRUE) %>%
                        rename(mean = fn1, sd = fn2) %>%
                        ungroup()

df.fluo.mean <- df.fluo %>%
                select(power, time, delta) %>%
                group_by(time, power) %>%
                summarise_all(list(mean, sd), na.rm = TRUE) %>%
                rename(mean = fn1, sd = fn2) %>%
                ungroup()

time.list <- c(0, 20, 80)
df.fluo.crop <- df.fluo %>%
  select(power, delta, time) %>%
  filter(time %in% time.list) %>%
  group_by(time) %>%
  unique() %>%
  ungroup()

df.fluo.stat.group <- df.fluo.crop %>%
  group_by(time) %>%
  kruskal_test(delta ~ power) %>%
  add_significance("p")

# Bonferoni-Holm adjastment
# https://www.statisticshowto.com/holm-bonferroni-method/
# https://rpkgs.datanovia.com/rstatix/
df.fluo.stat.pair <- df.fluo.crop %>%
  group_by(time) %>%
  pairwise_wilcox_test(delta ~ power,
                       p.adjust.method = 'BH') %>%
  add_xy_position(fun = "mean_sd", scales = 'free')

df.fluo.summary <- df.fluo %>%
                  filter(time == 0) %>%
                  group_by(power) %>%
                  summarise(N = n(),
                            Mean = round(mean(delta), digits = 2),
                            IQR = round(IQR(delta), digits = 2)) %>%
                  rename('Power(%)' = power, 'Max ΔF/F0' = Mean)


##### GLOBAL PLOT OPTIONS #####
font_family <- "ubuntu mono"
font_size <- 15

line_size <- 1
line_size_light <- 1.5

# https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2
# '#121212', '#FF1C1C', '#0051F2', '#FFE20A', '#219900', '#FF7B08', '#008F8C', '#B02179', '#A6A6A6'
palette <- c('#121212', '#FF1C1C', '#0051F2', '#FFE20A', '#219900', '#FF7B08', '#008F8C', '#B02179', '#A6A6A6')

power.list <- levels(df.fluo$power)
col.list <- palette[1:length(power.list)]
names(col.list) <- power.list


##### PLOT #####
# line plot
line.plt <- ggplot(df.fluo.mean) +
  geom_line(aes(x = time, y = mean, colour = power),
            size = line_size) +
  geom_ribbon(aes(x = time,
                  ymin = mean-sd,
                  ymax = mean+sd,
                  fill = power,
                  colour = power),
              alpha=.15,
              size=0) +
  geom_vline(xintercept = time.list) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(-20, 90),
                     breaks = seq(-100, 100, 10),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'ΔF/F0',
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  labs(colour = 'Power (%)', fill='Power (%)') +
  theme_classic() +
  theme(legend.position = c(.90, .95),
        legend.justification = c("right", "top"),
        text = element_text(family=font_family, face="bold", size=font_size))


# abs intensity line plot
abs.plt <- ggplot(df.fluo.abs.baseline) +
  geom_line(aes(x = time, y = mean, colour = power),
            size = line_size) +
  geom_ribbon(aes(x = time,
                  ymin = mean-sd,
                  ymax = mean+sd,
                  fill = power,
                  colour = power),
              alpha=.15,
              size=0) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(-20, 0),
                     breaks = seq(-100, 100, 5),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'a.u.',
                     # limits = c(50, 250),
                     breaks = seq(-100, 500, 25),
                     expand = c(0,0)) +
  labs(colour = 'Power (%)', fill='Power (%)') +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(family=font_family, face="bold", size=font_size))


# https://www.datanovia.com/en/blog/how-to-add-p-values-to-ggplot-facets/
# boxplot          
point.box <- ggplot(df.fluo.crop, aes(x = power, y = delta)) +
  geom_boxplot(aes(fill = power), alpha = .6) +
  geom_point() +
  facet_grid(cols = vars(as.factor(time))) +
  labs(fill = 'Power (%)',
       x = 'Power (%)',
       y = 'ΔF/F0') +
  scale_y_continuous(breaks = seq(-100, 100, 1)) +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(family=font_family, face="bold", size=font_size)) +
  stat_pvalue_manual(df.fluo.stat.pair, label = 'p.adj.signif',
                     family=font_family, size = 5,
                     tip.length = 0,
                     hide.ns = TRUE) +
  stat_compare_means(family=font_family, size=3,
                     label.x = 1.7)


# dose-dep plot
dd.plt <- ggplot(df.fluo.mean %>%
       filter(time == '0') %>%
       select(power, mean, sd) %>%
       mutate(power = as.numeric(as.character(power))),
       aes(x = power, y = mean)) +
  geom_line(size = line_size, colour = '#A6A6A6') +
  geom_point() +
  geom_errorbar(aes(ymin = mean-sd,
                    ymax = mean+sd),
                size=.4,
                width=3) +
  scale_x_continuous(name = 'Power (%)',
                     limits = c(0, 102),
                     breaks = seq(0, 100, 10),
                     expand = c(0,0)) +
  scale_y_continuous(name = 'ΔF/F0',
                     limits = c(0, 4.75),
                     breaks = seq(-100, 100, .5),
                     expand = c(0,0)) +
  theme_classic() +
  theme(text = element_text(family=font_family, face="bold", size=font_size))

# summary table
summ.table <- ggtexttable(df.fluo.summary,
                          rows = NULL,
                          theme = ttheme("light"))
              # tab_add_title(text = 'Groups size \n and maximal ΔF/F0',
              #               face = "bold",
              #               size = 10)


# combine plot with cowplot
fin_draw <- ggdraw() +
  draw_plot(line.plt, x = 0, y = .66, width = .8, height = .33) +
  draw_plot(abs.plt, x = .75, y = .66, width = .2, height = .33) +
  draw_plot(point.box, x = 0, y = .33, width = 1, height = .33) +
  draw_plot(dd.plt, x = 0, y = 0, width = .75, height = .33) +
  draw_plot(summ.table, x = .75, y = 0, width = .25, height = .33) +
  draw_plot_label(label = c('A', 'B', 'C', 'D', 'E'), size = 15,
                  x = c(0, .75, 0, 0, .75), y = c(1, 1, .66, .33, .33))

ggsave('fin_draw.png',
       fin_draw,
       width = 10,
       height = 10)


# # combine plot with ggpubr
# line.box <- ggarrange(line.plt, point.box,  # abs.plt 
#                       labels = c("A", "B"),
#                       ncol = 1, nrow = 2)
# 
# dd.table <- ggarrange(dd.plt, summ.table,
#                     labels = c("C", "D"),
#                     ncol = 2, nrow = 1,
#                     widths=c(0.75, 0.25))
# 
# fin.plt <- ggarrange(line.box, dd.table,
#                      ncol = 1, nrow = 2,
#                      heights=c(0.8, 0.2))
# 
# ggsave('fin_plt.png',
#        fin.plt,
#        width = 10,
#        height = 10)
