# Flyo-4 and Fura Red data analysis, uncaging with different stimulation power
# Copyright © 2020 Borys Olifirov

require(ggplot2)
require(dplyr)
# require(magrittr)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/12_12_2020_without_np/scripts')

df <- read.csv('results_8.csv')
f_0 <- mean(df$int[df$time >= 1 & df$time <= 15 & df$cell == 'cell8_01'])

deltaF <- function(f, f0){return((f-f0)/f0)}

df.01 <- subset(df, cell == 'cell8_01')
df.02 <- subset(df, cell == 'cell8_02')
df.03 <- subset(df, cell == 'cell8_03')
df.04 <- subset(df, cell == 'cell8_04')

df.01$delta <- deltaF(df.01$int, f_0)
df.02$delta <- deltaF(df.02$int, f_0)
df.03$delta <- deltaF(df.03$int, f_0)
df.04$delta <- deltaF(df.04$int, f_0)

l.size <- 1

ggplot() +
  geom_line(data = df.01, aes(x = time, y = delta, color = cell), size = l.size) +
  geom_line(data = df.02, aes(x = time, y = delta, color = cell), size = l.size) +
  geom_line(data = df.03, aes(x = time, y = delta, color = cell), size = l.size) +
  geom_line(data = df.04, aes(x = time, y = delta, color = cell), size = l.size) +
  scale_x_continuous(limits = c(-10, 180),
                     breaks = seq(0, 180, 20)) +
  labs(y ='ΔF/F0', x = 'Time (s)',
       title = 'Cell 8 ΔF/F0') +
  theme_minimal(base_size = 18)