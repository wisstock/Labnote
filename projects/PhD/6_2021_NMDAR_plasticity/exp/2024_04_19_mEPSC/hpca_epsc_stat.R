# mEPSC analysis, -70 and -40 LF stim with HPCA
# Copyright © 2024 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)
require(dunn.test)

setwd('/home/wisstock/bio/note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_04_19_mEPSC/-70_-40_HPCA')

a <- read.csv('260313_s7.csv') %>%
  mutate(TimRef = as.numeric(TimRef),
         Ref = as.numeric(Ref),
         TimPeak = as.numeric(TimPeak),
         Peak = as.numeric(Peak),
         TimEnd = as.numeric(TimEnd),
         End = as.numeric(End),
         IEI = as.numeric(IEI),
         Amp = as.numeric(Amp),
         Area = as.numeric(Area))

df.full <- bind_rows(read.csv('180213_s50.csv'), #  # r 
                     read.csv('110213_s23.csv'),
                     read.csv('180313_s30.csv'),
                     read.csv('190213_s22.csv'),
                     read.csv('190313_s49.csv'),
                     read.csv('210213_s22.csv'),
                     a) %>%
  mutate(Sample = as.factor(Sample),
         Group = as.factor(Group)) %>%
  mutate(Group = factor(Group, levels = c("-70ctrl", "-70mid", "-40mid", "-40post", "-70post"))) %>%
  filter(Amp <= 90.0)
remove(a)

##### AMPLITUDE ABS #####
df.amp <- df.full %>%
  select(Amp, Sample, Group, TimPeak)

# overview
ggplot(data = df.amp,
       aes(x = TimPeak, y = Amp, color = Group, alpha = .75)) +
  geom_point() +
  facet_wrap(facets = vars(Sample), nrow = nlevels(df.full$Sample), strip.position = 'right')

ggplot() +
  geom_boxplot(data = df.amp,
               aes(x = Group, y = Amp, fill = Sample,
                   group = interaction(Group, Sample))) +
  stat_summary(data = df.amp,
               fun = median,
               geom = 'line',
               aes(x = Group, y = Amp, color = Sample,
                   group = Sample),
               position = position_dodge(width = 0.75)) +
  stat_summary(data = df.amp,
               fun = median,
               geom = 'point',
               aes(x = Group, y = Amp, color = Sample,
                   group = Sample),
               alpha = 0.75,
               position = position_dodge(width = 0.75)) +
facet_wrap(facets = vars(Sample), ncol = nlevels(df.full$Sample), strip.position = 'right')

# stat by group
df.amp.test <- df.amp %>%
  filter((Group == '-70ctrl')|(Group == '-70mid')|(Group == '-70post'))
sample.list <- levels(df.amp.test$Sample)
for (selected.sample in sample.list) {
  print(selected.sample)
  dunn.test(df.amp.test$Amp[df.amp.test$Sample == selected.sample],
            df.amp.test$Group[df.amp.test$Sample == selected.sample],
            method = 'bh')
}

##### ANPLITUDE NORM #####
# norm overview
df.amp.norm <- df.amp %>%
  filter((Sample == '190313_s49')|  # (Sample == '210213_s22')|
         (Sample == '180313_s30')|
         (Sample == '190213_s22')) %>%
  droplevels() %>%
  group_by(Sample) %>%
  mutate(Amed = Amp / median(Amp[Group == '-70ctrl'])) %>%
  ungroup()

ggplot() +
  geom_boxplot(data = df.amp.norm,
               aes(x = Group, y = Amed, fill = Sample,
                   group = interaction(Group, Sample))) +
  stat_summary(data = df.amp.norm,
               fun = median,
               geom = 'line',
               aes(x = Group, y = Amed, color = Sample,
                   group = Sample),
               position = position_dodge(width = 0.75)) +
  stat_summary(data = df.amp.norm,
               fun = median,
               geom = 'point',
               aes(x = Group, y = Amed, color = Sample,
                   group = Sample),
               alpha = 0.75,
               position = position_dodge(width = 0.75))


df.amp.norm.test <- df.amp.norm %>%
  group_by(Sample, Group) %>%
  filter((Group == '-70ctrl')|(Group == '-70mid')|(Group == '-70post')) %>%
  summarise(med = mean(Amed), iqr = IQR(Amed), .groups="keep") %>%
  droplevels() %>%
  ungroup()

dunn.test(df.amp.norm.test$med,
          df.amp.norm.test$Group,
          method = 'bh')

df.amp.norm.stat <- df.amp.norm.summ %>%
  filter((Group == '-70ctrl')|(Group == '-70mid')) %>%
  droplevels() %>%
  select(-Sample, -iqr) %>%
  wilcox_test(med~Group)


df.amp.70 <- df.amp %>%
  select(-TimPeak) %>%
  filter((Group == '-70ctrl') | (Group == '-70post')) %>%
  droplevels()
df.amp.40 <- df.amp %>%
  select(-TimPeak) %>%
  filter((Group == '-40mid') | (Group == '-40post')) %>%
  select(Amp, Sample, Group) %>%
  droplevels()




# box stat
df.mid.amp.stat <- df.amp %>%
  filter((Group == '-70ctrl') | (Group == '-70mid')) %>%
  group_by(Sample) %>%
  wilcox_test(Amp~Group) %>%
  add_significance() %>%
  add_y_position(fun = 'median_iqr', scales = 'free')

df.post.amp.stat <- df.amp %>%
  filter((Group == '-70ctrl') | (Group == '-70post')) %>%
  group_by(Sample) %>%
  wilcox_test(Amp~Group) %>%
  add_significance() %>%
  add_y_position(fun = 'median_iqr', scales = 'free')

df.in.amp.stat <- df.amp %>%
  filter((Group == '-70mid') | (Group == '-70post')) %>%
  group_by(Sample) %>%
  wilcox_test(Amp~Group) %>%
  add_significance() %>%
  add_y_position(fun = 'median_iqr', scales = 'free')

ggplot(data = df.full,
       aes(x = Group, y = Amp, fill = Group)) +
  geom_boxplot() +
  facet_wrap(facets = vars(Sample), ncol = nlevels(df.full$Sample), strip.position = 'right')
  

##### IEI STAT #####
df.iei <- df.full %>%
  select(IEI, Sample, Group, TimPeak)


# prof
ggplot() +
  geom_point(data = df.iei,
             aes(x = TimPeak, y = IEI, color = Group, alpha = .75))


ggplot() +
  geom_boxplot(data = df.full,
               aes(x = Group, y = IEI, fill = Group))

##### AREA STAT #####