require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)
require(ggalt)

setwd('/home/astria/bio/lab/code/glt_anal')

##### DATA PREPROCESSING #####
df.exp <- read.csv('results.csv')
df.exp$day[df$day == ' 7d'] <- '7d'
df.exp$day[df$day == ' 14d'] <- '14d'
df.exp <- df.exp %>%
      mutate(group = as.factor(group)) %>%
      mutate(day = as.factor(day)) %>%
      subset((group == 'TBI' &
                day == '3d' &
                !(num %in% c(7, 12, 15, 2, 1, 14, 6, 13))) |
             (group == 'cef' & day == '3d' & !(num %in% c(9, 4, 1, 5))) |
             (group == 'TBI' & day == '14d' & !(num %in% c(3, 6))) |
             (group == 'cef' & day == '14d' & !(num %in% c(11, 9, 2))) |
             (group == 'cef' & day == '7d' & !(num %in% c(8, 1))) |
             (group == 'TBI' & day == '7d'))
df.ctrl <- read.csv('ctrl.csv') %>%
           subset(!(num %in% c(10, 11, 14)))

df <- bind_rows(df.exp, df.ctrl) %>%
      mutate(day =  ordered(day, levels = c('1d', '3d', '7d', '14d')))
remove(df.exp, df.ctrl)


##### DAY COMPARISON #####
days <- ggplot(df %>%
       filter(group != 'cont') %>%
       mutate(day = recode(day, '3d'='День 3', '7d'='День 7', '14d'='День 14')),
       aes(y = cluster, x = group, fill = group)) +
  geom_boxplot() +
  geom_point()+
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c('cef', 'TBI'))) +
  facet_grid(. ~ day) +
  scale_x_discrete(labels = c('Цеф.','ЧМТ')) +
  labs(fill = 'Група',
       x = 'Група',
       y = 'Кількість кластерів') +
  theme(legend.position = 'none',
        text = element_text(family='Times New Roman',
                            face="bold", size=20))

ggsave(sprintf('days.png', sel.group),
       days,
       width = 9,
       height = 4)


##### CTRL COMPARISON #####
sel.group <- 'cef'
comp <- ggplot(filter(df, group == 'cont' | group == sel.group),
       aes(x = day, y = cluster, fill = day)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c('1d', '3d'),
                                        c('1d', '7d'),
                                        c('1d', '14d'),
                                        c('3d', '7d'),
                                        c('3d', '14d'),
                                        c('7d', '14d'))) +
  scale_x_discrete(labels = c('Контроль','День 3','День 7', 'День 14')) +
  labs(fill = 'Група',
       x = 'Група',
       y = 'Кількість кластерів') +
  theme(legend.position = 'none',
        text = element_text(family='Times New Roman',
                            face="bold", size=20))

ggsave(sprintf('ctrl_%s.png', sel.group),
       comp,
       width = 5,
       height = 7)


##### INT vs CLUSTER #####
df.day_num <- df %>%
              mutate(day_num = recode(day, '1d'=1, '3d'=3, '7d'=7, '14d'=14))
i_c <- ggplot(df.day_num,
       aes(y = int, x = cluster, colour = group, fill=group, label = num)) +
  geom_point(size = 2) +
  geom_encircle(data = df.day_num %>% filter(group != "TBI"), alpha = .4) +
  geom_smooth(data = df.day_num %>% filter(group == 'TBI'),
              method = 'lm', alpha=.2, size = 1.5) + # geom_text(hjust = -0.5) +
  labs(fill = 'Група',
       colour = 'Група',
       x = 'Кількість кластерів',
       y = 'Сумма інтенсивності кластерів, a.u.') +
  scale_color_manual(labels = c('Цеф.', 'Контроль', 'ЧМТ'),
                     values = c('#F8766D', '#7CAE00', '#00BFC4')) +
  scale_fill_manual(labels = c('Цеф.', 'Контроль', 'ЧМТ'),
                    values = c('#F8766D', '#7CAE00', '#00BFC4')) +
  theme(text = element_text(family='Times New Roman',
                            face="bold", size=20))

ggsave('int_clst.png',
       i_c,
       width = 10,
       height = 7)


##### LM COMPARISON #####
lm.ctrl <- lm(int ~ cluster, data = filter(df, group == 'cont'))
lm.TBI <- lm(int ~ cluster, data = filter(df, group == 'TBI'))
lm.cef <- lm(int ~ cluster, data = filter(df, group == 'cef'))
summary(lm.cef)

b.TBI <- summary(lm.TBI)$coefficients[2,1]
se.TBI <- summary(lm.TBI)$coefficients[2,2]
b.cef <- summary(lm.cef)$coefficients[2,1]
se.cef <- summary(lm.cef)$coefficients[2,2]

compare.coeff <- function(b1,se1,b2,se2){
  return((b1-b2)/sqrt(se1^2+se2^2))
}
p_value = 2*pnorm(-abs(compare.coeff(b.TBI,se.TBI,b.cef,se.cef)))


##### MEDIAN #####
df.median <- df %>%
             subset(select = c('group', 'day', 'int', 'cluster')) %>%
             group_by(group, day) %>%
             summarise_all(list(mean, sd), na.rm = TRUE) %>%
             rename(int = int_fn1, int_sd = int_fn2,
                    clst = cluster_fn1, clst_sd = cluster_fn2)
write.csv(df.median, 'stat.csv')


