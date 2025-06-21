require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
library(gridExtra)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_06_19_Sasa_HPCA-PSD/data')
df


f0_idx <- seq(15,30)

#LOAD FILE
all_cells <- read.csv('E:/Data_new/all_cells_hpca+psd95/all_cells_hpca_psd95.csv')
spine <- read.csv('E:/Data_new/all_cells_hpca+psd95/sorted/spines/cell15_ch0_spine.csv')%>%
  filter(time > 15) %>%
  mutate(roi = as.factor(roi)) %>%
  group_by(roi) %>%
  mutate(
    f0 = mean(abs_int[time %in% f0_idx], na.rm = TRUE),  # F₀ per ROI
    df = (abs_int - f0) / f0,
    app = "0.5")
spine <- spine %>%
  mutate(roi = as.factor(roi),
         rel_time = time-30) %>%
  group_by(roi) 


#parameters for plot
font.size <- 15
font.fam <- 'Arial'

#PLOTS
#checking ROIs
ggplot(data = raw_shaft, aes(x = time, y = df, color = roi)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~ roi, scales = "free_y") +
  theme(legend.position = "none") +
  labs(title = "ROI traces", x = "Time", y = "df")

#Plot standart for poster
ggplot(spine, aes(x=rel_time, y=df, color = roi))+
  annotate('rect',
           xmin = 0, xmax = 0.5,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  geom_hline(yintercept = 0,  lty = 2) +
  geom_line(size = 1.5) +
  geom_point() +
  geom_line(size = 2) +
  scale_color_manual(values = c('1' = 'magenta',
                                '2' = 'yellow3',
                                '3' = 'orange3')) +
  labs(title = "HPCA_mBaoJin_sec_0.5app_spine", x = "Time", y = expression(Delta * F / F[0]), color = 'ROIs') +
  theme_classic() 
  theme(legend.position = c(0.9, 0.9),
        legend.title = element_text(size = font.size),
        text=element_text(size = font.size, family = font.fam, face = 'bold')) +
  facet_wrap(~id, nrow = 4, scales = 'free_y', strip.position = 'right')

ggsave("HPCA_mBaoJin_sec.jpg", plot = df1, width = 20, height = 15)

all_cells <- bind_rows(all_cells, df4_60_1, df4_60_2)
df4_60_1$app <- as.numeric(df4_60_1$app)
#DIST VS AMP
ggplot(
  all_cells %>% filter(lab_id != "psd" & id == "cell7_ch0")%>%
    filter((time >= 50 & time <= 60) |
             (time >= 70 & time <= 120)|
             (time >= 200 & time <= 230)) %>%
    mutate(group_amp = case_when(
      time >= 50 & time <= 60 ~ "before",
      time >= 70 & time <= 120 ~ "peak",
      time >= 200 & time <= 230 ~ "after"))) +
  geom_hline(yintercept = 0,  lty = 2) +
  geom_point(size = 3, alpha = 0.7, aes(x = dist, y = df, color = lab_id, shape = lab_id, alpha = 0.1)) +
  geom_smooth(method = 'lm', aes(x = dist, y = df, color = lab_id, fill = lab_id)) +
  facet_wrap(~ index) +
  scale_y_continuous(limits = c(-.2,0.3)) +
  theme(legend.position = "none", text=element_text(size = font.size, family = font.fam, face = 'bold')) +
  labs(title = "DIST vs AMP_", x = "dist", y = expression(Delta * F / F[0])) 

dist_peaks <- peaks %>%
  mutate(dist_group = case_when(
    dist < 100 ~ 'max',
    dist >= 100 & dist < 300 ~ 'mid',
    dist >= 300 & dist ~ 'min',
    TRUE ~ NA_character_  # in case there are distances > 400
  ))

dist_peaks <- dist_peaks %>%
  group_by(id, roi, dist_group) %>%
  mutate(median_int_peak = median(int_peak, na.rm = TRUE)) %>%
  ungroup()

ggplot(data = dist_peaks %>% filter(app == '60') %>%
         mutate(app = as.numeric(as.character(app))),  
       aes(x = dist_group, y = median_int_peak,
           color = lab_id, group = lab_id)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(fun = median, geom = 'line', size = 1.25) +
  stat_summary(fun = median, geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75) +
  scale_color_manual(values = c('oreol' = 'yellow3', 'shaft' = 'blue3', 'psd' = 'magenta')) +
  theme_classic() +
  theme(legend.position = 'right',
        text = element_text(size = font.size+3, family = font.fam, face = 'bold'),
        plot.caption = element_text(size = font.size+4),
        legend.title = element_text(size = font.size - 2)) +
  facet_wrap(~id, nrow = 3, scales = 'free_y', strip.position = 'right') +
  labs(title = "HPCA_mBaoJin_amp_peaks_at_60sec_app_time",
       color = 'ROItype',
       x = 'Distance group',
       y = expression(Delta * F / F[0]))

#PEAKS
#dataframe
apps <- c('0.5', '2.5', '5', '10', '20', '30', '60', '120')

peaks <- all_cells %>%
  select(id, lab_id, int_peak, app, roi, dist) %>%
  filter(!is.na(int_peak), app %in% apps) %>%  
  mutate(app = factor(app, levels = apps, ordered = TRUE)) %>%
  group_by(id, app, lab_id, roi) %>%
  slice(1) %>%  
  ungroup() %>%
  droplevels()


peaks <- peaks %>% 
  filter(int_peak >= 0.01)

#stat
kruskal.test(int_peak~app,  data = peaks %>% filter(lab_id == 'shaft'))

peak_stat_df_oreol<- peaks %>%
  ungroup() %>%
  filter(lab_id == 'oreol') %>%
  pairwise_wilcox_test(int_peak~app, p.adjust.method = 'BH',
                       group.by = "lab_id") %>%
  add_xy_position()
#boxplot
ggplot(peaks %>% filter(lab_id == 'oreol'), aes(x=app, y=int_peak))+
  geom_boxplot(aes(fill = app), position = position_dodge(width = 0.5), alpha = 0.8) +  
  geom_point(aes(x = app, y = int_peak), 
             size = 2, color = 'grey25') + 
  geom_line(aes(x = app, y = int_peak), 
            size = 0.75, color = 'grey25') +     
  stat_pvalue_manual(data = peak_stat_df_oreol, label = 'p.adj.signif',
                     hide.ns = TRUE, size = 8) +
  labs(title = "HPCA_mBaoJin_app_peaks_at_diff_app_time", x = "App_time", y = expression(Delta * F / F[0])) +
  theme_minimal()+
  theme(legend.position = "none", text=element_text(size = font.size, family = font.fam, face = 'bold'))

#points
ggplot(data = peaks %>% filter(lab_id != 'psd') %>%
         mutate(app = as.numeric(as.character(app))),  
       aes(x = app, y = int_peak,
           color = lab_id, group = lab_id)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(fun = median, geom = 'line', size = 1.25) +
  stat_summary(fun = median, geom = 'point', size = 2.5) +
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median,
               geom = 'errorbar', width = .1, size = 0.75) +
  scale_color_manual(values = c('oreol' = 'yellow3', 'shaft' = 'blue3')) +
  scale_x_log10(breaks = c(0.5, 2.5, 5, 10, 20, 30, 60, 120)) +
  theme_classic() +
  theme(legend.position = 'right',
        text = element_text(size = font.size, family = font.fam),
        plot.caption = element_text(size = font.size - 4),
        legend.title = element_text(size = font.size - 2)) +
  labs(color = 'ROItype',
       x = 'Application duration',
       y = expression(Delta * F / F[0]))

#DECAY

decay_0.5sec <- all_cells %>% filter(app == "0.5") %>%
  filter(time %in% seq(31, 41)) %>%
  group_by(roi, lab_id) %>%
  select(time, mean_df, lab_id, roi) %>%
  mutate(idx = row_number()) %>%
  distinct() %>%
  droplevels()


decay_0.5_shaft <- lm(log(mean_df) ~ idx, data = decay_0.5sec %>% filter(lab_id == "shaft"))

#SETTINGS
#change name of raw 
names(df_)[names(df_) == "n"] <- "n1"
#change meanings in row 
df$app[df$app == 10] <- 20
#filter_something
raw_shaft<- raw_shaft %>% 
  filter(!roi %in% c())