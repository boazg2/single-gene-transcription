#all for 10 kb
library(zoo)
library(gridExtra)
setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

k_weak = 0.008
k_mod = 0.05
k_strong = 0.2

opt_gyr = 5.625e-05
opt_top = 0.000225

weak_df = import_data("weak-mod/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_df = import_data("mod-mod/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_df = import_data("strong-mod/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_close_df = import_data("weak-close/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_close_df = import_data("mod-close/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_close_df = import_data("strong-close/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_low_df = import_data("weak-close_low-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_low_df = import_data("mod-close_low-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_low_df = import_data("strong-close_low-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_far_df = import_data("weak-far/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_far_df = import_data("mod-far/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_far_df = import_data("strong-far/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_far_hi_df = import_data("weak-far_high-topo/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_far_hi_df = import_data("mod-far_high-topo/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_far_hi_df = import_data("strong-far_high-topo/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_mod_low_df = import_data("weak-mod_low-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_mod_low_df = import_data("mod-mod_low-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_mod_low_df = import_data("strong-mod_low-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)
# 2024-12-09
# Calculate elongation speeds

weak_df = mutate(weak_df, speed = 1000 / mean_elong)
weak_df$speed[is.nan(weak_df$speed)] <- 0
mod_df = mutate(mod_df, speed = 1000 / mean_elong)
mod_df$speed[is.nan(mod_df$speed)] <- 0
strong_df = mutate(strong_df, speed = 1000 / mean_elong)
strong_df$speed[is.nan(strong_df$speed)] <- 0

weak_close_df = mutate(weak_close_df, speed = 1000 / mean_elong)
weak_close_df$speed[is.nan(weak_close_df$speed)] <- 0
mod_close_df = mutate(mod_close_df, speed = 1000 / mean_elong)
mod_close_df$speed[is.nan(mod_close_df$speed)] <- 0
strong_close_df = mutate(strong_close_df, speed = 1000 / mean_elong)
strong_close_df$speed[is.nan(strong_close_df$speed)] <- 0

weak_low_df = mutate(weak_low_df, speed = 1000 / mean_elong)
weak_low_df$speed[is.nan(weak_low_df$speed)] <- 0
mod_low_df = mutate(mod_low_df, speed = 1000 / mean_elong)
mod_low_df$speed[is.nan(mod_low_df$speed)] <- 0
strong_low_df = mutate(strong_low_df, speed = 1000 / mean_elong)
strong_low_df$speed[is.nan(strong_low_df$speed)] <- 0



#rate
get_transcription_grid(weak_df, 
                                   plot_title = "weak 10 kb", 
                                   fillvar = "normed_rate",
                                   x_label = "Gyrase Activity (Lk/kb/s)",
                                   y_label = "TopoI Activity (Lk/kb/s)",
                                   fill_label = "Normed Rate",
                                   gyr_scale = 4e3,
                                   topoI_scale = 1e3,
                                   fill_min = 0,
                                   fill_max = 1,
                                   opt_x = opt_gyr,
                                   opt_y = opt_top) +
  theme_classic()

get_transcription_grid(mod_df, 
                       plot_title = "moderate 10 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()

get_transcription_grid(strong_df, 
                       plot_title = "strong 10 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()


#elongation
get_transcription_grid(weak_df, 
                       plot_title = "weak 10 kb", 
                       fillvar = "mean_elong",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 40,
                       fill_max = 310,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       low_color = "white",
                       high_color = "red") +
  theme_classic()

get_transcription_grid(mod_df, 
                       plot_title = "moderate 10 kb", 
                       fillvar = "mean_elong",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 40,
                       fill_max = 310,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       low_color = "white",
                       high_color = "red") +
  theme_classic()

get_transcription_grid(strong_df, 
                       plot_title = "strong 10 kb", 
                       fillvar = "mean_elong",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 40,
                       fill_max = 310,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       low_color = "white",
                       high_color = "red") +
  theme_classic()

# 2024-12-09 speed
get_transcription_grid(weak_df, 
                       plot_title = "weak 10 kb", 
                       fillvar = "speed",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Elongation Speed (bp/s)",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 25,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       low_color = "black",
                       high_color = "white") +
  theme_classic()

get_transcription_grid(mod_df, 
                       plot_title = "moderate 10 kb", 
                       fillvar = "speed",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Elongation Speed (bp/s)",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 25,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       low_color = "black",
                       high_color = "white") +
  theme_classic()

get_transcription_grid(strong_df, 
                       plot_title = "strong 10 kb", 
                       fillvar = "speed",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Elongation Speed (bp/s)",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 25,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       low_color = "black",
                       high_color = "white") +
  theme_classic()



#transcription, close
get_transcription_grid(weak_close_df, 
                       plot_title = "weak 1 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10) +
  theme_classic()

get_transcription_grid(mod_close_df, 
                       plot_title = "moderate 1 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10) +
  theme_classic()

get_transcription_grid(strong_close_df, 
                       plot_title = "strong 1 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10) +
  theme_classic()


#transcription, close, low topo
get_transcription_grid(weak_low_df, 
                       plot_title = "weak 1 kb low topo", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       opt = F) +
  theme_classic()

get_transcription_grid(mod_low_df, 
                       plot_title = "moderate 1 kb low topo", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       opt = F) +
  theme_classic()

get_transcription_grid(strong_low_df, 
                       plot_title = "strong 1 kb low topo", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       opt = F) +
  theme_classic()

#2024-09-03: plotting steric clash
weak_mod_steric = import_steric_data("weak-mod/data", start_time = 1800, ki = k_weak)
mod_mod_steric = import_steric_data("mod-mod/data", start_time = 1800, ki = k_mod)
strong_mod_steric = import_steric_data("strong-mod/data", start_time = 1800, ki = k_strong)


get_transcription_grid(weak_mod_steric,
                       plot_title = "weak 10 kb steric",
                       fillvar = "free_prom",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Free Promoter",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()

get_transcription_grid(mod_mod_steric,
                       plot_title = "mod 10 kb steric",
                       fillvar = "free_prom",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Free Promoter",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()

get_transcription_grid(strong_mod_steric,
                       plot_title = "strong 10 kb steric",
                       fillvar = "free_prom",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Free Promoter",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()

# Box plots of drug treatment

lansG_list = unique(weak_df$lansG)
lansG_list_low = unique(weak_mod_low_df$lansG)

lansT_list = unique(weak_df$lansT)
lansT_list_low = unique(weak_mod_low_df$lansT)

lansG1 = lansG_list[6]
lansG2 = lansG_list[2]
lansG1_low = lansG_list_low[6]
lansG2_low = lansG_list_low[2]

lansT1 = lansT_list[6]
lansT2 = lansT_list[2]
lansT1_low = lansT_list_low[6]
lansT2_low = lansT_list_low[2]


change_box_df = data.frame(
  promoter = rep(c('Weak', 'Moderate', 'Strong'), 4),
  topo_level = rep(c('Low', 'Low', 'Low', 'Normal', 'Normal', 'Normal'), 2),
  inhibited = c('Gyrase','Gyrase','Gyrase','Gyrase','Gyrase','Gyrase','TopoI','TopoI','TopoI','TopoI','TopoI','TopoI'),
  log2_change = log2(c(
    filter(weak_mod_low_df, lansT == lansT1_low, lansG == lansG2_low)$mean_rate/filter(weak_mod_low_df, lansT == lansT1_low, lansG == lansG1_low)$mean_rate,
    filter(mod_mod_low_df, lansT == lansT1_low, lansG == lansG2_low)$mean_rate/filter(mod_mod_low_df, lansT == lansT1_low, lansG == lansG1_low)$mean_rate,
    filter(strong_mod_low_df, lansT == lansT1_low, lansG == lansG2_low)$mean_rate/filter(strong_mod_low_df, lansT == lansT1_low, lansG == lansG1_low)$mean_rate,
    filter(weak_df, lansT == lansT1, lansG == lansG2)$mean_rate/filter(weak_df, lansT == lansT1, lansG == lansG1)$mean_rate,
    filter(mod_df, lansT == lansT1, lansG == lansG2)$mean_rate/filter(mod_df, lansT == lansT1, lansG == lansG1)$mean_rate,
    filter(strong_df, lansT == lansT1, lansG == lansG2)$mean_rate/filter(strong_df, lansT == lansT1, lansG == lansG1)$mean_rate,
    filter(weak_mod_low_df, lansT == lansT2_low, lansG == lansG1_low)$mean_rate/filter(weak_mod_low_df, lansT == lansT1_low, lansG == lansG1_low)$mean_rate,
    filter(mod_mod_low_df, lansT == lansT2_low, lansG == lansG1_low)$mean_rate/filter(mod_mod_low_df, lansT == lansT1_low, lansG == lansG1_low)$mean_rate,
    filter(strong_mod_low_df, lansT == lansT2_low, lansG == lansG1_low)$mean_rate/filter(strong_mod_low_df, lansT == lansT1_low, lansG == lansG1_low)$mean_rate,
    filter(weak_df, lansT == lansT2, lansG == lansG1)$mean_rate/filter(weak_df, lansT == lansT1, lansG == lansG1)$mean_rate,
    filter(mod_df, lansT == lansT2, lansG == lansG1)$mean_rate/filter(mod_df, lansT == lansT1, lansG == lansG1)$mean_rate,
    filter(strong_df, lansT == lansT2, lansG == lansG1)$mean_rate/filter(strong_df, lansT == lansT1, lansG == lansG1)$mean_rate
  ))
)

change_box_df$promoter = factor(change_box_df$promoter, level = c("Weak", "Moderate", "Strong"))

ggplot(data = filter(change_box_df, inhibited == 'TopoI', topo_level == "Low"), aes(x = promoter, y = log2_change)) + 
  stat_summary(fun = mean, geom = "bar", fill = "blue") +
  labs(x = "Promoter Strength", y = "log2 Change", title = "Seconeolitsine, Low Topoisomerases") +
  theme_classic() +
  ylim(-3, 2)

ggplot(data = filter(change_box_df, inhibited == 'TopoI', topo_level == "Normal"), aes(x = promoter, y = log2_change)) + 
  stat_summary(fun = mean, geom = "bar", fill = "blue") +
  labs(x = "Promoter Strength", y = "log2 Change", title = "Seconeolitsine") +
  theme_classic() +
  ylim(-3, 2)

ggplot(data = filter(change_box_df, inhibited == 'Gyrase', topo_level == "Low"), aes(x = promoter, y = log2_change)) + 
  stat_summary(fun = mean, geom = "bar", fill = "blue") +
  labs(x = "Promoter Strength", y = "log2 Change", title = "Novobiocin, Low Topoisomerases") +
  theme_classic() +
  ylim(-3, 2)

ggplot(data = filter(change_box_df, inhibited == 'Gyrase', topo_level == "Normal"), aes(x = promoter, y = log2_change)) + 
  stat_summary(fun = mean, geom = "bar", fill = "blue") +
  labs(x = "Promoter Strength", y = "log2 Change", title = "Novobiocin") +
  theme_classic() +
  ylim(-3, 2)





#2024-12-09
#log2 fold changes resulting from reducing barrier distances
weak_low_reduction <- weak_df %>%
  inner_join(weak_low_df, by = c("lansG", "lansT"), suffix = c("_mod", "_close_low")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    mod_rate = mean_rate_mod,
    close_low_rate = mean_rate_close_low,
    reduction = log2(close_low_rate / mod_rate)
  )
mod_low_reduction <- mod_df %>%
  inner_join(mod_low_df, by = c("lansG", "lansT"), suffix = c("_mod", "_close_low")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    mod_rate = mean_rate_mod,
    close_low_rate = mean_rate_close_low,
    reduction = log2(close_low_rate / mod_rate)
  )
strong_low_reduction <- strong_df %>%
  inner_join(strong_low_df, by = c("lansG", "lansT"), suffix = c("_mod", "_close_low")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    mod_rate = mean_rate_mod,
    close_low_rate = mean_rate_close_low,
    reduction = log2(close_low_rate / mod_rate)
  )

get_transcription_grid_tricolor(weak_low_reduction, 
                       plot_title = "weak 10 kb to weak 1 kb low topo", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()
get_transcription_grid_tricolor(mod_low_reduction, 
                       plot_title = "mod 10 kb to mod 1 kb low topo", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()
get_transcription_grid_tricolor(strong_low_reduction, 
                       plot_title = "strong 10 kb to strong 1 kb low topo", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()



# for the 10x topo case (e.g., weak_mod to weak_close, rather than weak_close_low-topo)
# Sort both dataframes by lansT and lansG
weak_df_sorted <- weak_df %>% arrange(lansT, lansG)
weak_close_df_sorted <- weak_close_df %>% arrange(lansT, lansG)
mod_df_sorted <- mod_df %>% arrange(lansT, lansG)
mod_close_df_sorted <- mod_close_df %>% arrange(lansT, lansG)
strong_df_sorted <- strong_df %>% arrange(lansT, lansG)
strong_close_df_sorted <- strong_close_df %>% arrange(lansT, lansG)

weak_reduction <- bind_cols(
  weak_df_sorted %>% ungroup() %>% select(lansG, lansT, mod_rate = mean_rate),
  weak_close_df_sorted %>% ungroup() %>% select(close_rate = mean_rate)
) %>%
  mutate(
    reduction = log2(close_rate / mod_rate)
  )
mod_reduction <- bind_cols(
  mod_df_sorted %>% ungroup() %>% select(lansG, lansT, mod_rate = mean_rate),
  mod_close_df_sorted %>% ungroup() %>% select(close_rate = mean_rate)
) %>%
  mutate(
    reduction = log2(close_rate / mod_rate)
  )
strong_reduction <- bind_cols(
  strong_df_sorted %>% ungroup() %>% select(lansG, lansT, mod_rate = mean_rate),
  strong_close_df_sorted %>% ungroup() %>% select(close_rate = mean_rate)
) %>%
  mutate(
    reduction = log2(close_rate / mod_rate)
  )



get_transcription_grid_tricolor(weak_reduction, 
                       plot_title = "weak 10 kb to weak 1 kb", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()
get_transcription_grid_tricolor(mod_reduction, 
                       plot_title = "mod 10 kb to mod 1 kb", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()
get_transcription_grid_tricolor(strong_reduction, 
                       plot_title = "strong 10 kb to strong 1 kb", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr,
                       opt_y = opt_top) +
  theme_classic()



# repeat for the moderate to far barrier transition:

far_weak_hi_increase <- weak_df %>%
  inner_join(weak_far_hi_df, by = c("lansG", "lansT"), suffix = c("_mod", "_far_hi")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    mod_rate = mean_rate_mod,
    far_hi_rate = mean_rate_far_hi,
    reduction = log2(far_hi_rate / mod_rate)
  )

far_mod_hi_increase <- mod_df %>%
  inner_join(mod_far_hi_df, by = c("lansG", "lansT"), suffix = c("_mod", "_far_hi")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    mod_rate = mean_rate_mod,
    far_hi_rate = mean_rate_far_hi,
    reduction = log2(far_hi_rate / mod_rate)
  )
far_strong_hi_increase <- strong_df %>%
  inner_join(strong_far_hi_df, by = c("lansG", "lansT"), suffix = c("_mod", "_far_hi")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    mod_rate = mean_rate_mod,
    far_hi_rate = mean_rate_far_hi,
    reduction = log2(far_hi_rate / mod_rate)
  )

get_transcription_grid_tricolor(far_weak_hi_increase, 
                       plot_title = "weak 10 kb to weak 100 kb high topo", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()
get_transcription_grid_tricolor(far_mod_hi_increase, 
                       plot_title = "mod 10 kb to mod 100 kb high topo", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()
get_transcription_grid_tricolor(far_strong_hi_increase, 
                       plot_title = "strong 10 kb to strong 100 kb high topo", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()


weak_far_df_sorted <- weak_far_df %>% arrange(lansT, lansG)
mod_far_df_sorted <- mod_far_df %>% arrange(lansT, lansG)
strong_far_df_sorted <- strong_far_df %>% arrange(lansT, lansG)
#note that mod barrier cases have already been done

far_weak_reduction <- bind_cols(
  weak_far_df_sorted %>% ungroup() %>% select(lansG, lansT, mod_rate = mean_rate),
  weak_df_sorted %>% ungroup() %>% select(close_rate = mean_rate)
) %>%
  mutate(
    reduction = log2(mod_rate / close_rate)
  )
far_mod_reduction <- bind_cols(
  mod_far_df_sorted %>% ungroup() %>% select(lansG, lansT, mod_rate = mean_rate),
  mod_df_sorted %>% ungroup() %>% select(close_rate = mean_rate)
) %>%
  mutate(
    reduction = log2(mod_rate / close_rate)
  )
far_strong_reduction <- bind_cols(
  strong_far_df_sorted %>% ungroup() %>% select(lansG, lansT, mod_rate = mean_rate),
  strong_df_sorted %>% ungroup() %>% select(close_rate = mean_rate)
) %>%
  mutate(
    reduction = log2(mod_rate / close_rate)
  )

get_transcription_grid_tricolor(far_weak_reduction, 
                       plot_title = "weak 10 kb to weak 100 kb", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()
get_transcription_grid_tricolor(far_mod_reduction, 
                       plot_title = "mod 10 kb to mod 100 kb", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()
get_transcription_grid_tricolor(far_strong_reduction, 
                       plot_title = "strong 10 kb to strong 100 kb", 
                       fillvar = "reduction",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = -3,
                       fill_max = 3,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()



