#all for 10 kb

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
                       plot_title = "strong 10 ksb", 
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


# 2024-08-09: for supplemental figure on effects of the further barrier distance
weak_far_df = import_data("weak-far/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_far_df = import_data("mod-far/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_far_df = import_data("strong-far/data", start_time = 3600) %>%
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

# far
get_transcription_grid(weak_far_df, 
                       plot_title = "weak 100 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()

get_transcription_grid(mod_far_df, 
                       plot_title = "moderate 100 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()

get_transcription_grid(strong_far_df, 
                       plot_title = "strong 100 kb", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr/10,
                       opt_y = opt_top/10) +
  theme_classic()

# mod low-topo

get_transcription_grid(weak_mod_low_df, 
                       plot_title = "weak 10 kb low topo", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       opt = F) +
  theme_classic()

get_transcription_grid(mod_mod_low_df, 
                       plot_title = "moderate 10 kb low topo", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       opt = F) +
  theme_classic()

get_transcription_grid(strong_mod_low_df, 
                       plot_title = "strong 10 kb low topo", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Normed Rate",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top,
                       opt = F) +
  theme_classic()
