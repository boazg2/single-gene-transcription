
# far high analysis

weak_far_hi_df = import_data("weak-far_high-topo/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_far_hi_df = import_data("mod-far_high-topo/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_far_hi_df = import_data("strong-far_high-topo/data", start_time = 3600) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

#rate
get_transcription_grid(weak_far_hi_df, 
                       plot_title = "weak 100 kb high topo", 
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

get_transcription_grid(mod_far_hi_df, 
                       plot_title = "moderate 100 kb high topo", 
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

get_transcription_grid(strong_far_hi_df, 
                       plot_title = "strong 100 kb high topo", 
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

weak_far_hi_df = mutate(weak_far_hi_df, speed = 1000 / mean_elong)
weak_far_hi_df$speed[is.nan(weak_far_hi_df$speed)] <- 0
mod_far_hi_df = mutate(mod_far_hi_df, speed = 1000 / mean_elong)
mod_far_hi_df$speed[is.nan(mod_far_hi_df$speed)] <- 0
strong_far_hi_df = mutate(strong_far_hi_df, speed = 1000 / mean_elong)
strong_far_hi_df$speed[is.nan(strong_far_hi_df$speed)] <- 0

get_transcription_grid(weak_far_hi_df, 
                       plot_title = "weak 100 kb high topo", 
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

get_transcription_grid(mod_far_hi_df, 
                       plot_title = "moderate 100 kb high topo", 
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

get_transcription_grid(strong_far_hi_df, 
                       plot_title = "strong 100 kb high topo", 
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

#promoter availability
weak_far_steric = import_steric_data("weak-far_high-topo/data", start_time = 1800, ki = k_weak)
mod_far_steric = import_steric_data("mod-far_high-topo/data", start_time = 1800, ki = k_mod)
strong_far_steric = import_steric_data("strong-far_high-topo/data", start_time = 1800, ki = k_strong)


get_transcription_grid(weak_far_steric,
                       plot_title = "weak 100 kb high topo steric",
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

get_transcription_grid(mod_far_steric,
                       plot_title = "mod 100 kb high topo steric",
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

get_transcription_grid(strong_far_steric,
                       plot_title = "strong 100 kb high topo steric",
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

