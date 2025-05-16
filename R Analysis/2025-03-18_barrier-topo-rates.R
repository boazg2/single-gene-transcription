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

weak_far_high_df = import_data("weak-far_high-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_far_high_df = import_data("mod-far_high-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_far_high_df = import_data("strong-far_high-topo/data", start_time = 1800) %>%
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

weak_mod_high_df = import_data("weak-mod_high-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_mod_high_df = import_data("mod-mod_high-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_mod_high_df = import_data("strong-mod_high-topo/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)


#log2 fold changes resulting from reducing barrier distances
weak_mod_to_weak_close_low <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_low_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_mod_to_mod_close_low <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_low_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_mod_to_strong_close_low  <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_low_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))


weak_mod_to_weak_close <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_close_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
  ) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_mod_to_mod_close <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_close_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_mod_to_strong_close <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_close_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

weak_mod_to_weak_far_high <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_far_high_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_mod_to_mod_far_high <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_far_high_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_mod_to_strong_far_high <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_far_high_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

weak_mod_to_weak_far <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_far_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_mod_to_mod_far <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_far_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_mod_to_strong_far <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_far_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))


weak_mod_to_weak_mod_low <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_mod_low_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_mod_to_mod_mod_low <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_mod_low_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_mod_to_strong_mod_low  <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_mod_low_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))


weak_mod_to_weak_mod_high <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_mod_high_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_mod_to_mod_mod_high <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_mod_high_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_mod_to_strong_mod_high  <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_mod_high_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

# plot grids


get_transcription_grid_tricolor(weak_mod_to_weak_close_low, 
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

get_transcription_grid_tricolor(mod_mod_to_mod_close_low, 
                                plot_title = "mod 10 kb to weak 1 kb low topo", 
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

get_transcription_grid_tricolor(strong_mod_to_strong_close_low, 
                                plot_title = "mod 10 kb to weak 1 kb low topo", 
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


get_transcription_grid_tricolor(weak_mod_to_weak_close, 
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

get_transcription_grid_tricolor(mod_mod_to_mod_close, 
                                plot_title = "mod 10 kb to weak 1 kb", 
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

get_transcription_grid_tricolor(strong_mod_to_strong_close, 
                                plot_title = "mod 10 kb to weak 1 kb", 
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


get_transcription_grid_tricolor(weak_mod_to_weak_far_high, 
                                plot_title = "weak 10 kb to weak 100 kb high topo", 
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

get_transcription_grid_tricolor(mod_mod_to_mod_far_high, 
                                plot_title = "mod 10 kb to mod 100 kb high topo", 
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

get_transcription_grid_tricolor(strong_mod_to_strong_far_high, 
                                plot_title = "strong 10 kb to strong 100 kb high topo", 
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


get_transcription_grid_tricolor(weak_mod_to_weak_far, 
                                plot_title = "weak 10 kb to weak 100 kb", 
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

get_transcription_grid_tricolor(mod_mod_to_mod_far, 
                                plot_title = "mod 10 kb to mod 100 kb", 
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

get_transcription_grid_tricolor(strong_mod_to_strong_far, 
                                plot_title = "strong 10 kb to strong 100 kb", 
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


# mod to mod low

get_transcription_grid_tricolor(weak_mod_to_weak_mod_low, 
                                plot_title = "weak 10 kb to weak 10 kb low topo", 
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

get_transcription_grid_tricolor(mod_mod_to_mod_mod_low, 
                                plot_title = "mod 10 kb to mod 10 kb low topo", 
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

get_transcription_grid_tricolor(strong_mod_to_strong_mod_low, 
                                plot_title = "strong 10 kb to strong 10 kb low topo", 
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

# mod to mod high

get_transcription_grid_tricolor(weak_mod_to_weak_mod_high, 
                                plot_title = "weak 10 kb to weak 10 kb high topo", 
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

get_transcription_grid_tricolor(mod_mod_to_mod_mod_high, 
                                plot_title = "mod 10 kb to mod 10 kb high topo", 
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

get_transcription_grid_tricolor(strong_mod_to_strong_mod_high, 
                                plot_title = "strong 10 kb to strong 10 kb high topo", 
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


# Line plots
# scale to translate into kilobases and account for gyrase processivity
top_scale = 1e3
gyr_scale = 4e3
topo_slice = 0.2380952

# Define promoter strength labels
weak_mod_to_weak_close_low$Promoter <- "Weak"
mod_mod_to_mod_close_low$Promoter <- "Moderate"
strong_mod_to_strong_close_low$Promoter <- "Strong"

data_close_low <- rbind(weak_mod_to_weak_close_low, mod_mod_to_mod_close_low, strong_mod_to_strong_close_low) %>%
  filter(lansT != 0, lansG != 0) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)

weak_mod_to_weak_close$Promoter <- "Weak"
mod_mod_to_mod_close$Promoter <- "Moderate"
strong_mod_to_strong_close$Promoter <- "Strong"

data_close <- rbind(weak_mod_to_weak_close, mod_mod_to_mod_close, strong_mod_to_strong_close) %>%
  filter(lansT != 0, lansG != 0) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)

weak_mod_to_weak_far_high$Promoter <- "Weak"
mod_mod_to_mod_far_high$Promoter <- "Moderate"
strong_mod_to_strong_far_high$Promoter <- "Strong"

data_far_high <- rbind(weak_mod_to_weak_far_high, mod_mod_to_mod_far_high, strong_mod_to_strong_far_high) %>%
  filter(lansT != 0, lansG != 0) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)

weak_mod_to_weak_far$Promoter <- "Weak"
mod_mod_to_mod_far$Promoter <- "Moderate"
strong_mod_to_strong_far$Promoter <- "Strong"

data_far <- rbind(weak_mod_to_weak_far, mod_mod_to_mod_far, strong_mod_to_strong_far) %>%
  filter(lansT != 0, lansG != 0) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)

weak_mod_to_weak_mod_low$Promoter <- "Weak"
mod_mod_to_mod_mod_low$Promoter <- "Moderate"
strong_mod_to_strong_mod_low$Promoter <- "Strong"

data_mod_low <- rbind(weak_mod_to_weak_mod_low, mod_mod_to_mod_mod_low, strong_mod_to_strong_mod_low) %>%
  filter(lansT != 0, lansG != 0) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)


weak_mod_to_weak_mod_high$Promoter <- "Weak"
mod_mod_to_mod_mod_high$Promoter <- "Moderate"
strong_mod_to_strong_mod_high$Promoter <- "Strong"

data_mod_high <- rbind(weak_mod_to_weak_mod_high, mod_mod_to_mod_mod_high, strong_mod_to_strong_mod_high) %>%
  filter(lansT != 0, lansG != 0) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)

# Create plots with adjusted font sizes
custom_theme <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),
  axis.title = element_text(size = 0.125 * 72),
  legend.position = "none"
)

plot_list <- list(
  ggplot(subset(data_close_low, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    ylab("10 kb to 1 kb\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5) +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_close_low, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_low, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    ylab("1/10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_low, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_close, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    ylab("10 kb to 1 kb, 10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_close, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_far_high, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    ylab("10 kb to 100 kb\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_far_high, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_high, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    ylab("10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_high, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_far, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    xlab("TopoI Activity (Lk/kb/s)") +
    ylab("10 kb to 100 kb, 1/10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme,
  
  ggplot(subset(data_far, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point() +
    xlab("Gyrase Activity (Lk/kb/s)") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(legend.position = c(1, 1.1), 
          legend.justification = c("right", "top"), 
          legend.direction = "horizontal", 
          legend.title.position = "top",
          legend.title = element_text(size = 0.125 * 72),
          legend.text = element_text(size = 0.105 * 72))
  
  
)

# Arrange in a 5x2 grid
combined_plot <- arrangeGrob(grobs = plot_list, ncol = 2, nrow = 6)

ggsave("figures/barrier_plots.png", combined_plot, dpi = 300, width = 5.2, height = 10, units = "in")
ggsave("figures/barrier_plots.eps", combined_plot, dpi = 300, width = 5.2, height = 10, units = "in")

# Create plots with adjusted font sizes
custom_theme <- theme(
  axis.text = element_text(size = 0.262 * 72, color = "black"),
  axis.title = element_text(size = 0.262 * 72),
  legend.position = "none"
)

point_size = 3.5

plot_list <- list(
  ggplot(subset(data_close_low, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    ylab("10 kb to 1 kb\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5) +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_close_low, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_low, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    ylab("1/10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_low, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_close, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    ylab("10 kb to 1 kb, 10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_close, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_far_high, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    ylab("10 kb to 100 kb\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_far_high, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_high, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    ylab("10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_mod_high, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()),
  
  ggplot(subset(data_far, gyr_activity == topo_slice), aes(x = top_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    xlab("TopoI Activity (Lk/kb/s)") +
    ylab("10 kb to 100 kb, 1/10x Topo\nLog2 Fold Change") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme,
  
  ggplot(subset(data_far, top_activity == topo_slice), aes(x = gyr_activity, y = reduction, color = Promoter)) +
    geom_line() + geom_point(size = point_size) +
    xlab("Gyrase Activity (Lk/kb/s)") +
    scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
    ylim(-2.6, 2.6) +
    xlim(0, 0.5)  +
    theme_classic() + custom_theme +
    theme(axis.title.y = element_blank()) +
    theme(legend.position = c(1, 1.1), 
          legend.justification = c("right", "top"), 
          legend.direction = "horizontal", 
          legend.title.position = "top",
          legend.title = element_text(size = 0.295 * 72),
          legend.text = element_text(size = 0.262 * 72))
  
  
)
combined_plot <- arrangeGrob(grobs = plot_list, ncol = 2, nrow = 6)

ggsave("figures/barrier_plots_burf.png", combined_plot, dpi = 300, width = 12.667, height = 16.4, units = "in")
ggsave("figures/barrier_plots_burf.eps", combined_plot, dpi = 300, width = 12.667, height = 16.4, units = "in")
