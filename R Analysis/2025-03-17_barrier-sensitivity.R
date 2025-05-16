library(zoo)
library(gridExtra)
setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2025-03-14_barrier-sensitivity")

k_weak = 0.008
k_mod = 0.05
k_strong = 0.2

opt_gyr = 0.0023809525
opt_top = 0.007440475


# data import
weak_base = import_data("weak_up-320_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_base = import_data("mod_up-320_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_base = import_data("strong_up-320_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_up1044 = import_data("weak_up-1044_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_up1044 = import_data("mod_up-1044_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_up1044 = import_data("strong_up-1044_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_up3408 = import_data("weak_up-3408_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_up3408 = import_data("mod_up-3408_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_up3408 = import_data("strong_up-3408_down-250", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_down900= import_data("weak_up-320_down-900", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_down900 = import_data("mod_up-320_down-900", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_down900 = import_data("strong_up-320_down-900", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_down3200= import_data("weak_up-320_down-3200", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_down3200 = import_data("mod_up-320_down-3200", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_down3200 = import_data("strong_up-320_down-3200", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)


# plot base case

weak_base_rate <- get_transcription_grid(weak_base, 
                       plot_title = "weak 320 up 250 down", 
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

#ggsave("plots/weak_base_rate.eps", weak_base_rate)
#ggsave("plots/weak_base_rate.png", weak_base_rate, dpi=300)


mod_base_rate <- get_transcription_grid(mod_base, 
                       plot_title = "mod 320 up 250 down", 
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

#ggsave("plots/mod_base_rate.eps", mod_base_rate)
#ggsave("plots/mod_base_rate.png", mod_base_rate, dpi=300)


strong_base_rate <- get_transcription_grid(strong_base, 
                       plot_title = "strong 320 up 250 down", 
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

#ggsave("plots/strong_base_rate.eps", strong_base_rate)
#ggsave("plots/strong_base_rate.png", strong_base_rate, dpi=300)



# calculate fold changes

weak_up1044_change <- weak_base %>%
  inner_join(weak_up1044, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

mod_up1044_change <- mod_base %>%
  inner_join(mod_up1044, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

strong_up1044_change <- strong_base %>%
  inner_join(strong_up1044, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

weak_up3408_change <- weak_base %>%
  inner_join(weak_up3408, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

mod_up3408_change <- mod_base %>%
  inner_join(mod_up3408, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

strong_up3408_change <- strong_base %>%
  inner_join(strong_up3408, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

weak_down900_change <- weak_base %>%
  inner_join(weak_down900, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

mod_down900_change <- mod_base %>%
  inner_join(mod_down900, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

strong_down900_change <- strong_base %>%
  inner_join(strong_down900, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

weak_down3200_change <- weak_base %>%
  inner_join(weak_down3200, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

mod_down3200_change <- mod_base %>%
  inner_join(mod_down3200, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

strong_down3200_change <- strong_base %>%
  inner_join(strong_down3200, by = c("lansG", "lansT"), suffix = c("_base", "_new")) %>%
  transmute(
    lansG = lansG,
    lansT = lansT,
    base_rate = mean_rate_base,
    new_rate = mean_rate_new,
    change = log2(new_rate / base_rate)
  )

# Plot fold change dataframes

ggplot_weak_up1044 <- get_transcription_grid_tricolor(weak_up1044_change, 
                                               plot_title = "Weak Up1044 Fold Change", 
                                               fillvar = "change",
                                               x_label = "Gyrase Activity (Lk/kb/s)",
                                               y_label = "TopoI Activity (Lk/kb/s)",
                                               fill_label = "Log2 Fold Change",
                                               gyr_scale = 4e3,
                                               topoI_scale = 1e3,
                                               fill_min = -2,
                                               fill_max = 2,
                                               opt_x = opt_gyr,
                                               opt_y = opt_top) +
  theme_classic()
#ggsave("plots/weak_up1044_change.eps", ggplot_weak_up1044)
#ggsave("plots/weak_up1044_change.png", ggplot_weak_up1044, dpi = 300)

ggplot_mod_up1044 <- get_transcription_grid_tricolor(mod_up1044_change, 
                                              plot_title = "Moderate Up1044 Fold Change", 
                                              fillvar = "change",
                                              x_label = "Gyrase Activity (Lk/kb/s)",
                                              y_label = "TopoI Activity (Lk/kb/s)",
                                              fill_label = "Log2 Fold Change",
                                              gyr_scale = 4e3,
                                              topoI_scale = 1e3,
                                              fill_min = -2,
                                              fill_max = 2,
                                              opt_x = opt_gyr,
                                              opt_y = opt_top) +
  theme_classic()
#ggsave("plots/mod_up1044_change.eps", ggplot_mod_up1044)
#ggsave("plots/mod_up1044_change.png", ggplot_mod_up1044, dpi = 300)

ggplot_strong_up1044 <- get_transcription_grid_tricolor(strong_up1044_change, 
                                                 plot_title = "Strong Up1044 Fold Change", 
                                                 fillvar = "change",
                                                 x_label = "Gyrase Activity (Lk/kb/s)",
                                                 y_label = "TopoI Activity (Lk/kb/s)",
                                                 fill_label = "Log2 Fold Change",
                                                 gyr_scale = 4e3,
                                                 topoI_scale = 1e3,
                                                 fill_min = -2,
                                                 fill_max = 2,
                                                 opt_x = opt_gyr,
                                                 opt_y = opt_top) +
  theme_classic()
#ggsave("plots/strong_up1044_change.eps", ggplot_strong_up1044)
#ggsave("plots/strong_up1044_change.png", ggplot_strong_up1044, dpi = 300)

ggplot_weak_up3408 <- get_transcription_grid_tricolor(weak_up3408_change, 
                                               plot_title = "Weak Up3408 Fold Change", 
                                               fillvar = "change",
                                               x_label = "Gyrase Activity (Lk/kb/s)",
                                               y_label = "TopoI Activity (Lk/kb/s)",
                                               fill_label = "Log2 Fold Change",
                                               gyr_scale = 4e3,
                                               topoI_scale = 1e3,
                                               fill_min = -2,
                                               fill_max = 2,
                                               opt_x = opt_gyr,
                                               opt_y = opt_top) +
  theme_classic()
#ggsave("plots/weak_up3408_change.eps", ggplot_weak_up3408)
#ggsave("plots/weak_up3408_change.png", ggplot_weak_up3408, dpi = 300)

ggplot_mod_up3408 <- get_transcription_grid_tricolor(mod_up3408_change, 
                                              plot_title = "Moderate Up3408 Fold Change", 
                                              fillvar = "change",
                                              x_label = "Gyrase Activity (Lk/kb/s)",
                                              y_label = "TopoI Activity (Lk/kb/s)",
                                              fill_label = "Log2 Fold Change",
                                              gyr_scale = 4e3,
                                              topoI_scale = 1e3,
                                              fill_min = -2,
                                              fill_max = 2,
                                              opt_x = opt_gyr,
                                              opt_y = opt_top) +
  theme_classic()
#ggsave("plots/mod_up3408_change.eps", ggplot_mod_up3408)
#ggsave("plots/mod_up3408_change.png", ggplot_mod_up3408, dpi = 300)

ggplot_strong_up3408 <- get_transcription_grid_tricolor(strong_up3408_change, 
                                                 plot_title = "Strong Up3408 Fold Change", 
                                                 fillvar = "change",
                                                 x_label = "Gyrase Activity (Lk/kb/s)",
                                                 y_label = "TopoI Activity (Lk/kb/s)",
                                                 fill_label = "Log2 Fold Change",
                                                 gyr_scale = 4e3,
                                                 topoI_scale = 1e3,
                                                 fill_min = -2,
                                                 fill_max = 2,
                                                 opt_x = opt_gyr,
                                                 opt_y = opt_top) +
  theme_classic()
#ggsave("plots/strong_up3408_change.eps", ggplot_strong_up3408)
#ggsave("plots/strong_up3408_change.png", ggplot_strong_up3408, dpi = 300)

ggplot_weak_down900 <- get_transcription_grid_tricolor(weak_down900_change, 
                                                plot_title = "Weak Down900 Fold Change", 
                                                fillvar = "change",
                                                x_label = "Gyrase Activity (Lk/kb/s)",
                                                y_label = "TopoI Activity (Lk/kb/s)",
                                                fill_label = "Log2 Fold Change",
                                                gyr_scale = 4e3,
                                                topoI_scale = 1e3,
                                                fill_min = -2,
                                                fill_max = 2,
                                                opt_x = opt_gyr,
                                                opt_y = opt_top) +
  theme_classic()
#ggsave("plots/weak_down900_change.eps", ggplot_weak_down900)
#ggsave("plots/weak_down900_change.png", ggplot_weak_down900, dpi = 300)

ggplot_mod_down900 <- get_transcription_grid_tricolor(mod_down900_change, 
                                               plot_title = "Moderate Down900 Fold Change", 
                                               fillvar = "change",
                                               x_label = "Gyrase Activity (Lk/kb/s)",
                                               y_label = "TopoI Activity (Lk/kb/s)",
                                               fill_label = "Log2 Fold Change",
                                               gyr_scale = 4e3,
                                               topoI_scale = 1e3,
                                               fill_min = -2,
                                               fill_max = 2,
                                               opt_x = opt_gyr,
                                               opt_y = opt_top) +
  theme_classic()
#ggsave("plots/mod_down900_change.eps", ggplot_mod_down900)
#ggsave("plots/mod_down900_change.png", ggplot_mod_down900, dpi = 300)

ggplot_strong_down900 <- get_transcription_grid_tricolor(strong_down900_change, 
                                                  plot_title = "Strong Down900 Fold Change", 
                                                  fillvar = "change",
                                                  x_label = "Gyrase Activity (Lk/kb/s)",
                                                  y_label = "TopoI Activity (Lk/kb/s)",
                                                  fill_label = "Log2 Fold Change",
                                                  gyr_scale = 4e3,
                                                  topoI_scale = 1e3,
                                                  fill_min = -2,
                                                  fill_max = 2,
                                                  opt_x = opt_gyr,
                                                  opt_y = opt_top) +
  theme_classic()
#ggsave("plots/strong_down900_change.eps", ggplot_strong_down900)
#ggsave("plots/strong_down900_change.png", ggplot_strong_down900, dpi = 300)

ggplot_weak_down3200 <- get_transcription_grid_tricolor(weak_down3200_change, 
                                                plot_title = "Weak Down3200 Fold Change", 
                                                fillvar = "change",
                                                x_label = "Gyrase Activity (Lk/kb/s)",
                                                y_label = "TopoI Activity (Lk/kb/s)",
                                                fill_label = "Log2 Fold Change",
                                                gyr_scale = 4e3,
                                                topoI_scale = 1e3,
                                                fill_min = -2,
                                                fill_max = 2,
                                                opt_x = opt_gyr,
                                                opt_y = opt_top) +
  theme_classic()
#ggsave("plots/weak_down3200_change.eps", ggplot_weak_down3200)
#ggsave("plots/weak_down3200_change.png", ggplot_weak_down3200, dpi = 300)

ggplot_mod_down3200 <- get_transcription_grid_tricolor(mod_down3200_change, 
                                               plot_title = "Moderate Down3200 Fold Change", 
                                               fillvar = "change",
                                               x_label = "Gyrase Activity (Lk/kb/s)",
                                               y_label = "TopoI Activity (Lk/kb/s)",
                                               fill_label = "Log2 Fold Change",
                                               gyr_scale = 4e3,
                                               topoI_scale = 1e3,
                                               fill_min = -2,
                                               fill_max = 2,
                                               opt_x = opt_gyr,
                                               opt_y = opt_top) +
  theme_classic()
#ggsave("plots/mod_down3200_change.eps", ggplot_mod_down3200)
#ggsave("plots/mod_down3200_change.png", ggplot_mod_down3200, dpi = 300)

ggplot_strong_down3200 <- get_transcription_grid_tricolor(strong_down3200_change, 
                                                  plot_title = "Strong Down3200 Fold Change", 
                                                  fillvar = "change",
                                                  x_label = "Gyrase Activity (Lk/kb/s)",
                                                  y_label = "TopoI Activity (Lk/kb/s)",
                                                  fill_label = "Log2 Fold Change",
                                                  gyr_scale = 4e3,
                                                  topoI_scale = 1e3,
                                                  fill_min = -2,
                                                  fill_max = 2,
                                                  opt_x = opt_gyr,
                                                  opt_y = opt_top) +
  theme_classic()
#ggsave("plots/strong_down3200_change.eps", ggplot_strong_down3200)
#ggsave("plots/strong_down3200_change.png", ggplot_strong_down3200, dpi = 300)



# 2025-03-26
# plotting over barrier distance
setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2025-03-25_barrier-sensitivity-fine")

weak_fine_df = import_data("weak", start_time = 1800) %>%
  summarize_data() %>%
  mutate(Promoter = "Weak")

weak_base_rate = weak_fine_df$mean_rate[weak_fine_df$up == 320 & weak_fine_df$down == 250]

weak_fine_df = weak_fine_df %>% mutate(
  sensitivity = mean_rate / weak_base_rate
)

mod_fine_df = import_data("mod", start_time = 1800) %>%
  summarize_data() %>%
  mutate(Promoter = "Moderate")

mod_base_rate = mod_fine_df$mean_rate[mod_fine_df$up == 320 & mod_fine_df$down == 250]

mod_fine_df = mod_fine_df %>% mutate(
  sensitivity = mean_rate / mod_base_rate
)

strong_fine_df = import_data("strong", start_time = 1800) %>%
  summarize_data() %>%
  mutate(Promoter = "Strong")

strong_base_rate = strong_fine_df$mean_rate[strong_fine_df$up == 320 & strong_fine_df$down == 250]

strong_fine_df = strong_fine_df %>% mutate(
  sensitivity = mean_rate / strong_base_rate
)


barrier_data = rbind(weak_fine_df, mod_fine_df, strong_fine_df)

# Custom theme without legend (used in upstream plot)
custom_theme_no_legend <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
)

# Custom theme with legend in top-right corner (used in downstream plot)
custom_theme_with_legend <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.title = element_tend(size = 0.125 * 72),
  legend.text = element_text(size = 0.105 * 72),
  legend.position = c(1, 1),
  legend.justification = c("right", "top")
)

# Upstream plot
upstream_plot <- ggplot(subset(barrier_data, down == 250), aes(x = up, y = sensitivity, color = Promoter)) +
  geom_line() +
  geom_point(shape = 16) +
  xlab("Upstream Distance") +
  ylab("Sensitivity") +
  scale_y_continuous(limits = c(0.8, 1.8)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme_no_legend# +
  #scale_x_log10()

# Downstream plot
downstream_plot <- ggplot(subset(barrier_data, up == 320), aes(x = down, y = sensitivity, color = Promoter)) +
  geom_line() +
  geom_point(shape = 16) +
  xlab("Downstream Distance") +
  scale_y_continuous(limits = c(0.8, 1.8)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme_with_legend +
  theme(axis.title.y = element_blank())# +
  #scale_x_log10()

# Arrange side by side
combined_plot <- arrangeGrob(upstream_plot, downstream_plot, ncol = 2, nrow = 1)
grid::grid.draw(combined_plot)

# Save to file
ggsave("sensitivity_plots_linear.png", combined_plot, dpi = 300, width = 5.2, height = 3, units = "in")
ggsave("sensitivity_plots_linear.eps", combined_plot, dpi = 300, width = 5.2, height = 3, units = "in")

# 2025-04-14
# Repeat plots over barrier distance for kG = kG*
setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2025-04-13_barrier-sensitivity-fine-2")

weak_fine_df = import_data("weak", start_time = 1800) %>%
  summarize_data() %>%
  mutate(Promoter = "Weak")

weak_base_rate = weak_fine_df$mean_rate[weak_fine_df$up == 320 & weak_fine_df$down == 250]

weak_fine_df = weak_fine_df %>% mutate(
  sensitivity = mean_rate / weak_base_rate
)

mod_fine_df = import_data("mod", start_time = 1800) %>%
  summarize_data() %>%
  mutate(Promoter = "Moderate")

mod_base_rate = mod_fine_df$mean_rate[mod_fine_df$up == 320 & mod_fine_df$down == 250]

mod_fine_df = mod_fine_df %>% mutate(
  sensitivity = mean_rate / mod_base_rate
)

strong_fine_df = import_data("strong", start_time = 1800) %>%
  summarize_data() %>%
  mutate(Promoter = "Strong")

strong_base_rate = strong_fine_df$mean_rate[strong_fine_df$up == 320 & strong_fine_df$down == 250]

strong_fine_df = strong_fine_df %>% mutate(
  sensitivity = mean_rate / strong_base_rate
)


barrier_data = rbind(weak_fine_df, mod_fine_df, strong_fine_df)

# Custom theme without legend (used in upstream plot)
custom_theme_no_legend <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
)

# Custom theme with legend in top-right corner (used in downstream plot)
custom_theme_with_legend <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  #legend.title = element_tend(size = 0.125 * 72),
  legend.text = element_text(size = 0.105 * 72),
  legend.position = c(1, 1),
  legend.justification = c("right", "top")
)

# Upstream plot
upstream_plot <- ggplot(subset(barrier_data, down == 250), aes(x = up, y = sensitivity, color = Promoter)) +
  geom_line() +
  geom_point(shape = 16) +
  xlab("Upstream Distance") +
  ylab("Sensitivity") +
  scale_y_continuous(limits = c(0.79, 1.8)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme_no_legend

# Downstream plot
downstream_plot <- ggplot(subset(barrier_data, up == 320), aes(x = down, y = sensitivity, color = Promoter)) +
  geom_line() +
  geom_point(shape = 16) +
  xlab("Downstream Distance") +
  scale_y_continuous(limits = c(0.79, 1.8)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme_with_legend +
  theme(axis.title.y = element_blank())

# Arrange side by side
combined_plot <- arrangeGrob(upstream_plot, downstream_plot, ncol = 2, nrow = 1)
grid::grid.draw(combined_plot)

# Save to file
ggsave("sensitivity_plots_linear.png", combined_plot, dpi = 300, width = 5.2, height = 3, units = "in")
ggsave("sensitivity_plots_linear.eps", combined_plot, dpi = 300, width = 5.2, height = 3, units = "in")

