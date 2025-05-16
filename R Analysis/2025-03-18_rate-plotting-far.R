#all for 10 kb
library(zoo)
library(gridExtra)
setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

k_weak = 0.008
k_mod = 0.05
k_strong = 0.2

opt_gyr = 5.625e-06
opt_top = 0.0000225

weak_far_df = import_data("weak-far/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_far_df = import_data("mod-far/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_far_df = import_data("strong-far/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)


# Calculate elongation speeds

weak_far_df = mutate(weak_far_df, speed = 1000 / mean_elong)
weak_far_df$speed[is.nan(weak_far_df$speed)] <- 0
mod_far_df = mutate(mod_far_df, speed = 1000 / mean_elong)
mod_far_df$speed[is.nan(mod_far_df$speed)] <- 0
strong_far_df = mutate(strong_far_df, speed = 1000 / mean_elong)
strong_far_df$speed[is.nan(strong_far_df$speed)] <- 0

# get steric data

weak_far_steric = import_steric_data("weak-far/data", start_time = 1800, ki = k_weak)
mod_far_steric = import_steric_data("mod-far/data", start_time = 1800, ki = k_mod)
strong_far_steric = import_steric_data("strong-far/data", start_time = 1800, ki = k_strong)


# scale to translate into kilobases and account for gyrase processivity
top_scale = 1e3
gyr_scale = 4e3
topo_slice = 0.023809520

# Combine dataframes with promoter strength labels
weak_far_df$Promoter <- "Weak"
mod_far_df$Promoter <- "Moderate"
strong_far_df$Promoter <- "Strong"

far_data_combined <- rbind(weak_far_df, mod_far_df, strong_far_df) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)

weak_far_steric$Promoter <- "Weak"
mod_far_steric$Promoter <- "Moderate"
strong_far_steric$Promoter <- "Strong"

far_data_steric_combined <- rbind(weak_far_steric, mod_far_steric, strong_far_steric) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale,
         free_prom = ifelse(free_prom > 1, 1, free_prom))

# Create plots with adjusted font sizes
custom_theme <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),  # Convert inches to points, set color to black
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"  # Remove legend
)

lansT_rate <- ggplot(subset(far_data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  ylab("Normed Rate") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.x = element_blank())

lansG_rate <- ggplot(subset(far_data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())

lansT_speed <- ggplot(subset(far_data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  ylab("Speed (bp/s)") +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.x = element_blank())

lansG_speed <- ggplot(subset(far_data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())

lansT_steric <- ggplot(subset(far_data_steric_combined, gyr_activity == topo_slice), aes(x = top_activity, y = free_prom, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Free Promoter Fraction") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme

lansG_steric <- ggplot(subset(far_data_steric_combined, top_activity == topo_slice), aes(x = gyr_activity, y = free_prom, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("Gyrase Activity (Lk/kb/s)") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"))

# Arrange in a 3x2 grid
combined_plot <- arrangeGrob(lansT_rate, lansG_rate, lansT_speed, lansG_speed, lansT_steric, lansG_steric, ncol = 2, nrow = 3)

ggsave("figures/rate_plots_far.png", combined_plot, dpi = 300, width = 5.2, height = 7.5, units = "in")
ggsave("figures/rate_plots_far.eps", combined_plot, dpi = 300, width = 5.2, height = 7.5, units = "in")


# Custom theme with consistent formatting
custom_theme <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
)

# Common x-axis tick values
x_breaks <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05)

# lansT_rate
lansT_rate <- ggplot(subset(far_data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Normed Rate") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansG_rate
lansG_rate <- ggplot(subset(far_data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("Gyrase Activity (Lk/kb/s)") +
  ylab("Normed Rate") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansT_speed
lansT_speed <- ggplot(subset(far_data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Speed (bp/s)") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansG_speed
lansG_speed <- ggplot(subset(far_data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("Gyrase Activity (Lk/kb/s)") +
  ylab("Speed (bp/s)") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansT_steric
lansT_steric <- ggplot(subset(far_data_steric_combined, gyr_activity == topo_slice), aes(x = top_activity, y = free_prom, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Free Promoter Fraction") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansG_steric
lansG_steric <- ggplot(subset(far_data_steric_combined, top_activity == topo_slice), aes(x = gyr_activity, y = free_prom, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("Gyrase Activity (Lk/kb/s)") +
  ylab("Free Promoter Fraction") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme +
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"))

# Transposed layout: 2 rows × 3 columns
combined_plot_transposed <- arrangeGrob(
  lansT_rate, lansT_speed, lansT_steric,
  lansG_rate, lansG_speed, lansG_steric,
  ncol = 3, nrow = 2
)

ggsave("figures/far_rate_plots_transposed.png", combined_plot_transposed, dpi = 300, width = 7.5, height = 5.2, units = "in")
ggsave("figures/far_rate_plots_transposed.eps", combined_plot_transposed, dpi = 300, width = 7.5, height = 5.2, units = "in")


# transposed grids

get_transcription_grid(weak_far_df, 
                       plot_title = "Weak Far Rate", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(mod_far_df, 
                       plot_title = "Mod Far Rate", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(strong_far_df, 
                       plot_title = "Strong Far Rate", 
                       fillvar = "normed_rate",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(weak_far_df, 
                       plot_title = "Weak Far Speed", 
                       fillvar = "speed",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 25,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(mod_far_df, 
                       plot_title = "Mod Far Speed", 
                       fillvar = "speed",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 25,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(strong_far_df, 
                       plot_title = "Strong Far Speed", 
                       fillvar = "speed",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 25,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(weak_far_steric, 
                       plot_title = "Weak Far free_prom", 
                       fillvar = "free_prom",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(mod_far_steric, 
                       plot_title = "Mod Far free_prom", 
                       fillvar = "free_prom",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top)

get_transcription_grid(strong_far_steric, 
                       plot_title = "Strong Far free_prom", 
                       fillvar = "free_prom",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Difference Log2 Fold Change",
                       gyr_scale = 4e3,
                       topoI_scale = 1e3,
                       fill_min = 0,
                       fill_max = 1,
                       opt_x = opt_gyr,
                       opt_y = opt_top)
