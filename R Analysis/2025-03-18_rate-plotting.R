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


# Calculate elongation speeds

weak_df = mutate(weak_df, speed = 1000 / mean_elong)
weak_df$speed[is.nan(weak_df$speed)] <- 0
mod_df = mutate(mod_df, speed = 1000 / mean_elong)
mod_df$speed[is.nan(mod_df$speed)] <- 0
strong_df = mutate(strong_df, speed = 1000 / mean_elong)
strong_df$speed[is.nan(strong_df$speed)] <- 0

# get steric data

weak_mod_steric = import_steric_data("weak-mod/data", start_time = 1800, ki = k_weak)
mod_mod_steric = import_steric_data("mod-mod/data", start_time = 1800, ki = k_mod)
strong_mod_steric = import_steric_data("strong-mod/data", start_time = 1800, ki = k_strong)


# scale to translate into kilobases and account for gyrase processivity
top_scale = 1e3
gyr_scale = 4e3
topo_slice = 0.23809520

# Combine dataframes with promoter strength labels
weak_df$Promoter <- "Weak"
mod_df$Promoter <- "Moderate"
strong_df$Promoter <- "Strong"

data_combined <- rbind(weak_df, mod_df, strong_df) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale)

weak_mod_steric$Promoter <- "Weak"
mod_mod_steric$Promoter <- "Moderate"
strong_mod_steric$Promoter <- "Strong"

data_steric_combined <- rbind(weak_mod_steric, mod_mod_steric, strong_mod_steric) %>%
  mutate(top_activity = lansT * top_scale, gyr_activity = lansG * gyr_scale,
         free_prom = ifelse(free_prom > 1, 1, free_prom))

# Create plots with adjusted font sizes
custom_theme <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),  # Convert inches to points, set color to black
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"  # Remove legend
)

lansT_rate <- ggplot(subset(data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  ylab("Normed Rate") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.x = element_blank())

lansG_rate <- ggplot(subset(data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())

lansT_speed <- ggplot(subset(data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  ylab("Speed (bp/s)") +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.x = element_blank())

lansG_speed <- ggplot(subset(data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())

lansT_steric <- ggplot(subset(data_steric_combined, gyr_activity == topo_slice), aes(x = top_activity, y = free_prom, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Free Promoter Fraction") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red"), breaks = c("Weak", "Moderate", "Strong")) +
  theme_classic() +
  custom_theme

lansG_steric <- ggplot(subset(data_steric_combined, top_activity == topo_slice), aes(x = gyr_activity, y = free_prom, color = Promoter)) +
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

ggsave("figures/rate_plots.png", combined_plot, dpi = 300, width = 5.2, height = 7.5, units = "in")
ggsave("figures/rate_plots.eps", combined_plot, dpi = 300, width = 5.2, height = 7.5, units = "in")

# wide plot
# Arrange in a 2x3 grid (transpose of original 3x2)
# Custom theme with appropriate font sizes
custom_theme <- theme(
  axis.text = element_text(size = 0.105 * 72, color = "black"),
  axis.title = element_text(size = 0.125 * 72),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "none"
)

# lansT_rate: TopoI Activity vs. Normed Rate
lansT_rate <- ggplot(subset(data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Normed Rate") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansG_rate: Gyrase Activity vs. Normed Rate
lansG_rate <- ggplot(subset(data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = normed_rate, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("Gyrase Activity (Lk/kb/s)") +
  ylab("Normed Rate") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansT_speed: TopoI Activity vs. Speed
lansT_speed <- ggplot(subset(data_combined, gyr_activity == topo_slice), aes(x = top_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Speed (bp/s)") +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansG_speed: Gyrase Activity vs. Speed
lansG_speed <- ggplot(subset(data_combined, top_activity == topo_slice), aes(x = gyr_activity, y = speed, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("Gyrase Activity (Lk/kb/s)") +
  ylab("Speed (bp/s)") +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansT_steric: TopoI Activity vs. Free Promoter
lansT_steric <- ggplot(subset(data_steric_combined, gyr_activity == topo_slice), aes(x = top_activity, y = free_prom, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("TopoI Activity (Lk/kb/s)") +
  ylab("Free Promoter Fraction") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme

# lansG_steric: Gyrase Activity vs. Free Promoter
lansG_steric <- ggplot(subset(data_steric_combined, top_activity == topo_slice), aes(x = gyr_activity, y = free_prom, color = Promoter)) +
  geom_line() +
  geom_point() +
  xlab("Gyrase Activity (Lk/kb/s)") +
  ylab("Free Promoter Fraction") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Weak" = "blue", "Moderate" = "green", "Strong" = "red")) +
  theme_classic() +
  custom_theme +
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"))

combined_plot_transposed <- arrangeGrob(
  lansT_rate, lansT_speed, lansT_steric,  # First row (formerly first column)
  lansG_rate, lansG_speed, lansG_steric,  # Second row (formerly second column)
  ncol = 3, nrow = 2
)

ggsave("figures/rate_plots_transposed.png", combined_plot_transposed, dpi = 300, width = 7.5, height = 5.2, units = "in")
ggsave("figures/rate_plots_transposed.eps", combined_plot_transposed, dpi = 300, width = 7.5, height = 5.2, units = "in")
