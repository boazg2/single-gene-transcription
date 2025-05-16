#all for 10 kb
library(zoo)
library(gridExtra)
setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

k_weak = 0.008
k_mod = 0.05
k_strong = 0.2

opt_gyr = 5.625e-05
opt_top = 0.000225

# data import
weak_df = import_data("weak-mod/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_df = import_data("mod-mod/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_mod)
strong_df = import_data("strong-mod/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_strong)

weak_novo_df = import_data("weak-mod_low-gyrase/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_novo_df = import_data("mod-mod_low-gyrase/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
strong_novo_df = import_data("strong-mod_low-gyrase/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)

weak_seco_df = import_data("weak-mod_low-topoI/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
mod_seco_df = import_data("mod-mod_low-topoI/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)
strong_seco_df = import_data("strong-mod_low-topoI/data", start_time = 1800) %>%
  summarize_data() %>%
  mutate(normed_rate = mean_rate / k_weak)

# calculate comparisons
weak_novo_change <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_novo_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_novo_change <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_novo_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_novo_change <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_novo_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

weak_seco_change <- bind_cols(
  weak_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  weak_seco_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

mod_seco_change <- bind_cols(
  mod_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  mod_seco_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))

strong_seco_change <- bind_cols(
  strong_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(lansG, lansT, start_rate = mean_rate),
  strong_seco_df %>% ungroup() %>% arrange(lansG, lansT) %>% select(end_rate = mean_rate)
) %>%
  mutate(reduction = log2(end_rate / start_rate))


# plot comparisons

fill_min = -4
fill_max = 4
weak_novo_plot <- get_transcription_grid_tricolor(weak_novo_change, 
                                                  plot_title = "Weak 1/5x Gyrase", 
                                                  fillvar = "reduction",
                                                  x_label = "Gyrase Activity (Lk/kb/s)",
                                                  y_label = "TopoI Activity (Lk/kb/s)",
                                                  fill_label = "Log2 Fold Change",
                                                  gyr_scale = 4e3,
                                                  topoI_scale = 1e3,
                                                  fill_min = fill_min,
                                                  fill_max = fill_max,
                                                  opt_x = opt_gyr,
                                                  opt_y = opt_top) +
  theme_classic()

mod_novo_plot <- get_transcription_grid_tricolor(mod_novo_change, 
                                                  plot_title = "Moderate 1/5x Gyrase", 
                                                  fillvar = "reduction",
                                                  x_label = "Gyrase Activity (Lk/kb/s)",
                                                  y_label = "TopoI Activity (Lk/kb/s)",
                                                  fill_label = "Log2 Fold Change",
                                                  gyr_scale = 4e3,
                                                  topoI_scale = 1e3,
                                                  fill_min = fill_min,
                                                  fill_max = fill_max,
                                                  opt_x = opt_gyr,
                                                  opt_y = opt_top) +
  theme_classic()

strong_novo_plot <- get_transcription_grid_tricolor(strong_novo_change, 
                                                 plot_title = "Strong 1/5x Gyrase", 
                                                 fillvar = "reduction",
                                                 x_label = "Gyrase Activity (Lk/kb/s)",
                                                 y_label = "TopoI Activity (Lk/kb/s)",
                                                 fill_label = "Log2 Fold Change",
                                                 gyr_scale = 4e3,
                                                 topoI_scale = 1e3,
                                                 fill_min = fill_min,
                                                 fill_max = fill_max,
                                                 opt_x = opt_gyr,
                                                 opt_y = opt_top) +
  theme_classic()

weak_seco_plot <- get_transcription_grid_tricolor(weak_seco_change, 
                                                  plot_title = "Weak 1/5x TopoI", 
                                                  fillvar = "reduction",
                                                  x_label = "Gyrase Activity (Lk/kb/s)",
                                                  y_label = "TopoI Activity (Lk/kb/s)",
                                                  fill_label = "Log2 Fold Change",
                                                  gyr_scale = 4e3,
                                                  topoI_scale = 1e3,
                                                  fill_min = fill_min,
                                                  fill_max = fill_max,
                                                  opt_x = opt_gyr,
                                                  opt_y = opt_top) +
  theme_classic()

mod_seco_plot <- get_transcription_grid_tricolor(mod_seco_change, 
                                                 plot_title = "Moderate 1/5x TopoI", 
                                                 fillvar = "reduction",
                                                 x_label = "Gyrase Activity (Lk/kb/s)",
                                                 y_label = "TopoI Activity (Lk/kb/s)",
                                                 fill_label = "Log2 Fold Change",
                                                 gyr_scale = 4e3,
                                                 topoI_scale = 1e3,
                                                 fill_min = fill_min,
                                                 fill_max = fill_max,
                                                 opt_x = opt_gyr,
                                                 opt_y = opt_top) +
  theme_classic()

strong_seco_plot <- get_transcription_grid_tricolor(strong_seco_change, 
                                                    plot_title = "Strong 1/5x TopoI", 
                                                    fillvar = "reduction",
                                                    x_label = "Gyrase Activity (Lk/kb/s)",
                                                    y_label = "TopoI Activity (Lk/kb/s)",
                                                    fill_label = "Log2 Fold Change",
                                                    gyr_scale = 4e3,
                                                    topoI_scale = 1e3,
                                                    fill_min = fill_min,
                                                    fill_max = fill_max,
                                                    opt_x = opt_gyr,
                                                    opt_y = opt_top) +
  theme_classic()


combined_antibiotic_plot <- arrangeGrob(grobs = list(
  weak_novo_plot,
  weak_seco_plot,
  mod_novo_plot,
  mod_seco_plot,
  strong_novo_plot,
  strong_seco_plot
), ncol = 2, nrow = 3)

ggsave("figures/barrier_plots.png", combined_antibiotic_plot, dpi = 300, width = 7.5, height = 5.625, units = "in")
ggsave("figures/barrier_plots.eps", combined_antibiotic_plot, dpi = 300, width = 7.5, height = 5.625, units = "in")



# plot box plots of individual simulations
# Define the filter values
lansG_plot <- 1.190476e-05*5
lansT_plot <- 2.380952e-04

# Function to filter and label data
filter_and_label <- function(df, strength, treatment) {
  df %>% 
    filter(lansG == lansG_plot, lansT == lansT_plot) %>% 
    mutate(strength = strength, treatment = treatment)
}

# Combine all data into one data frame
combined_data <- bind_rows(
  filter_and_label(weak_novo_change, "Weak", "Novobiocin"),
  filter_and_label(mod_novo_change, "Moderate", "Novobiocin"),
  filter_and_label(strong_novo_change, "Strong", "Novobiocin"),
  filter_and_label(weak_seco_change, "Weak", "Seconeolitsine"),
  filter_and_label(mod_seco_change, "Moderate", "Seconeolitsine"),
  filter_and_label(strong_seco_change, "Strong", "Seconeolitsine")
)

# Combine all data into one data frame
combined_data <- bind_rows(
  filter_and_label(weak_novo_change, "Weak", "Novobiocin"),
  filter_and_label(mod_novo_change, "Moderate", "Novobiocin"),
  filter_and_label(strong_novo_change, "Strong", "Novobiocin"),
  filter_and_label(weak_seco_change, "Weak", "Seconeolitsine"),
  filter_and_label(mod_seco_change, "Moderate", "Seconeolitsine"),
  filter_and_label(strong_seco_change, "Strong", "Seconeolitsine")
)

# Relabel strength levels only after filtering
novo_data <- combined_data %>% 
  filter(treatment == "Novobiocin") %>% 
  mutate(strength = factor(strength, levels = c("Weak", "Moderate", "Strong"), labels = c("W", "M", "S")))

seco_data <- combined_data %>% 
  filter(treatment == "Seconeolitsine") %>% 
  mutate(strength = factor(strength, levels = c("Weak", "Moderate", "Strong"), labels = c("W", "M", "S")))

# Generate individual plots
novo_bar <- ggplot(novo_data, aes(x = strength, y = reduction)) +
  geom_col(fill = "blue") +
  labs(title = "Novobiocin", x = "Promoter Strength", y = "Reduction") +
  scale_y_continuous(limits = c(-3, 4)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = rel(0.125 / 0.3528)),
    axis.text = element_text(size = rel(0.105 / 0.3528)),
    legend.position = "none"
  )

seco_bar <- ggplot(seco_data, aes(x = strength, y = reduction)) +
  geom_col(fill = "blue") +
  labs(title = "Seconeolitsine", x = "Promoter Strength", y = "Reduction") +
  scale_y_continuous(limits = c(-3, 4)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = rel(0.125 / 0.3528)),
    axis.text = element_text(size = rel(0.105 / 0.3528)),
    legend.position = "none"
  )

# Arrange the two plots side by side
combined_bars <- grid.arrange(novo_bar, seco_bar, ncol = 2)
combined_bars

ggsave("figures/antibiotic_bars_steady-state.png", combined_bars, dpi = 300, width = 5.2, height = 3, units = "in")
ggsave("figures/antibiotic_bars_steady-state-gyr.eps", combined_bars, dpi = 300, width = 5.2, height = 3, units = "in")





# 2025-04-12: difference between strong and weak genes

novo_difference <- bind_cols(
  weak_novo_change %>% select(lansG, lansT, weak_reduction = reduction),
  strong_novo_change %>% select(strong_reduction = reduction)
) %>%
  mutate(difference = strong_reduction - weak_reduction)

seco_difference <- bind_cols(
  weak_seco_change %>% select(lansG, lansT, weak_reduction = reduction),
  strong_seco_change %>% select(strong_reduction = reduction)
) %>%
  mutate(difference = strong_reduction - weak_reduction)


novo_difference_plot <- get_transcription_grid_tricolor(novo_difference, 
                                                         plot_title = "Difference 1/5x Gyrase", 
                                                         fillvar = "difference",
                                                         x_label = "Gyrase Activity (Lk/kb/s)",
                                                         y_label = "TopoI Activity (Lk/kb/s)",
                                                         fill_label = "Difference Log2 Fold Change",
                                                         gyr_scale = 4e3,
                                                         topoI_scale = 1e3,
                                                         fill_min = fill_min,
                                                         fill_max = fill_max,
                                                         opt_x = opt_gyr,
                                                         opt_y = opt_top) +
  theme_classic()

ggsave("figures/novo_difference.png", novo_difference_plot, dpi = 300, width = 5.2, height = 3, units = "in")
ggsave("figures/novo_difference.eps", novo_difference_plot, dpi = 300, width = 5.2, height = 3, units = "in")

seco_difference_plot <- get_transcription_grid_tricolor(seco_difference, 
                                                        plot_title = "Difference 1/5x TopoI", 
                                                        fillvar = "difference",
                                                        x_label = "Gyrase Activity (Lk/kb/s)",
                                                        y_label = "TopoI Activity (Lk/kb/s)",
                                                        fill_label = "Difference Log2 Fold Change",
                                                        gyr_scale = 4e3,
                                                        topoI_scale = 1e3,
                                                        fill_min = fill_min,
                                                        fill_max = fill_max,
                                                        opt_x = opt_gyr,
                                                        opt_y = opt_top) +
  theme_classic()


ggsave("figures/seco_difference.png", seco_difference_plot, dpi = 300, width = 5.2, height = 3, units = "in")
ggsave("figures/seco_difference.eps", seco_difference_plot, dpi = 300, width = 5.2, height = 3, units = "in")