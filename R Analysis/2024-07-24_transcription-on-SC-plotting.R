library(ggridges)
library(ggplot2)
library(dplyr)

setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

t0 = 1800
sigma_s = 0.061666666666667

barrier_len = 10000
gene_len = 1000

dir = "lansG-0.0000595238_lansT-0.0002380952"


#data imports
notrans_traj = read.table(file.path("none-mod/data", dir, "traj_promoter.txt"), header = TRUE, sep = "")
weak_traj = read.table(file.path("weak-mod/data", dir, "traj_promoter.txt"), header = TRUE, sep = "")
mod_traj = read.table(file.path("mod-mod/data", dir, "traj_promoter.txt"), header = TRUE, sep = "")
strong_traj = read.table(file.path("strong-mod/data", dir, "traj_promoter.txt"), header = TRUE, sep = "")

combined_data <- data.frame(
  promoter = factor(c(rep("None", length(notrans_traj$sigma_avg)),
            rep("Weak", length(weak_traj$sigma_avg)),
            rep("Moderate", length(mod_traj$sigma_avg)),
            rep("Strong", length(strong_traj$sigma_avg))),
            levels = c("None", "Weak", "Moderate", "Strong")),
  sigma_avg = c(notrans_traj$sigma_avg, weak_traj$sigma_avg, mod_traj$sigma_avg, strong_traj$sigma_avg),
  sigma_up = c(notrans_traj$sigma_up, weak_traj$sigma_up, mod_traj$sigma_up, strong_traj$sigma_up),
  sigma_gene = c(notrans_traj$sigma_gene, weak_traj$sigma_gene, mod_traj$sigma_gene, strong_traj$sigma_gene),
  sigma_down = c(notrans_traj$sigma_down, weak_traj$sigma_down, mod_traj$sigma_down, strong_traj$sigma_down)

)

combined_summary = combined_data %>%
  group_by(promoter) %>%
  summarise(avg_mean = mean(sigma_avg),
            avg_sd = sd(sigma_avg),
            up_mean = mean(sigma_up),
            up_sd = sd(sigma_up),
            gene_mean = mean(sigma_gene),
            gene_sd = sd(sigma_gene),
            down_mean = mean(sigma_down),
            down_sd = sd(sigma_down))


#ridges
ggplot(combined_data, aes(x = sigma_avg, y = promoter, fill = promoter)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.0005) +
  theme_ridges() +
  labs(
    x = "σ",
    y = "Promoter",
    title = "Domain Supercoiling Density"
  ) +
  scale_fill_manual(values = c("None" = "blue", "Weak" = "green", "Moderate" = "orange", "Strong" = "red"))+
  xlim(c(-.1, .07)) +
  guides(fill = FALSE, color = FALSE) +
  geom_segment(data = combined_summary, aes(x = avg_mean, xend = avg_mean, 
                                            y = as.numeric(promoter), 
                                            yend = as.numeric(promoter) + 0.75, 
                                            color = promoter),
               size = 1.5) +
  scale_color_manual(values = c("None" = "black", "Weak" = "black", "Moderate" = "black", "Strong" = "black")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12)
  )


ggplot(combined_data, aes(x = sigma_up, y = promoter, fill = promoter)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.0005) +
  theme_ridges() +
  labs(
    x = "σ",
    y = "Promoter",
    title = "Upstream Supercoiling Density"
  ) +
  scale_fill_manual(values = c("None" = "blue", "Weak" = "green", "Moderate" = "orange", "Strong" = "red"))+
  xlim(c(-.1, .07)) +
  guides(fill = FALSE, color = FALSE) +
  geom_segment(data = combined_summary, aes(x = up_mean, xend = up_mean, 
                                 y = as.numeric(promoter), 
                                 yend = as.numeric(promoter) + 0.75, 
                                 color = promoter),
                                 size = 1.5) +
  scale_color_manual(values = c("None" = "black", "Weak" = "black", "Moderate" = "black", "Strong" = "black")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12)
  )


ggplot(combined_data, aes(x = sigma_gene, y = promoter, fill = promoter)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.0005) +
  theme_ridges() +
  labs(
    x = "σ",
    y = "Promoter",
    title = "Intragene Supercoiling Density"
  ) +
  scale_fill_manual(values = c("None" = "blue", "Weak" = "green", "Moderate" = "orange", "Strong" = "red"))+
  xlim(c(-.1, .07)) +
  guides(fill = FALSE, color = FALSE) +
  geom_segment(data = combined_summary, aes(x = gene_mean, xend = gene_mean, 
                                            y = as.numeric(promoter), 
                                            yend = as.numeric(promoter) + 0.75, 
                                            color = promoter),
               size = 1.5) +
  scale_color_manual(values = c("None" = "black", "Weak" = "black", "Moderate" = "black", "Strong" = "black")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12)
  )

ggplot(combined_data, aes(x = sigma_down, y = promoter, fill = promoter)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, bandwidth = 0.0005) +
  theme_ridges() +
  labs(
    x = "σ",
    y = "Promoter",
    title = "Downstream Supercoiling Density"
  ) +
  scale_fill_manual(values = c("None" = "blue", "Weak" = "green", "Moderate" = "orange", "Strong" = "red"))+
  xlim(c(-.1, .07)) +
  guides(fill = FALSE, color = FALSE) +
  geom_segment(data = combined_summary, aes(x = down_mean, xend = down_mean, 
                                            y = as.numeric(promoter), 
                                            yend = as.numeric(promoter) + 0.75, 
                                            color = promoter),
               size = 1.5) +
  scale_color_manual(values = c("None" = "black", "Weak" = "black", "Moderate" = "black", "Strong" = "black")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12)
  )

#supercoiling level grid
notrans_prom_data = import_prom_data("none-mod/data", start_time = t0) 
notrans_prom_data = notrans_prom_data %>%
  filter(lansG > 0, lansT > 0)

get_transcription_grid(notrans_prom_data, 
                       plot_title = paste("Average Supercoiling Density Without Transcription"), 
                       fillvar = "up_mean",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "σ",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,
                       low_color = "red",
                       high_color = "blue",
                       opt_x = optimal_x/2,
                       opt_y = optimal_y/2,
                       fill_min = -0.06,
                       fill_max = -0.019999) +
  theme_classic()


#SC plots
get_SC_plot(file.path("none-mod/data", dir), plot_title = "No Promoter", SC_min = -0.12, SC_max = 0.07, start_time = 1800, end_time = 2100) +
  theme_classic() +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "black")
get_SC_plot(file.path("weak-mod/data", dir), plot_title = "Weak Promoter", SC_min = -0.12, SC_max = 0.07, start_time = 1800, end_time = 2100) +
  theme_classic() +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "black")
get_SC_plot(file.path("mod-mod/data", dir), plot_title = "Moderate Promoter", SC_min = -0.12, SC_max = 0.07, start_time = 1800, end_time = 2100) +
  theme_classic() +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "black")
get_SC_plot(file.path("strong-mod/data", dir), plot_title = "Strong Promoter", SC_min = -0.12, SC_max = 0.07, start_time = 1800, end_time = 2100) +
  theme_classic() +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "black")
