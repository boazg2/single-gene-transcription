library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

# Constants

opt_gyr = 1.1905e-5/2
opt_top = 4.7619e-5/2

# Data import
weak_close_copy = import_copy_data(file.path("weak-close/data"), start_time = 1800, end_time = 45000)
weak_mod_copy = import_copy_data(file.path("weak-mod/data"), start_time = 1800, end_time = 45000)
weak_far_copy = import_copy_data(file.path("weak-far/data"), start_time = 1800, end_time = 45000)

mod_close_copy = import_copy_data(file.path("mod-close/data"), start_time = 1800, end_time = 45000)
mod_mod_copy = import_copy_data(file.path("mod-mod/data"), start_time = 1800, end_time = 45000)
mod_far_copy = import_copy_data(file.path("mod-far/data"), start_time = 1800, end_time = 45000)

strong_close_copy = import_copy_data(file.path("strong-close/data"), start_time = 1800, end_time = 45000)
strong_mod_copy = import_copy_data(file.path("strong-mod/data"), start_time = 1800, end_time = 45000)
strong_far_copy = import_copy_data(file.path("strong-far/data"), start_time = 1800, end_time = 45000)

weak_close_copy = mutate(weak_close_copy, log_fano = log(fano, base = 2))
weak_mod_copy = mutate(weak_mod_copy, log_fano = log(fano, base = 2))
weak_far_copy = mutate(weak_far_copy, log_fano = log(fano, base = 2))

mod_close_copy = mutate(mod_close_copy, log_fano = log(fano, base = 2))
mod_mod_copy = mutate(mod_mod_copy, log_fano = log(fano, base = 2))
mod_far_copy = mutate(mod_far_copy, log_fano = log(fano, base = 2))

strong_close_copy = mutate(strong_close_copy, log_fano = log(fano, base = 2))
strong_mod_copy = mutate(strong_mod_copy, log_fano = log(fano, base = 2))
strong_far_copy = mutate(strong_far_copy, log_fano = log(fano, base = 2))

weak_close_copy$log_fano[weak_close_copy$lansG == 0 | weak_close_copy$lansT == 0] <- NaN
weak_mod_copy$log_fano[weak_mod_copy$lansG == 0 | weak_mod_copy$lansT == 0] <- NaN
weak_far_copy$log_fano[weak_far_copy$lansG == 0 | weak_far_copy$lansT == 0] <- NaN

mod_close_copy$log_fano[mod_close_copy$lansG == 0 | mod_close_copy$lansT == 0] <- NaN
mod_mod_copy$log_fano[mod_mod_copy$lansG == 0 | mod_mod_copy$lansT == 0] <- NaN
mod_far_copy$log_fano[mod_far_copy$lansG == 0 | mod_far_copy$lansT == 0] <- NaN

strong_close_copy$log_fano[strong_close_copy$lansG == 0 | strong_close_copy$lansT == 0] <- NaN
strong_mod_copy$log_fano[strong_mod_copy$lansG == 0 | strong_mod_copy$lansT == 0] <- NaN
strong_far_copy$log_fano[strong_far_copy$lansG == 0 | strong_far_copy$lansT == 0] <- NaN

weak_close_copy$noise[weak_close_copy$lansG == 0 | weak_close_copy$lansT == 0] <- NaN
weak_mod_copy$noise[weak_mod_copy$lansG == 0 | weak_mod_copy$lansT == 0] <- NaN
weak_far_copy$noise[weak_far_copy$lansG == 0 | weak_far_copy$lansT == 0] <- NaN

mod_close_copy$noise[mod_close_copy$lansG == 0 | mod_close_copy$lansT == 0] <- NaN
mod_mod_copy$noise[mod_mod_copy$lansG == 0 | mod_mod_copy$lansT == 0] <- NaN
mod_far_copy$noise[mod_far_copy$lansG == 0 | mod_far_copy$lansT == 0] <- NaN

strong_close_copy$noise[strong_close_copy$lansG == 0 | strong_close_copy$lansT == 0] <- NaN
strong_mod_copy$noise[strong_mod_copy$lansG == 0 | strong_mod_copy$lansT == 0] <- NaN
strong_far_copy$noise[strong_far_copy$lansG == 0 | strong_far_copy$lansT == 0] <- NaN


# Plotting

get_transcription_grid(weak_close_copy, 
                       plot_title = paste("Weak Close"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 100,
                       opt_y = opt_top * 100,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()
get_transcription_grid(weak_mod_copy, 
                       plot_title = paste("Weak Mod"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()
get_transcription_grid(weak_far_copy, 
                       plot_title = paste("Weak Far"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 1,
                       opt_y = opt_top * 1,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()

get_transcription_grid(mod_close_copy, 
                       plot_title = paste("mod Close"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 100,
                       opt_y = opt_top * 100,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()
get_transcription_grid(mod_mod_copy, 
                       plot_title = paste("mod Mod"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()
get_transcription_grid(mod_far_copy, 
                       plot_title = paste("mod Far"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 1,
                       opt_y = opt_top * 1,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()

get_transcription_grid(strong_close_copy, 
                       plot_title = paste("strong Close"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 100,
                       opt_y = opt_top * 100,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()
get_transcription_grid(strong_mod_copy, 
                       plot_title = paste("strong Mod"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()
get_transcription_grid(strong_far_copy, 
                       plot_title = paste("strong Far"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 1,
                       opt_y = opt_top * 1,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5) +
  theme_classic()

# 2024-09-12 Noise 

get_transcription_grid(weak_close_copy, 
                       plot_title = paste("Weak Close"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 100,
                       opt_y = opt_top * 100,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("weak-close_noise.png")

get_transcription_grid(weak_mod_copy, 
                       plot_title = paste("Weak Mod"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("weak-mod_noise.png")

get_transcription_grid(weak_far_copy, 
                       plot_title = paste("Weak Far"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 1,
                       opt_y = opt_top * 1,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("weak-far_noise.png")

get_transcription_grid(mod_close_copy, 
                       plot_title = paste("mod Close"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 100,
                       opt_y = opt_top * 100,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("mod-close_noise.png")

get_transcription_grid(mod_mod_copy, 
                       plot_title = paste("mod Mod"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("mod-mod_noise.png")

get_transcription_grid(mod_far_copy, 
                       plot_title = paste("mod Far"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 1,
                       opt_y = opt_top * 1,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("mod-far_noise.png")

get_transcription_grid(strong_close_copy, 
                       plot_title = paste("strong Close"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 100,
                       opt_y = opt_top * 100,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("strong-close_noise.png")

get_transcription_grid(strong_mod_copy, 
                       plot_title = paste("strong Mod"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 10,
                       opt_y = opt_top * 10,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("strong-mod_noise.png")

get_transcription_grid(strong_far_copy, 
                       plot_title = paste("strong Far"), 
                       fillvar = "noise",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Noise",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,                       
                       opt_x = opt_gyr * 1,
                       opt_y = opt_top * 1,
                       high_color = "blue",
                       low_color = "white") +
  theme_classic()
ggsave("strong-far_noise.png")



#2024-12-09: log-mean log-fano plotting

all_data <- bind_rows(
  weak_close_copy %>% mutate(Promoter = "weak", Distance = "close"),
  weak_mod_copy %>% mutate(Promoter = "weak", Distance = "mod"),
  weak_far_copy %>% mutate(Promoter = "weak", Distance = "far"),
  
  mod_close_copy %>% mutate(Promoter = "mod", Distance = "close"),
  mod_mod_copy %>% mutate(Promoter = "mod", Distance = "mod"),
  mod_far_copy %>% mutate(Promoter = "mod", Distance = "far"),
  
  strong_close_copy %>% mutate(Promoter = "strong", Distance = "close"),
  strong_mod_copy %>% mutate(Promoter = "strong", Distance = "mod"),
  strong_far_copy %>% mutate(Promoter = "strong", Distance = "far")
) %>% filter(lansT != 0, lansG != 0)

all_data$Promoter = factor(all_data$Promoter, levels = c("weak", "mod", "strong"))
all_data$Distance = factor(all_data$Distance, levels = c("close", "mod", "far"))


# Create the plot
ggplot(all_data, aes(x = mean, y = fano, shape = Promoter, color = Distance)) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Mean (log scale)",
    y = "Fano (log scale)",
    title = "Mean vs. Fano by Promoter and Barrier Distance",
    color = "Barrier Distance",
    shape = "Promoter"
  ) +
  scale_shape_manual(values = c("weak" = 15, "mod" = 17, "strong" = 19)) + # Square, Triangle, Circle
  scale_color_manual(values = c("close" = "red", "mod" = "green", "far" = "blue")) +
  scale_y_log10(breaks = c(0.5, 1, 2, 4, 8))


noise_model = lm(log10(noise) ~ log10(mean), data = all_data)
x_pois = seq(0.025, 20, 0.005)
y_pois = 1/x_pois
# repeat for noise
ggplot(all_data, aes(x = mean, y = noise, shape = Promoter, color = Distance)) +
  geom_point(size = 1.7) +
  geom_line(data = data.frame(mean = x_pois, noise = y_pois, Promoter = "mod", Distance = "mod"), color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Mean (log scale)",
    y = "Noise (log scale)",
    title = "Mean vs. Noise by Promoter and Barrier Distance",
    color = "Barrier Distance",
    shape = "Promoter"
  ) +
  scale_shape_manual(values = c("weak" = 15, "mod" = 17, "strong" = 19)) + # Square, Triangle, Circle
  scale_color_manual(values = c("close" = "red", "mod" = "green", "far" = "blue"))# +
  #scale_y_log10(breaks = c(0.1, 1, 10)) +
  #scale_x_log10(breaks = c(0.1, 1, 10))


