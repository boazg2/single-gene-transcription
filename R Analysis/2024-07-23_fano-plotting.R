library(dplyr)
library(tidyr)
library(ggplot2)

setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

# Constants

opt_gyr = 1.1905e-5/2
opt_top = 4.7619e-5/2
# Data import
weak_close_copy = import_copy_data(file.path("weak-close/data"), start_time = 1800, end_time = 45000)
weak_mod_copy = import_copy_data(file.path("weak-mod/data"), start_time = 1800, end_time = 45000)
weak_far_copy = import_copy_data(file.path("weak-far/data"), start_time = 1800, end_time = 45000)
weak_close_low_copy = import_copy_data(file.path("weak-close_low-topo/data"), start_time = 1800, end_time = 45000)
weak_mod_low_copy = import_copy_data(file.path("weak-mod_low-topo/data"), start_time = 1800, end_time = 45000)

mod_close_copy = import_copy_data(file.path("mod-close/data"), start_time = 1800, end_time = 45000)
mod_mod_copy = import_copy_data(file.path("mod-mod/data"), start_time = 1800, end_time = 45000)
mod_far_copy = import_copy_data(file.path("mod-far/data"), start_time = 1800, end_time = 45000)
mod_close_low_copy = import_copy_data(file.path("mod-close_low-topo/data"), start_time = 1800, end_time = 45000)
mod_mod_low_copy = import_copy_data(file.path("mod-mod_low-topo/data"), start_time = 1800, end_time = 45000)

strong_close_copy = import_copy_data(file.path("strong-close/data"), start_time = 1800, end_time = 45000)
strong_mod_copy = import_copy_data(file.path("strong-mod/data"), start_time = 1800, end_time = 45000)
strong_far_copy = import_copy_data(file.path("strong-far/data"), start_time = 1800, end_time = 45000)
strong_close_low_copy = import_copy_data(file.path("strong-close_low-topo/data"), start_time = 1800, end_time = 45000)
strong_mod_low_copy = import_copy_data(file.path("strong-mod_low-topo/data"), start_time = 1800, end_time = 45000)

weak_close_copy = mutate(weak_close_copy, log_fano = log(fano, base = 2))
weak_mod_copy = mutate(weak_mod_copy, log_fano = log(fano, base = 2))
weak_far_copy = mutate(weak_far_copy, log_fano = log(fano, base = 2))
weak_close_low_copy = mutate(weak_close_low_copy, log_fano = log(fano, base = 2))
weak_mod_low_copy = mutate(weak_mod_low_copy, log_fano = log(fano, base = 2))

mod_close_copy = mutate(mod_close_copy, log_fano = log(fano, base = 2))
mod_mod_copy = mutate(mod_mod_copy, log_fano = log(fano, base = 2))
mod_far_copy = mutate(mod_far_copy, log_fano = log(fano, base = 2))
mod_close_low_copy = mutate(mod_close_low_copy, log_fano = log(fano, base = 2))
mod_mod_low_copy = mutate(mod_mod_low_copy, log_fano = log(fano, base = 2))

strong_close_copy = mutate(strong_close_copy, log_fano = log(fano, base = 2))
strong_mod_copy = mutate(strong_mod_copy, log_fano = log(fano, base = 2))
strong_far_copy = mutate(strong_far_copy, log_fano = log(fano, base = 2))
strong_close_low_copy = mutate(strong_close_low_copy, log_fano = log(fano, base = 2))
strong_mod_low_copy = mutate(strong_mod_low_copy, log_fano = log(fano, base = 2))

weak_close_copy$log_fano[weak_close_copy$lansG == 0 | weak_close_copy$lansT == 0] <- NaN
weak_mod_copy$log_fano[weak_mod_copy$lansG == 0 | weak_mod_copy$lansT == 0] <- NaN
weak_far_copy$log_fano[weak_far_copy$lansG == 0 | weak_far_copy$lansT == 0] <- NaN
weak_close_low_copy$log_fano[weak_close_low_copy$lansG == 0 | weak_close_low_copy$lansT == 0] <- NaN
weak_mod_low_copy$log_fano[weak_mod_low_copy$lansG == 0 | weak_mod_low_copy$lansT == 0] <- NaN

mod_close_copy$log_fano[mod_close_copy$lansG == 0 | mod_close_copy$lansT == 0] <- NaN
mod_mod_copy$log_fano[mod_mod_copy$lansG == 0 | mod_mod_copy$lansT == 0] <- NaN
mod_far_copy$log_fano[mod_far_copy$lansG == 0 | mod_far_copy$lansT == 0] <- NaN
mod_close_low_copy$log_fano[mod_close_low_copy$lansG == 0 | mod_close_low_copy$lansT == 0] <- NaN
mod_mod_low_copy$log_fano[mod_mod_low_copy$lansG == 0 | mod_mod_low_copy$lansT == 0] <- NaN

strong_close_copy$log_fano[strong_close_copy$lansG == 0 | strong_close_copy$lansT == 0] <- NaN
strong_mod_copy$log_fano[strong_mod_copy$lansG == 0 | strong_mod_copy$lansT == 0] <- NaN
strong_far_copy$log_fano[strong_far_copy$lansG == 0 | strong_far_copy$lansT == 0] <- NaN
strong_close_low_copy$log_fano[strong_close_low_copy$lansG == 0 | strong_close_low_copy$lansT == 0] <- NaN
strong_mod_low_copy$log_fano[strong_mod_low_copy$lansG == 0 | strong_mod_low_copy$lansT == 0] <- NaN




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
get_transcription_grid(weak_close_low_copy, 
                       plot_title = paste("Weak Close Low Topo"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5,
                       opt = F) +
  theme_classic()
get_transcription_grid(weak_mod_low_copy, 
                       plot_title = paste("Weak Mod Low Topo"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5,
                       opt = F) +
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
get_transcription_grid(mod_close_low_copy, 
                       plot_title = paste("mod Close Low Topo"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5,
                       opt = F) +
  theme_classic()
get_transcription_grid(mod_mod_low_copy, 
                       plot_title = paste("mod Mod Low Topo"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5,
                       opt = F) +
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
get_transcription_grid(strong_close_low_copy, 
                       plot_title = paste("strong Close Low Topo"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5,
                       opt = F) +
  theme_classic()
get_transcription_grid(strong_mod_low_copy, 
                       plot_title = paste("strong Mod Low Topo"), 
                       fillvar = "log_fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "log2 Fano",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3,
                       high_color = "orange",
                       low_color = "white",
                       fill_min = 0,
                       fill_max = 2.5,
                       opt = F) +
  theme_classic()
