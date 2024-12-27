library(ggridges)

setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

top_dir = "lansG-0.00001190476_lansT-0.0002380952" #mean 1.427475, fano 2.00111
cont_dir = "lansG-0.0000595238_lansT-0.0002380952" #mean 3.618158, fano 1.600801
gyr_dir = "lansG-0.0000595238_lansT-0.0000952381" #mean 4.754099, fano 0.8794382

sigma_s = 0.061666666666667

opt_gyr = 1.1905e-5/2
opt_top = 4.7619e-5/2

top_lansG = 0.000001190476
top_lansT = 0.00002380952
cont_lansG = 0.00000595238
cont_lansT = 0.00002380952
gyr_lansG = 0.00000595238
gyr_lansT = 0.00000952381

#data import

copy_data = import_copy_data(file.path("mod-mod/data"))
copy_data$fano[copy_data$lansG == 0 | copy_data$lansT == 0] <- NaN

top_traj = read.table(file.path("mod-mod/data", top_dir, "traj_promoter.txt"), header = TRUE, sep = "")
cont_traj = read.table(file.path("mod-mod/data", cont_dir, "traj_promoter.txt"), header = TRUE, sep = "")
gyr_traj = read.table(file.path("mod-mod/data", gyr_dir, "traj_promoter.txt"), header = TRUE, sep = "")

top_copy = read.csv(file.path("mod-mod/data", top_dir, "copy_nums_.csv"), header = TRUE, sep = ",")
top_copy = top_copy %>%
  mutate(sigma_up = top_traj$sigma_up[1:nrow(top_copy)])
cont_copy = read.csv(file.path("mod-mod/data", cont_dir, "copy_nums_.csv"), header = TRUE, sep = ",")
cont_copy = cont_copy %>%
  mutate(sigma_up = cont_traj$sigma_up[1:nrow(cont_copy)])
gyr_copy = read.csv(file.path("mod-mod/data", gyr_dir, "copy_nums_.csv"), header = TRUE, sep = ",")
gyr_copy = gyr_copy %>%
  mutate(sigma_up = gyr_traj$sigma_up[1:nrow(gyr_copy)])



#fano grid
get_transcription_grid(copy_data, 
                       plot_title = paste("Fano Factors for Moderate Promoter"), 
                       fillvar = "fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Fano Factor",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3, 
                       low_color = "white",
                       high_color = "orange",
                       opt_x = opt_gyr*10,
                       opt_y = opt_top*10,
                       fill_min = 1,
                       fill_max = 3) +
  theme_classic()

#poissonian comparisons
get_pois_plot(file.path("mod-mod/data", top_dir), plot_title = "High topoI, low Gyrase", start_time = 1800, end_time = 45000, mean = F) +
  theme_classic()
get_pois_plot(file.path("mod-mod/data", cont_dir), plot_title = "Optimal", start_time = 1800, end_time = 45000, mean = F) +
  theme_classic()
get_pois_plot(file.path("mod-mod/data", gyr_dir), plot_title = "High Gyrase, low topoI", start_time = 1800, end_time = 45000, mean = F) +
  theme_classic()

#time courses
get_SC_plot(file.path("mod-mod/data", top_dir), plot_title = "High topoI, low Gyrase", SC_min = -0.12, SC_max = 0.062, start_time = 1800, end_time = 5400, interval = 20) +
  theme_classic() +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "green") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "green")
get_SC_plot(file.path("mod-mod/data", cont_dir), plot_title = "Optimal", SC_min = -0.12, SC_max = 0.062, start_time = 1800, end_time = 5400, interval = 20) +
  theme_classic() +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "green") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "green")
get_SC_plot(file.path("mod-mod/data", gyr_dir), plot_title = "High Gyrase, low topoI", SC_min = -0.12, SC_max = 0.062, start_time = 1800, end_time = 5400, interval = 20) +
  theme_classic() +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "green") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "green")

copy_timecourse(file.path("mod-mod/data", top_dir), plot_title = "High topoI, low gyrase", start_time = 1800, end_time = 5400, copy_max = 15)
copy_timecourse(file.path("mod-mod/data", cont_dir), plot_title = "Optimal", start_time = 1800, end_time = 5400, copy_max = 15)
copy_timecourse(file.path("mod-mod/data", gyr_dir), plot_title = "low topoI, high gyrase", start_time = 1800, end_time = 5400, copy_max = 15)
