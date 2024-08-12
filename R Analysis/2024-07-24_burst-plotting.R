library(dplyr)
library(tidyr)
library(ggridges)

setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

top_dir = "lansG-0.000001190476_lansT-0.00004285714" #mean 3.57437076, fano 4.7456208
cont_dir = "lansG-0.00000595238_lansT-0.00002380952" #mean 17.87291215, fano 0.9434300
gyr_dir = "lansG-0.00001071428_lansT-0.000004761905" #mean 7.04744718, fano 1.4842676

top_lansG = 0.000001190476
top_lansT = 0.00004285714
cont_lansG = 0.00000595238
cont_lansT = 0.00002380952
gyr_lansG = 0.00001071428
gyr_lansT = 0.000004761905


#data import

strong_far_copy = import_copy_data(file.path("strong-far/data"), start_time = 1800, end_time = 45000)
strong_far_copy = strong_far_copy %>% mutate(fano_all = fano)
strong_far_copy$fano[strong_far_copy$lansG == 0 | strong_far_copy$lansT == 0] <- NaN

top_traj = read.table(file.path("strong-far/data", top_dir, "traj_promoter.txt"), header = TRUE, sep = "")
cont_traj = read.table(file.path("strong-far/data", cont_dir, "traj_promoter.txt"), header = TRUE, sep = "")
gyr_traj = read.table(file.path("strong-far/data", gyr_dir, "traj_promoter.txt"), header = TRUE, sep = "")

top_copy = read.csv(file.path("strong-far/data", top_dir, "copy_nums_.csv"), header = TRUE, sep = ",")
top_copy = top_copy %>%
  mutate(sigma_up = top_traj$sigma_up[1:nrow(top_copy)])
cont_copy = read.csv(file.path("strong-far/data", cont_dir, "copy_nums_.csv"), header = TRUE, sep = ",")
cont_copy = cont_copy %>%
  mutate(sigma_up = cont_traj$sigma_up[1:nrow(cont_copy)])
gyr_copy = read.csv(file.path("strong-far/data", gyr_dir, "copy_nums_.csv"), header = TRUE, sep = ",")
gyr_copy = gyr_copy %>%
  mutate(sigma_up = gyr_traj$sigma_up[1:nrow(gyr_copy)])


#fano grid
get_transcription_grid(strong_far_copy, 
                       plot_title = paste("Fano Factors for strong-far"), 
                       fillvar = "fano",
                       x_label = "Gyrase Activity (Lk/kb/s)",
                       y_label = "TopoI Activity (Lk/kb/s)",
                       fill_label = "Fano Factor",
                       topoI_scale = 1e3,
                       gyr_scale = 4e3, 
                       low_color = "white",
                       high_color = "orange",
                       opt_x = optimal_x/20,
                       opt_y = optimal_y/20,
                       fill_min = 1,
                       fill_max = 5) +
  theme_classic()



#poissonian comparisons
get_pois_plot(file.path("strong-far/data", top_dir), plot_title = "High topoI, low Gyrase", start_time = 1800, end_time = 45000, mean = F) +
  theme_classic()
get_pois_plot(file.path("strong-far/data", cont_dir), plot_title = "Optimal", start_time = 1800, end_time = 45000, mean = F) +
  theme_classic()
get_pois_plot(file.path("strong-far/data", gyr_dir), plot_title = "High Gyrase, low topoI", start_time = 1800, end_time = 45000, mean = F) +
  theme_classic()

#time courses
get_SC_plot(file.path("strong-far/data", top_dir), plot_title = "High topoI, low Gyrase", SC_min = -0.12, SC_max = 0.062, start_time = 1800, end_time = 7200, interval = 20) +
  theme_classic()
get_SC_plot(file.path("strong-far/data", cont_dir), plot_title = "Optimal", SC_min = -0.12, SC_max = 0.062, start_time = 1800, end_time = 7200, interval = 20) +
  theme_classic()
get_SC_plot(file.path("strong-far/data", gyr_dir), plot_title = "High Gyrase, low topoI", SC_min = -0.12, SC_max = 0.062, start_time = 1800, end_time = 7200, interval = 20) +
  theme_classic()

copy_timecourse(file.path("strong-far/data", top_dir), plot_title = "High topoI, low gyrase", start_time = 1800, end_time = 7200)
copy_timecourse(file.path("strong-far/data", cont_dir), plot_title = "Optimal", start_time = 1800, end_time = 7200)
copy_timecourse(file.path("strong-far/data", gyr_dir), plot_title = "low topoI, high gyrase", start_time = 1800, end_time = 7200)


# bar plots

strong_close_copy = import_copy_data(file.path("strong-close/data"), start_time = 1800, end_time = 45000) %>% 
  mutate(fano_all = fano)
strong_close_copy$fano[strong_close_copy$lansG == 0 | strong_close_copy$lansT == 0] <- NaN

strong_mod_copy = import_copy_data(file.path("strong-mod/data"), start_time = 1800, end_time = 45000) %>% 
  mutate(fano_all = fano)
strong_mod_copy$fano[strong_mod_copy$lansG == 0 | strong_mod_copy$lansT == 0] <- NaN


mod_far_copy = import_copy_data(file.path("mod-far/data"), start_time = 1800, end_time = 45000) %>%
  mutate(fano_all = fano)
mod_far_copy$fano[mod_far_copy$lansG == 0 | mod_far_copy$lansT == 0] <- NaN

weak_far_copy = import_copy_data(file.path("weak-far/data"), start_time = 1800, end_time = 45000) %>%
  mutate(fano_all = fano)
weak_far_copy$fano[weak_far_copy$lansG == 0 | weak_far_copy$lansT == 0] <- NaN

top_box_df = data.frame(
  promoter = factor(c("Strong", "Strong", "Strong", "Moderate", "Weak"), levels = c("Weak", "Moderate", "Strong")),
  barrier = as.factor(c(100, 10, 1, 100, 100)),
  fano = c(filter(strong_far_copy, lansG == top_lansG, lansT == top_lansT)$fano,
           filter(strong_mod_copy, lansG == top_lansG * 10, lansT == top_lansT * 10)$fano,
           filter(strong_close_copy, lansG == top_lansG * 100, lansT == top_lansT * 100)$fano,
           filter(mod_far_copy, lansG == top_lansG, lansT == top_lansT)$fano,
           filter(weak_far_copy, lansG == top_lansG, lansT == top_lansT)$fano)
)

cont_box_df = data.frame(
  promoter = factor(c("Strong", "Strong", "Strong", "Moderate", "Weak"), levels = c("Weak", "Moderate", "Strong")),
  barrier = as.factor(c(100, 10, 1, 100, 100)),
  fano = c(filter(strong_far_copy, lansG == cont_lansG, lansT == cont_lansT)$fano,
           filter(strong_mod_copy, lansG == cont_lansG * 10, lansT == cont_lansT * 10)$fano,
           filter(strong_close_copy, lansG == cont_lansG * 100, lansT == cont_lansT * 100)$fano,
           filter(mod_far_copy, lansG == cont_lansG, lansT == cont_lansT)$fano,
           filter(weak_far_copy, lansG == cont_lansG, lansT == cont_lansT)$fano)
)

gyr_box_df = data.frame(
  promoter = factor(c("Strong", "Strong", "Strong", "Moderate", "Weak"), levels = c("Weak", "Moderate", "Strong")),
  barrier = as.factor(c(100, 10, 1, 100, 100)),
  fano = c(filter(strong_far_copy, lansG == gyr_lansG, lansT == gyr_lansT)$fano,
           filter(strong_mod_copy, abs(lansG - gyr_lansG * 10) < 1e-8, abs(lansT - gyr_lansT * 10) < 1e-8)$fano,
           filter(strong_close_copy, abs(lansG - gyr_lansG * 100) < 1e-8, abs(lansT - gyr_lansT * 100) < 1e-8)$fano,
           filter(mod_far_copy, lansG == gyr_lansG, lansT == gyr_lansT)$fano,
           filter(weak_far_copy, lansG == gyr_lansG, lansT == gyr_lansT)$fano)
)

ggplot(data = filter(top_box_df, promoter == "Strong"), aes(x = barrier, y = fano)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Barrier Distance (kb)") +
  ylab("Fano Factor") +
  ggtitle("High TopoI, Low Gyrase") +
  ylim(0, 5) +
  geom_hline(yintercept = 1, color = "red")

ggplot(data = filter(top_box_df, barrier == 100), aes(x = promoter, y = fano)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Promoter") +
  ylab("Fano Factor") +
  ggtitle("High TopoI, Low Gyrase") +
  ylim(0, 5) +
  geom_hline(yintercept = 1, color = "red")

ggplot(data = filter(cont_box_df, promoter == "Strong"), aes(x = barrier, y = fano)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Barrier Distance (kb)") +
  ylab("Fano Factor") +
  ggtitle("Optimal") +
  ylim(0, 5) +
  geom_hline(yintercept = 1, color = "red")

ggplot(data = filter(cont_box_df, barrier == 100), aes(x = promoter, y = fano)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Promoter") +
  ylab("Fano Factor") +
  ggtitle("Optimal") +
  ylim(0, 5) +
  geom_hline(yintercept = 1, color = "red")


ggplot(data = filter(gyr_box_df, promoter == "Strong"), aes(x = barrier, y = fano)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Barrier Distance (kb)") +
  ylab("Fano Factor") +
  ggtitle("Low TopoI, High Gyrase") +
  ylim(0, 5) +
  geom_hline(yintercept = 1, color = "red")

ggplot(data = filter(gyr_box_df, barrier == 100), aes(x = promoter, y = fano)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Promoter") +
  ylab("Fano Factor") +
  ggtitle("Low topoI, High Gyrase") +
  ylim(0, 5) +
  geom_hline(yintercept = 1, color = "red")

