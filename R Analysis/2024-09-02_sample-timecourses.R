setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7/strong-mod/data")
library("dplyr")

sigma_s = 0.061666666666667
dir1 = "lansG-0.0000595238_lansT-0.0002380952"
dir2 = "lansG-0.0000595238_lansT-0.0000952381"
dir3 = "lansG-0.00001190476_lansT-0.0002380952"

get_SC_plot(dir1, plot_title = "Continuous Elongation", SC_min = -0.1, SC_max = 0.062, start_time = 1800, end_time = 3600, interval = 20) + 
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "green") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "green") +
  theme_classic()
get_SC_plot(dir2, plot_title = "Insufficient TopoI", SC_min = -0.1, SC_max = 0.062, start_time = 1800, end_time = 3600, interval = 20) +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "green") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "green") +
  theme_classic()
get_SC_plot(dir3, plot_title = "Insufficient Gyrase", SC_min = -0.1, SC_max = 0.062, start_time = 1800, end_time = 3600, interval = 20) +
  geom_hline(yintercept = sigma_s, linetype = "dashed", color = "green") +
  geom_hline(yintercept = -sigma_s, linetype = "dashed", color = "green") +
  theme_classic()

#2024-12-09 copy timecourses

copy_timecourse(dir1, plot_title = "Continuous Elongation", start_time = 1800, end_time = 3600, copy_max = 30)
copy_timecourse(dir2, plot_title = "Insufficient TopoI", start_time = 1800, end_time = 3600, copy_max = 30)
copy_timecourse(dir3, plot_title = "Insufficient Gyrase", start_time = 1800, end_time = 3600, copy_max = 30)
