library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(rle)

setwd("D:/!Xiao Lab/SC Modeling/single-gene-transcription/2024-07-03_processive-grids-7")

fix_mean = function(x_all, x_bad, n, i) {
  if(n > i){
    x_good = (n * x_all - i * x_bad) / (n - i)
  }
  else {
    x_good = x_bad
  }
  x_good
}

get_SC_plot_gene = function(dir, start_time = 0, end_time = Inf, SC_min = -0.08, SC_max = 0.08, plot_title="") {
  path = file.path(dir, 'traj_promoter.txt')
  
  SC_plot_data = read.table(path, header = TRUE, sep = "") %>%
    filter (time <= end_time, time >= start_time) %>%
    gather("location",
           "sigma",
           sigma_up,
           sigma_gene,
           sigma_down#,
           #sigma_avg
           ) %>%
    mutate(time = time)

  SC_plot = ggplot(data = SC_plot_data, aes(x = time, y = sigma, color = location, group = location)) +
    geom_line(size = 0.5) +
    ylab("Supercoiling level (σ)") + 
    xlab("Time (s)") +
    ggtitle(plot_title) +
    ylim(SC_min, SC_max) +
    xlim(start_time, min(end_time, max(SC_plot_data$time))) +
    scale_color_discrete("name" = "Location", labels = c("Downstream", "Gene", "Upstream")) +
    theme_bw() +
    theme(legend.title = element_text(size = 24),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 24),
          axis.title = element_text(size = 32),
          axis.text = element_text(color = "black", size = 24),
          legend.position = c(0.845, 0.8))

  SC_plot
}
copy_timecourse = function(dir, start_time = 0, end_time = Inf, copy_max = NaN, plot_title="", dt = 0.04) {
  path = file.path(dir, 'copy_nums_.csv')
  
  copy_data = read.csv(path, header = TRUE, sep = ",") %>%
    filter(time >= start_time, time <= end_time) %>%
    mutate(RNAP = 0)
  
  mean_path = file.path(dir, 'mean_properties.txt')
  
  mean_data = read.table(mean_path, header = TRUE, sep = "") %>%
    mutate(elong_time = mean_elong_time * transcripts_nb + lag(mean_elong_time, default = 0) * (1 - transcripts_nb)) %>%
    mutate(bind_time = time - elong_time,
           bind_index = pmax(pmin(ceiling((bind_time - start_time) / dt), nrow(copy_data)), 0),
           term_index = pmax(pmin((floor((time - start_time) / dt)), nrow(copy_data)), 0))
  
  for (i in 1:nrow(mean_data)) {
    if(mean_data$bind_index[i] < nrow(copy_data)){
      copy_data$RNAP[mean_data$bind_index[i]:mean_data$term_index[i]] = copy_data$RNAP[mean_data$bind_index[i]:mean_data$term_index[i]] + 1
    }
  }
  
  
  if(is.na(copy_max)) {
    copy_max = max(copy_data$copy_num, copy_data$RNAP)
  }
  
  copy_data = copy_data %>% pivot_longer(cols = c(copy_num, RNAP), names_to = "variable", values_to = "value")
  
  ggplot(copy_data, aes(x = time, y = value, color = variable)) +
    geom_line(size = 0.5) +
    ylab("Copy Number") + 
    xlab("Time (s)") +
    ggtitle(plot_title) +
    ylim(0, copy_max) +
    xlim(start_time, min(end_time, max(copy_data$time))) +
    scale_color_manual(values = c("green", "purple")) +
    theme_classic()
}

get_SC_plot = function(dir, start_time = 0, end_time = Inf, SC_min = -0.08, SC_max = 0.08, plot_title="", interval = 1) {
  path = file.path(dir, 'traj_promoter.txt')
  
  SC_plot_data = read.table(path, header = TRUE, sep = "") %>%
    filter (time <= end_time, time >= start_time) %>%
    gather("location",
           "sigma",
           sigma_up,
           #sigma_gene,
           sigma_down#,
           #sigma_avg
    ) %>%
    mutate(time = time) %>%
    slice(which(row_number() %% interval == 0))
  
  #SC_plot_data = SC_plot_data[seq(1, nrow(SC_plot_data), interval)]
  
  SC_plot = ggplot(data = SC_plot_data, aes(x = time, y = sigma, color = location, group = location)) +
    geom_line(size = 0.5) +
    ylab("σ") + 
    xlab("Time (s)") +
    ggtitle(plot_title) +
    ylim(SC_min, SC_max) +
    xlim(start_time, min(end_time, max(SC_plot_data$time))) +
    scale_color_discrete("name" = "Location", labels = c("Downstream", "Upstream")) +
    theme_bw() +
    theme(legend.title = element_text(size = 24),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 24),
          axis.title = element_text(size = 32),
          axis.text = element_text(color = "black", size = 24),
          legend.position = c(0.845, 0.8))
  
  SC_plot
}

get_copy_nums <- function(prod_times, end_time, dt = 0.04, RNA_hl = 120) {
  copy_df = data.frame(time = seq(0, end_time, dt), copy_num = 0)
  
  for (synth_time in prod_times) {
    RNA_lifetime = rexp(1, 0.008333333)
    decay_time = synth_time + RNA_lifetime
    
    copy_df = mutate(copy_df, 
                     copy_num = copy_num + (time >= synth_time)*(time <= decay_time))
  }
  
  return(copy_df)
}

save_copy_num <- function(datadir) {
  subdirs = list.files(path = datadir)
  # Iterate over subdirs
  for (subdir in subdirs) {
    print(paste0("starting ", subdir))
    # Construct the file path for mean_properties.txt
    file_path <- file.path(datadir, subdir, "mean_properties.txt")
    
    # Import dataframe from the CSV file
    df <- read.table(file_path, header = TRUE, sep = "")
    
    # Extract the "time" column
    time_column <- df$time
    end_time = max(time_column)
    if(end_time > 0) {
      
      # Get copy numbers using the specified function
      copy_nums_df <- get_copy_nums(time_column, end_time)
    }
    else {
      copy_nums_df <- get_copy_nums(time_column, 0)
    }
    
    # Save the resulting dataframe as a new CSV file in the current path
    new_file_path <- file.path(datadir, subdir, paste0("copy_nums_", ".csv"))
    write.csv(copy_nums_df, new_file_path, row.names = FALSE)
    
    cat("Copy numbers saved to:", new_file_path, "\n")
  }
}

get_pois_plot_df <- function(plot_data, plot_title = "", max_num = -1, mean = T) {
  
  # Extract copy_num column from the data frame
  copy_num <- plot_data$copy_num
  
  if (length(copy_num) > 0) {
    # Calculate the mean of the copy_num column
    mean_copy_num <- mean(copy_num)
    
    # Create a data frame for the histogram
    hist_data <- data.frame(breaks = seq(0.5, max(copy_num) + 0.5, 1))
    
    # Calculate normalized histogram values
    hist_data$counts <- hist(copy_num, breaks = seq(-0.5, max(copy_num) + 0.5, 1), plot = FALSE)$counts / length(copy_num)
    
    # Create a data frame for the Poisson distribution
    poisson_data <- data.frame(x = hist_data$breaks, y = dpois(hist_data$breaks - 0.5, lambda = mean_copy_num))
    
    # Create the plot
    plot <- ggplot() +
      geom_bar(data = poisson_data, aes(x = x, y = y), stat = "identity", width = 1, fill = "lightblue") +
      geom_bar(data = hist_data, aes(x = breaks, y = counts), stat = "identity", width = 0.5, fill = "black", color = 'black') +
      labs(title = plot_title,
           x = "Copy Number",
           y = "Proportion") +
      theme_minimal() +
      xlim(c(0, max_num + 1)) +
      theme(legend.title = element_text(size = 30),
            legend.text = element_text(size = 24),
            plot.title = element_text(size = 48),
            axis.title = element_text(size = 40),
            axis.text = element_text(color = "black", size = 30),
            legend.position = c(10.845, 0.6))
    
    if(mean){
      plot = plot + geom_vline(xintercept=mean_copy_num,lwd=1,colour="red")

    }
  }
  else {
    plot <- ggplot() +
      labs(title = plot_title,
           x = "Copy Number",
           y = "Proportion") +
      theme_minimal()
  }
  
  plot
}
get_pois_plot <- function(datadir, plot_title = "", max_num = -1, start_time = -Inf, end_time = Inf, mean = T){
  copy_data = read.csv(file.path(datadir, "copy_nums_.csv")) %>%
    filter(time >= start_time, time <= end_time)
  if(max_num == -1) {
    max_num = max(copy_data$copy_num)
  }
  get_pois_plot_df(copy_data, plot_title, max_num, mean)
}

get_copy_plot <- function(datadir, start_time=0, end_time=Inf, plot_title="") {
  plot_data <- read.csv(file.path(datadir, "copy_nums_.csv"))
  plot_data <- filter(plot_data, time >= start_time, time <= end_time)
  ggplot(plot_data, aes(x = time, y = copy_num)) +
    geom_point() +
    labs(title = plot_title,
         x = "Time (s)",
         y = "Copy Number")
}

get_elong_plot <- function(datadir, start_time=0, end_time=Inf, plot_title="") {
  plot_data <- read.table(file.path(datadir, "mean_properties.txt"), 
                          header = TRUE, 
                          sep = "")
  
  plot_data <- mutate(plot_data,
                      elong_time = transcripts_nb * mean_elong_time - (lag(transcripts_nb, default = 0) * lag(mean_elong_time, default = 0)))
  
  
  
  plot_data <- filter(plot_data, time >= start_time, time <= end_time)
  
  plot_data <- mutate(plot_data, time = time/60)
  ggplot(plot_data, aes(x = time, y = elong_time)) +
    geom_point() +
    labs(title = plot_title,
         x = "Time (min)",
         y = "Elongation time (s)") +
    theme_bw() +
    theme(legend.title = element_text(size = 24),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 24),
          axis.title = element_text(size = 24),
          axis.text = element_text(color = "black", size = 18),
          legend.position = c(0.845, 0.5))
}

calc_rate <- function(datadir, start_time = 0, end_time = 3600) {
  transcript_data <- read.table(file.path(datadir, "mean_properties.txt"), 
                          header = TRUE, 
                          sep = "")
  
  (count(filter(transcript_data, time >= start_time, time <= end_time)) / (end_time - start_time))
}

import_data <- function(datadir, start_time = 0, end_time = Inf) {
  subdirs = list.files(path = datadir)
  def_lansG = 0.0001
  def_LasG = 0
  def_lansT = 0.0001
  def_LasT = 0
  def_up = 0
  def_down = 0
  def_gene = 0
  def_barrier = 0
  def_k = 0
  def_gsa = -0.06166666666667
  def_gb = 1
  
  default_row = data.frame(transcripts_nb = 0,
                           time = NaN,
                           prod_rate = 0,
                           mean_prod_time = NaN,
                           mean_bind_time = NaN,
                           mean_ocf_time = NaN,
                           mean_esc_time = NaN,
                           mean_init_time = NaN,
                           mean_elong_time = NaN)
  
  data = data.frame()
  
  for (dir in subdirs){
    #get datas
    path = file.path(datadir, dir, 'mean_properties.txt')
    subdata = read.table(path, header = TRUE, sep = "") %>%
      filter(time >= start_time & time <= end_time)
    
    #assign default to empty data
    if(nrow(subdata) == 0){
      subdata = rbind(subdata, default_row)  
    }
    
    #set defaults
    subdata = subdata %>% 
      mutate(lansT = def_lansT,  
             LasT = def_LasT, 
             lansG = def_lansG, 
             LasG = def_LasG, 
             up = def_up,
             down = def_down,
             gene = def_gene,
             barrier = def_barrier,
             kb = def_k,
             ko = def_k,
             gsa = def_gsa,
             gb = def_gb)
    
    #set true values
    info = strsplit(dir, '_')
    info = info[[1]]
    for(param in info) {
      param = strsplit(param, '-')
      if(param[[1]][1] == "lansG") {
        subdata = mutate(subdata, lansG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasG") {
        subdata = mutate(subdata, LasG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "lansT") {
        subdata = mutate(subdata, lansT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasT") {
        subdata = mutate(subdata, LasT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "up") {
        subdata = mutate(subdata, up = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "down") {
        subdata = mutate(subdata, down = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gene") {
        subdata = mutate(subdata, gene = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "k") {
        subdata = mutate(subdata, kb = as.numeric(param[[1]][2]), ko = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gsa") {
        subdata = mutate(subdata, gsa = -as.numeric(param[[1]][2])) #note negative
      }
      else if(param[[1]][1] == "gb") {
        subdata = mutate(subdata, gb = as.numeric(param[[1]][2]))
      }
    }
    
    if(sum(dim(data)) == 0) {
      data = subdata
    } else {
      data = rbind(data, subdata)
    }
    
  } #data import
  
  data
}

import_prom_data <- function(datadir, start_time = 0, end_time = Inf) {
  subdirs = list.files(path = datadir)
  
  #defaults
  def_lansG = 0.0001
  def_LasG = 0
  def_lansT = 0.0001
  def_LasT = 0
  def_up = 0
  def_down = 0
  def_gene = 0
  def_barrier = 0
  def_k = 0
  def_gsa = -0.06166666666667
  def_gb = 1
  
  
  data = data.frame()
  
  for (dir in subdirs){
    print(paste0("Strating dir: ", dir))
    #get data
    path = file.path(datadir, dir, 'traj_promoter.txt')
    subdata = read.table(path, header = TRUE, sep = "") %>%
      filter(time >= start_time, time <= end_time) %>%
      summarize(up_mean = mean(sigma_up), 
                gene_mean = mean(sigma_gene), 
                down_mean = mean(sigma_down), 
                up_elong = mean(sigma_up >= -0.061666666666666667), 
                down_elong = mean(sigma_down <= 0.061666666666666667), 
                open = mean(sigma_up <= -0.05), 
                up_valid = mean(sigma_up >= -0.061666666666666667 & sigma_up <= -0.05), 
                valid = mean(sigma_up >= -0.061666666666666667 & 
                            sigma_up <= -0.05 &
                              sigma_down <= 0.061666666666666667),
                off_proportion = mean(t_upRNAPelong == 0),
                off_time = mean(with(rle(t_upRNAPelong == 0), lengths[values])))
    
    #set defaults
    subdata = subdata %>% 
      mutate(lansT = def_lansT,  
             LasT = def_LasT, 
             lansG = def_lansG, 
             LasG = def_LasG, 
             up = def_up,
             down = def_down,
             gene = def_gene,
             barrier = def_barrier,
             kb = def_k,
             ko = def_k,
             gsa = def_gsa)
    
    #set true values
    info = strsplit(dir, '_')
    info = info[[1]]
    for(param in info) {
      param = strsplit(param, '-')
      if(param[[1]][1] == "lansG") {
        subdata = mutate(subdata, lansG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasG") {
        subdata = mutate(subdata, LasG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "lansT") {
        subdata = mutate(subdata, lansT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasT") {
        subdata = mutate(subdata, LasT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "up") {
        subdata = mutate(subdata, up = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "down") {
        subdata = mutate(subdata, down = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gene") {
        subdata = mutate(subdata, gene = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "k") {
        subdata = mutate(subdata, kb = as.numeric(param[[1]][2]), ko = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gsa") {
        subdata = mutate(subdata, gsa = -as.numeric(param[[1]][2])) #note negative
      }  
      else if(param[[1]][1] == "gb") {
        subdata = mutate(subdata, gb = as.numeric(param[[1]][2]))
      }
    }
    
    if(sum(dim(data)) == 0) {
      data = subdata
    } else {
      data = rbind(data, subdata)
    }
    
    print(paste0("Finished dir: ", dir))
  } #data import
  
  data
}

import_steric_data <- function(datadir, start_time = 0, end_time = Inf, ki = 0.2, si = -0.04, bi = 0.005) {
  subdirs = list.files(path = datadir)
  
  #defaults
  def_lansG = 0.0001
  def_LasG = 0
  def_lansT = 0.0001
  def_LasT = 0
  def_up = 0
  def_down = 0
  def_gene = 0
  def_barrier = 0
  def_k = 0
  def_gsa = -0.06166666666667
  def_gb = 1
  
  
  data = data.frame()
  
  for (dir in subdirs){
    print(paste0("Starting dir: ", dir))
    #get data
    path_prom = file.path(datadir, dir, 'traj_promoter.txt')
    path_mean = file.path(datadir, dir, 'mean_properties.txt')
    
    mean_data = read.table(path_mean, header = TRUE, sep = "") %>%
      filter(time >= start_time, time <= end_time)
    
    if (nrow(mean_data) < 2){
      rate = 0
    } else {
      rate = fix_mean(mean_data$prod_rate[which.max(mean_data$transcripts_nb)], mean_data$prod_rate[which.min(mean_data$transcripts_nb)], max(mean_data$transcripts_nb), min(mean_data$transcripts_nb))
    }
    
    
    subdata = read.table(path_prom, header = TRUE, sep = "") %>%
      filter(time >= start_time, time <= end_time) %>%
      mutate(kiprime = ki / (1 + exp((sigma_up - si)/bi))) %>%
      summarize(kiprime_mean = mean(kiprime))
    
    #set defaults
    subdata = subdata %>% 
      mutate(lansT = def_lansT,  
             LasT = def_LasT, 
             lansG = def_lansG, 
             LasG = def_LasG, 
             up = def_up,
             down = def_down,
             gene = def_gene,
             barrier = def_barrier,
             kb = def_k,
             ko = def_k,
             gsa = def_gsa,
             prod_rate = rate,
             free_prom = prod_rate/kiprime_mean)
    
    
    
    #set true values
    info = strsplit(dir, '_')
    info = info[[1]]
    for(param in info) {
      param = strsplit(param, '-')
      if(param[[1]][1] == "lansG") {
        subdata = mutate(subdata, lansG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasG") {
        subdata = mutate(subdata, LasG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "lansT") {
        subdata = mutate(subdata, lansT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasT") {
        subdata = mutate(subdata, LasT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "up") {
        subdata = mutate(subdata, up = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "down") {
        subdata = mutate(subdata, down = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gene") {
        subdata = mutate(subdata, gene = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "k") {
        subdata = mutate(subdata, kb = as.numeric(param[[1]][2]), ko = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gsa") {
        subdata = mutate(subdata, gsa = -as.numeric(param[[1]][2])) #note negative
      }  
      else if(param[[1]][1] == "gb") {
        subdata = mutate(subdata, gb = as.numeric(param[[1]][2]))
      }
    }
    
    if(sum(dim(data)) == 0) {
      data = subdata
    } else {
      data = rbind(data, subdata)
    }
    
    print(paste0("Finished dir: ", dir))
  } #data import
  
  data
}

import_copy_data <- function(datadir, start_time = 0, end_time = Inf) {
  subdirs = list.files(path = datadir)
  
  #defaults
  def_lansG = 0.0001
  def_LasG = 0
  def_lansT = 0.0001
  def_LasT = 0
  def_up = 0
  def_down = 0
  def_gene = 0
  def_barrier = 0
  def_k = 0
  def_gsa = -0.06166666666667
  def_gb = 1
  
  
  data = data.frame()
  
  for (dir in subdirs){
    print(paste0("Strating dir: ", dir))
    #get data
    path = file.path(datadir, dir, 'copy_nums_.csv')
    subdata = read.csv(path, header = TRUE, sep = ",") %>%
      filter(time >= start_time, time <= end_time) %>%
      summarize(mean = mean(copy_num),
                var = var(copy_num),
                zero = mean(copy_num == 0),
                mean_off = 0.04 * mean(with(rle(copy_num == 0), lengths[values]))) %>%
      mutate(std = sqrt(var),
             fano = var/mean,
             noise = var/(mean**2))
    
    #set defaults
    subdata = subdata %>% 
      mutate(lansT = def_lansT,  
             LasT = def_LasT, 
             lansG = def_lansG, 
             LasG = def_LasG, 
             up = def_up,
             down = def_down,
             gene = def_gene,
             barrier = def_barrier,
             kb = def_k,
             ko = def_k,
             gsa = def_gsa,
             gb = def_gb)
    
    #set true values
    info = strsplit(dir, '_')
    info = info[[1]]
    for(param in info) {
      param = strsplit(param, '-')
      if(param[[1]][1] == "lansG") {
        subdata = mutate(subdata, lansG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasG") {
        subdata = mutate(subdata, LasG = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "lansT") {
        subdata = mutate(subdata, lansT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "LasT") {
        subdata = mutate(subdata, LasT = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "up") {
        subdata = mutate(subdata, up = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "down") {
        subdata = mutate(subdata, down = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gene") {
        subdata = mutate(subdata, gene = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "k") {
        subdata = mutate(subdata, kb = as.numeric(param[[1]][2]), ko = as.numeric(param[[1]][2]))
      }
      else if(param[[1]][1] == "gsa") {
        subdata = mutate(subdata, gsa = -as.numeric(param[[1]][2])) #note negative
      }
      else if(param[[1]][1] == "gb") {
        subdata = mutate(subdata, gb = as.numeric(param[[1]][2]))
      }
    }
    
    if(sum(dim(data)) == 0) {
      data = subdata
    } else {
      data = rbind(data, subdata)
    }
    
    print(paste0("Finished dir: ", dir))
  } #data import
  
  data
}

summarize_data <- function(data_df) {
  grouped_data = data_df %>% 
    group_by(lansT, LasT, lansG, LasG, up, gene, down, kb, gsa, gb)
  
  summary_data = grouped_data %>%
    summarize(mean_rate = fix_mean(prod_rate[which.max(transcripts_nb)], prod_rate[which.min(transcripts_nb)], max(transcripts_nb), min(transcripts_nb)),
              mean_elong = fix_mean(mean_elong_time[which.max(transcripts_nb)], mean_elong_time[which.min(transcripts_nb)], max(transcripts_nb), min(transcripts_nb)),
              n = n(),
              last_time = max(time)) 
              #remaining_time = end_time - last_time) #%>%
    #mutate(lansT = as.factor(lansT), lansG = as.factor(lansG))
  
  summary_data
}


get_transcription_plot = function(plot_data, xvar, yvar = "mean_rate", plot_title="", x_label = "", y_label = "Transcription Rate (s^-1)"){
  
  plot = ggplot(data = plot_data, aes_string(x = xvar, y = yvar)) +
    geom_point(size = 4, color = 'red') +
    ggtitle(plot_title) +
    ylab(y_label) +
    xlab(x_label)  +
    theme_bw() +
    theme(legend.title = element_text(size = 24),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 24),
          axis.title = element_text(size = 24),
          axis.text = element_text(color = "black", size = 10),
          legend.position = c(0.845, 0.6)) #+
    #ylim(0, 0.25)
  
  plot
}

get_normalized_transcription_plot = function(plot_data, xvar, plot_title=""){
  
  plot_data <- filter(plot_data, !is.na(mean_rate))
  
  plot_data$mean_rate <- scale(plot_data$mean_rate, center = FALSE, scale = max(plot_data$mean_rate))
  
  plot = ggplot(data = plot_data, aes_string(x = xvar, y = 'mean_rate')) +
    geom_point(size = 4, color = 'red') +
    ggtitle(plot_title) +
    ylab('Relative Transcription Rate') +
    xlab('Nonspecific Gyrase Activity') +
    ylim(0, 1)
  
  plot
}

get_transcription_grid <- function(plot_data, 
                                   xvar = "lansG", 
                                   yvar = "lansT", 
                                   fillvar = "mean_rate",
                                   plot_title="", 
                                   x_label = "Gyrase Activity (Lk/bp/s)", 
                                   y_label = "TopoI Activity (Lk/bp/s)",
                                   fill_label = "Transcription Rate (1/s)",
                                   gyr_scale = 1,
                                   topoI_scale = 1,
                                   low_color = "black",
                                   high_color = "white",
                                   fill_min = -Inf,
                                   fill_max = Inf,
                                   opt_x = 4.5e-5,
                                   opt_y = 4.5e-5,
                                   legend=TRUE,
                                   opt=TRUE){
  
  plot_data = mutate(plot_data, lansG = gyr_scale * lansG, lansT = topoI_scale * lansT) #%>%
  #mutate(lansG = as.factor(lansG), lansT = as.factor(lansT))
  
  if(fill_min == -Inf) {
    fill_min = min(plot_data[[fillvar]])
  }
  if(fill_max == Inf) {
    fill_max = max(plot_data[[fillvar]])
  }
  
  plot = ggplot(data = plot_data, aes_string(x = xvar, y = yvar, fill=fillvar)) +
    geom_tile() +
    #geom_hline(yintercept = opt_y * topo_scale, linetype = "dashed", color = "red", size = 1) +
    #geom_vline(xintercept = opt_x * topo_scale, linetype = "dashed", color = "red", size = 1) +
    scale_fill_gradient(low=low_color, high=high_color, limits = c(fill_min, fill_max), oob=squish) +
    ggtitle(plot_title) +
    xlab(x_label) +
    ylab(y_label) +
    labs(fill = fill_label)  +
    theme_bw() +
    theme(legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 24),
          axis.title = element_text(size = 24),
          axis.text = element_text(color = "black", size = 20))
  
  if(!legend) {
    plot = plot + theme(legend.position="none")
  }
  
  if(opt) {
    plot = plot +
      geom_hline(yintercept = opt_y * topoI_scale, color = "cyan", size = 0.8) +
      geom_vline(xintercept = opt_x * gyr_scale, color = "cyan", size = 0.8)
  }
  plot
}

get_transcription_grid_tricolor <- function(plot_data, 
                                   xvar = "lansG", 
                                   yvar = "lansT", 
                                   fillvar = "mean_rate",
                                   plot_title="", 
                                   x_label = "Gyrase Activity (Lk/bp/s)", 
                                   y_label = "TopoI Activity (Lk/bp/s)",
                                   fill_label = "Transcription Rate (1/s)",
                                   gyr_scale = 1,
                                   topoI_scale = 1,
                                   low_color = "blue",
                                   mid_color = "white",
                                   high_color = "red",
                                   fill_min = -Inf,
                                   fill_mid = 0,
                                   fill_max = Inf,
                                   opt_x = 4.5e-5,
                                   opt_y = 4.5e-5,
                                   legend=TRUE,
                                   opt=TRUE){
  
  plot_data = mutate(plot_data, lansG = gyr_scale * lansG, lansT = topoI_scale * lansT) #%>%
  #mutate(lansG = as.factor(lansG), lansT = as.factor(lansT))
  
  if(fill_min == -Inf) {
    fill_min = min(plot_data[[fillvar]])
  }
  if(fill_max == Inf) {
    fill_max = max(plot_data[[fillvar]])
  }
  
  plot = ggplot(data = plot_data, aes_string(x = xvar, y = yvar, fill=fillvar)) +
    geom_tile() +
    #geom_hline(yintercept = opt_y * topo_scale, linetype = "dashed", color = "red", size = 1) +
    #geom_vline(xintercept = opt_x * topo_scale, linetype = "dashed", color = "red", size = 1) +
    scale_fill_gradient2(low = low_color, 
                         mid = mid_color, 
                         high = high_color, 
                         midpoint = fill_mid, 
                         limits = c(fill_min, fill_max), 
                         oob = squish) +
    ggtitle(plot_title) +
    xlab(x_label) +
    ylab(y_label) +
    labs(fill = fill_label)  +
    theme_bw() +
    theme(legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 24),
          axis.title = element_text(size = 24),
          axis.text = element_text(color = "black", size = 20))
  
  if(!legend) {
    plot = plot + theme(legend.position="none")
  }
  
  if(opt) {
    plot = plot +
      geom_hline(yintercept = opt_y * topoI_scale, color = "cyan", size = 0.8) +
      geom_vline(xintercept = opt_x * gyr_scale, color = "cyan", size = 0.8)
  }
  plot
}
