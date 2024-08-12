library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

get_copy_nums <- function(prod_times, end_time, dt = 0.04, RNA_hl = 120) {
  
  copy_df = data.frame(time = seq(0, end_time, dt), copy_num = 0)
  
  if(length(prod_times) == 0) {
    return(copy_df)
  }
  
  last_index = end_time / dt + 1
  rna_df = data.frame(synth_index = 1 + prod_times / dt, lifetime = rgeom(length(prod_times), dt / RNA_hl)) %>%
    mutate(decay_index = pmin(synth_index + lifetime, end_time / dt + 1))
  
  for (i in 1:nrow(rna_df)) {
    copy_df$copy_num[rna_df$synth_index[i]:rna_df$decay_index[i]] = copy_df$copy_num[rna_df$synth_index[i]:rna_df$decay_index[i]] + 1
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
    end_time = 79999.96
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

save_copy_num("data")

