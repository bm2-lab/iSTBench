library(Seurat)
library(magrittr)
library(SeuratDisk)
library(anndata)
library(ggplot2)
library(dplyr)
library(getopt)
library(tidyr)
library(reticulate)

# Specify the Python environment to use (required by SeuratDisk for Python-R integration)
# use_python("/home/dongkj/anaconda3/envs/MultiSpatial/bin/python", required=TRUE)


# Command-line argument definition
# The following spec matrix defines the input arguments for the script: 
spec <- matrix(
  c(
    'input_file', 'i', 1, 'character', 'Input file path',   
    "slices", "s", 1, "character", "Slices"                 
  ),
  byrow = TRUE,
  ncol = 5
)

# Parse the command-line arguments
opt <- getopt(spec)

# Check if required arguments are provided
# If any of the required arguments 'input_file' or 'slices' is missing, stop execution and display an error message.
if (is.null(opt$input_file) | is.null(opt$slices)) {
  stop("Error: input_file parameter is required.")
}

# Assign parsed arguments to variables
input_file <- opt$input_file
slices <- opt$slices

# Function to analyze spatial data from multiple slices
# This function computes Accuracy and Ratio across slices and visualizes them.
analyze_spatial_data <- function(input_file, slices) {
  # Define the output directory and file names for the results
  output_file <- paste(input_file, "Metric", sep = "/")
  
  # List all files in the input directory, excluding irrelevant files
  fl <- list.files(input_file)
  fl <- fl[!grepl("Metric", fl)]  # Exclude files with 'Metric' in their names
  fl <- fl[!grepl("Dr.A", fl)]  # Exclude Dr.A related files
  
  # Initialize empty data frames to store the results
  domain_match_rate_all <- data.frame()
  domain_match_ratio_all <- data.frame()
  
  # Iterate over each file and process it
  for (i in 1:(length(fl) + 1)) {
    if(i == (length(fl) + 1)) {
      # Handle the "unCorrect" file
      f <- "unCorrect"
      f_file <- paste(input_file, "CellCharter", "Spatial_correct_data.h5ad", sep = "/")
      data <- read_h5ad(f_file)  # Read the .h5ad file
      
      metadata <- data$obs  # Extract metadata
      spatial <- metadata[, match(c("X", "Y"), colnames(metadata))]  # Extract spatial coordinates
    } else {
      # Process the other files
      f <- fl[i]
      f_file <- paste(input_file, f, "Spatial_correct_data.h5ad", sep = "/")
      data <- read_h5ad(f_file)  # Read the .h5ad file
      
      metadata <- data$obs  # Extract metadata
      spatial <- data$obsm$spatial  # Extract spatial coordinates
    }
    
    slices <- unlist(strsplit(slices, split = ","))
    
    # Initialize vectors to store matching rates and ratios for each slice pair
    domain_match_rate <- c()
    domain_match_ratio <- c()
    
    # Compute domain match rate and ratio for each consecutive slice pair
    for(s in 1:(length(slices) - 1)) {
      s1_metadata <- metadata[metadata$slice_name == slices[s], ]
      s2_metadata <- metadata[metadata$slice_name == slices[s + 1], ]
      
      s1_spatial <- spatial[which(metadata$slice_name == slices[s]), ]
      s2_spatial <- spatial[which(metadata$slice_name == slices[s + 1]), ]
      
      # Calculate distance matrix between spatial coordinates of slice pairs
      distance_matrix <- as.matrix(dist(rbind(s1_spatial, s2_spatial)))
      
      # Extract the relevant part of the distance matrix
      s1_length <- nrow(s1_spatial)
      s2_length <- nrow(s2_spatial)
      distance_matrix <- distance_matrix[1:s1_length, (s1_length + 1):(s1_length + s2_length)]
      min_distance_indices <- apply(distance_matrix, 1, which.min)
      
      s1_domain <- s1_metadata$original_domain
      s2_domain <- s2_metadata[min_distance_indices, ]$original_domain
      
      # Calculate domain match rate and ratio
      same_rate <- sum(s1_domain == s2_domain) / length(s1_domain)
      match_ratio <- min(length(s1_domain), length(s2_domain)) / length(unique(min_distance_indices))
      
      domain_match_rate <- c(domain_match_rate, same_rate)
      domain_match_ratio <- c(domain_match_ratio, match_ratio)
    }
    
    # Append results for this file to the overall data frames
    domain_match_rate <- c(domain_match_rate, f)
    domain_match_ratio <- c(domain_match_ratio, f)
    
    domain_match_rate_all <- rbind(domain_match_rate_all, domain_match_rate)
    domain_match_ratio_all <- rbind(domain_match_ratio_all, domain_match_ratio)
  }
  
  # Define column names for the result data frames
  colnames(domain_match_rate_all) <- c(paste(paste("s", (1: (length(slices) - 1)), sep = ""), paste("s", (2: (length(slices))), sep = ""), sep = "-"), "model")
  colnames(domain_match_ratio_all) <- colnames(domain_match_rate_all)
  
  # Convert domain match rate and ratio data to numeric
  domain_match_rate_all <- as.data.frame(domain_match_rate_all)
  model <- as.character(domain_match_rate_all$model)
  domain_match_rate_all <- apply(as.data.frame(domain_match_rate_all[,-which(colnames(domain_match_rate_all) == "model")]), 2, as.numeric)
  domain_match_rate_all <- as.data.frame(domain_match_rate_all)
  domain_match_rate_all$model <- model
  write.table(domain_match_rate_all, paste(output_file, "Accuracy.csv", sep = "/"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
  
  domain_match_ratio_all <- as.data.frame(domain_match_ratio_all)
  model <- as.character(domain_match_ratio_all$model)
  domain_match_ratio_all <- apply(as.data.frame(domain_match_ratio_all[,-which(colnames(domain_match_ratio_all) == "model")]), 2, as.numeric)
  domain_match_ratio_all <- as.data.frame(domain_match_ratio_all)
  domain_match_ratio_all$model <- model
  write.table(domain_match_ratio_all, paste(output_file, "Ratio.csv", sep = "/"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
  
  # Visualization: Generate bar and box plots for Accuracy and ratio
  for(i in 1:2) {
    if(i == 1) {
      # Bar plot for Accuracy
      df_long <- domain_match_rate_all %>%
        pivot_longer(cols = -model, names_to = "metric", values_to = "value")
      
      df_long$model <- factor(df_long$model, levels = c(c("Real,Slices_combind_data"), setdiff(unique(df_long$model), c("Real,Slices_combind_data"))))
      
      p1 <- ggplot(df_long, aes(x = metric, y = value, fill = model)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
        labs(x = "Slices", y = "Accuracy") +
        theme_minimal() +
        theme(legend.position = "top")
      ggsave(paste(output_file, "Accuracy_bar.pdf", sep = "/"), p1, height = 15, width = 35, units = "cm")
      
      # Box plot for Accuracy
      p2 <- ggplot(df_long, aes(x = model, y = value, fill = model)) +
        geom_boxplot() +
        labs(x = "Model", y = "Accuracy") +
        theme_minimal() +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
      ggsave(paste(output_file, "Accuracy_boxplot.pdf", sep = "/"), p2, height = 10, width = 20, units = "cm")
    } else {
      # Bar plot for Ratio
      df_long <- domain_match_ratio_all %>%
        pivot_longer(cols = -model, names_to = "metric", values_to = "value")
      
      df_long$model <- factor(df_long$model, levels = c(c("Real,Slices_combind_data"), setdiff(unique(df_long$model), c("Real,Slices_combind_data"))))
      
      p1 <- ggplot(df_long, aes(x = metric, y = value, fill = model)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
        labs(x = "Slices", y = "Ratio") +
        theme_minimal() +
        theme(legend.position = "top")
      ggsave(paste(output_file, "Ratio_bar.pdf", sep = "/"), p1, height = 15, width = 35, units = "cm")
      
      # Box plot for Ratio
      p2 <- ggplot(df_long, aes(x = model, y = value, fill = model)) +
        geom_boxplot() +
        labs(x = "Model", y = "Ratio") +
        theme_minimal() +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
      ggsave(paste(output_file, "Ratio_boxplot.pdf", sep = "/"), p2, height = 10, width = 20, units = "cm")
    }
  }
}

# Run the analysis function
analyze_spatial_data(input_file, slices)
