# Load necessary libraries for spatial transcriptomics visualization and metric calculation
library(Matrix)
library(Seurat)
library(magrittr)
library(SeuratDisk)
library(anndata)  
library(ggplot2)
library(ggpubr)
require(cowplot)
library(patchwork)
library(tidyverse)
library(aricode)
library(reticulate)

# Function to generate spatial plots and compute clustering metrics
generate_spatial_plots <- function(data_file, original_data, output_file) {
  # Parameters:
  # data_file: domain_plot.txt, is the path to store the results of the different methods
  # original_data: Path to the original data in .h5ad format to be used as a reference for spatial coordinates
  # output_file: Directory path where the output plots and metrics will be saved
  
  # Specify the Python environment to use (required by SeuratDisk for Python-R integration)
  #use_python("/home/dongkj/anaconda3/envs/MultiSpatial/bin/python", required=TRUE)
  
  # Read the spatial data file which contains information about the models and their data sources
  spatial_data_file <- read.table(data_file, header = FALSE)
  
  # Define a set of colors to be used for the plots (for predicted domains)
  color2 <- c("#E89C9A","#ACCBDE","#CF3732","#BAD691","#002FA7","#008C8C","#81D8CF","#B05923","#900021","#E60000","#FBD26A","#E85827","#432913")
  
  # Initialize metadata for storing information about the spatial and model data
  data <- read_h5ad(original_data)
  metadata <- matrix(,dim(data)[1],)
  data_meta <- data$obs
  metadata <- cbind(metadata, data_meta[,match(c("barcode","X","Y","slices","original_domain"), colnames(data_meta))]) 
  
  # Loop through each model's data file, extract metadata, and bind it together for further analysis
  for(i in 1:nrow(spatial_data_file)){
    model <- unlist(strsplit(spatial_data_file[i, 1], split = ":"))[1]
    ex <- unlist(strsplit(spatial_data_file[i, 1], split = ":"))[2]
    data <- read_h5ad(ex)
    
    # For subsequent models, match barcodes and add the predicted domains to the metadata
    data_meta <- data$obs
    index <- match(metadata[,2], data_meta$barcode)
    predicted_domain <- as.data.frame(as.character(data_meta[index,]$predicted_domain))
    colnames(predicted_domain) <- model
    metadata <- cbind(metadata, predicted_domain)
  }
  
  # Remove the first column of metadata (which was just for initialization)
  metadata <- metadata[, -1]
  metadata <- as.data.frame(metadata)
  
  # Read the original spatial data (used for spatial coordinates)
  spatial <- read_h5ad(original_data)
  
  # Match metadata barcodes with spatial barcodes and extract spatial coordinates
  index <- match(metadata$barcode, spatial$obs$barcode)
  spatial <- spatial$obsm$spatial[index,]
  metadata$X <- as.numeric(spatial[,1])
  metadata$Y <- as.numeric(spatial[,2])
  metadata <- na.omit(metadata)  # Remove any rows with missing values
  
  # Write the final metadata to a CSV file
  write.table(metadata, paste(output_file, "metaPredictedRe.csv", sep = "/"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
  
  # Initialize empty vectors to store ARI and NMI values for each model
  model_names <- colnames(metadata)[-c(1:5)]
  model_ari <- c()
  model_nmi <- c()
  
  # Compute ARI and NMI for each model compared to the ground truth "original_domain"
  for(m in model_names){
    o_domain <- metadata$original_domain
    p_domain <- metadata[,match(m, colnames(metadata))]
    
    # Calculate Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI)
    ari <- ARI(as.character(o_domain), as.character(p_domain))
    nmi <- NMI(as.character(o_domain), as.character(p_domain))
    model_ari <- c(model_ari, ari)
    model_nmi <- c(model_nmi, nmi)
  }
  
  # Combine the ARI and NMI results into a dataframe and write to CSV
  model_metric <- cbind(model_names, cbind(model_ari, model_nmi))
  model_metric <- as.data.frame(model_metric)
  model_metric$model_ari <- as.numeric(model_metric$model_ari)
  model_metric$model_nmi <- as.numeric(model_metric$model_nmi)
  write.table(model_metric, paste(output_file, "model_metric.csv", sep = "/"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
  
  # Create bar plots for ARI and NMI values for each model
  ari_plot <- ggplot(model_metric, aes(x = model_names, y = model_ari, fill = model_names)) +
    geom_bar(stat = "identity") +
    labs(title = "ARI Values", x = "Model", y = "ARI") +
    theme_minimal() + 
    theme(plot.title = element_text(size = 3),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 3),
          legend.title = element_text(size = 3),
          axis.title = element_text(size = 3),
          axis.text = element_text(size = 3),
          axis.text.x = element_text(angle = 30, hjust = 1)) 
  
  nmi_plot <- ggplot(model_metric, aes(x = model_names, y = model_nmi, fill = model_names)) +
    geom_bar(stat = "identity") +
    labs(title = "NMI Values", x = "Model", y = "NMI") +
    theme_minimal() + 
    theme(plot.title = element_text(size = 3),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 3),
          legend.title = element_text(size = 3),
          axis.title = element_text(size = 3),
          axis.text = element_text(size = 3),
          axis.text.x = element_text(angle = 30, hjust = 1)) 
  
  # Combine ARI and NMI plots and save as a PDF file
  metric_plot <- ari_plot + nmi_plot
  ggsave(paste(output_file, "metric.pdf", sep = "/"), metric_plot, width = 10, height = 5, units = "cm", limitsize = FALSE)
  
  # Generate spatial plots for each unique slice in the metadata
  sample_names <- unique(metadata$slices)
  spatial_plots <- lapply(sample_names, function(snm) {
    sample_data <- metadata[metadata$slices == snm,]
    
    # Initialize a list to store individual plots for each model and slice
    plot_list <- list()
    plot_list[[1]] <- ggplot(sample_data, aes(x = X, y = Y, col = original_domain)) +
      geom_point(size = 0.1) + 
      theme_classic() + 
      theme(
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
      labs(title = "Ground truth", y = paste("slices", as.character(snm))) +
      coord_equal() + 
      theme(plot.title = element_text(size = 6),
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.1, 'cm'),
            axis.title.y = element_text(size = 6),
            legend.title = element_text(size = 6)) 
    
    # Loop through each model and add predicted domain plots
    for(z in 1:length(model_names)){
      m <- model_names[z]
      ari <- model_metric[z, 2]
      
      df <- cbind.data.frame(x_index = sample_data$X, y_index = sample_data$Y, predicted_domain = sample_data[, model_names[z]])
      df$predicted_domain <- as.character(df$predicted_domain)
      
      plot_list[[z+1]] <- ggplot(df, aes(x = x_index, y = y_index, col = predicted_domain)) +
        geom_point(size = 0.1) + 
        theme_classic() + 
        theme(
          legend.position = "right",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
        labs(title =  sprintf("%s - ARI: %s", m, round(ari, 3))) +
        coord_equal() +
        theme(plot.title = element_text(size = 6),
              legend.key.size = unit(0.1, 'cm'),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 6)) 
    }
    
    # Combine all plots for a single slice
    p <- wrap_plots(plot_list, ncol = length(plot_list))
  })
  
  # Combine spatial plots for all slices and save to PDF
  p <- plot_grid(plotlist = spatial_plots, ncol = 1, byrow = TRUE)
  ggsave(paste(output_file, "domain.pdf", sep = "/"), p, width = 40, height = 25, units = "cm", limitsize = FALSE)
}

# Main script to run the function with arguments passed from the command line
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
original_data <- args[2]
output_file <- args[3]

# Call the function to generate spatial plots and metrics
generate_spatial_plots(data_file, original_data, output_file)
