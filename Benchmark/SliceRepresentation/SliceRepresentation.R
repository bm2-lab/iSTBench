# Import necessary libraries
library(Seurat)
library(magrittr)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(tidyr)
library(stats)
library(umap)
library(ggplot2)
library(patchwork)
library(aricode)
library(leidenAlg)
library(igraph)
library(reticulate)
library(getopt)

#find_optimal_resolution Function ----
# Function to find optimal resolution for clustering based on the Leiden algorithm
# This function aims to adjust the resolution parameter of the Leiden algorithm to achieve a target number of clusters.
# It uses a binary search approach to gradually adjust the resolution until the number of clusters meets the target.
find_optimal_resolution <- function(graph, target_clusters, tol = 1e-3, max_iter = 100) {
  lower_res <- 0.0
  upper_res <- 3.0
  optimal_res <- (lower_res + upper_res) / 2
  iter <- 0
  
  while (iter < max_iter) {
    iter <- iter + 1
    leiden_result <- leiden.community(graph, resolution = optimal_res)
    num_clusters <- length(unique(leiden_result[["membership"]]))
    
    if (abs(num_clusters - target_clusters) <= tol) {
      return(optimal_res)
    } else if (num_clusters < target_clusters) {
      lower_res <- optimal_res
    } else {
      upper_res <- optimal_res
    }
    
    optimal_res <- (lower_res + upper_res) / 2
  }
  
  return(optimal_res)
}

#slices_clustering Function ----
# Function to perform clustering of slices based on a given metadata and predicted domain
# This function calculates the abundance of each domain within each slice and performs clustering using either hierarchical clustering or Leiden community detection.
# It returns the clustering results along with various performance metrics, such as ARI and NMI, and a PCA plot.
slices_clustering <- function(metadata, predicted_domain, dist_method, cluster_method, cluster_number, slices_class, seed = 123){
  
  # Step 1: Calculate abundance matrix for each slice
  count_df <- metadata %>%
    group_by(slices, !!ensym(predicted_domain)) %>%
    summarise(count = n()) %>%
    ungroup()
  
  total_counts <- metadata %>%
    group_by(slices) %>%
    summarise(total_count = n())
  
  merged_df <- count_df %>%
    left_join(total_counts, by = "slices")
  
  abundance_df <- merged_df %>%
    mutate(proportion = count / total_count) %>%
    select(slices, !!ensym(predicted_domain), proportion)
  
  abundance_matrix <- abundance_df %>%
    pivot_wider(names_from = !!ensym(predicted_domain), values_from = proportion, values_fill = list(proportion = 0))
  
  # Step 2: Prepare abundance matrix
  abundance_matrix <- as.data.frame(abundance_matrix)
  rownames(abundance_matrix) <- abundance_matrix$slices
  abundance_matrix <- abundance_matrix[,-1]
  
  # Step 3: Clustering based on chosen method
  set.seed(seed)
  if(cluster_method == "hclust"){
    distance_matrix <- dist(abundance_matrix, method = dist_method)
    hclust_result <- hclust(distance_matrix, method = "complete")
    clusters <- as.data.frame(cutree(hclust_result, k = cluster_number))
  }else if(cluster_method == "leiden"){
    similarity_matrix <- cor(t(abundance_matrix))
    graph <- graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    optimal_resolution <- find_optimal_resolution(graph, cluster_number)
    
    leiden_result <- leiden.community(graph, resolution = optimal_resolution)
    clusters <- as.data.frame(leiden_result[["membership"]])
  }else if(cluster_method == "kmeans"){
    abundance_matrix_mat <- as.matrix(abundance_matrix)
    kmeans_result <- kmeans(abundance_matrix_mat, centers = cluster_number)
    clusters <- as.data.frame(kmeans_result$cluster  )
  }
    
  # Step 4: Add additional metadata and calculate ARI and NMI
  colnames(clusters) <- "predicted_class"
  clusters$slices_class <- metadata[match(rownames(clusters), metadata$slices) ,match(slices_class,colnames(metadata))]
  
  clusters$slices <- rownames(clusters)
  ari <- aricode::ARI(as.character(clusters$predicted_class), as.character(clusters$slices_class))
  nmi <- aricode::NMI(as.character(clusters$predicted_class), as.character(clusters$slices_class))
  
  # Step 5: Create an abundance matrix with slice classes
  abundance_matrix_re <- as.data.frame(abundance_matrix)
  abundance_matrix_re$slices_class <-  clusters[match(rownames(abundance_matrix_re), rownames(clusters)),]$slices_class
  abundance_matrix_re$slices <- rownames(abundance_matrix_re)
  
  # Return results: ARI, NMI, plot, and abundance matrix
  return(list(ari = ari, nmi = nmi, abundance = abundance_matrix_re))
}

#slices_calculation Function ----
# Function to perform clustering and calculate metrics across multiple datasets
# This function processes each dataset in a given directory, performs clustering, and stores results such as ARI, NMI, and PCA plots.
# It also generates and saves summary plots for ARI and NMI across all datasets.
slices_calculation <- function(original_file, file, model, cluster_number,cluster_method,predicted_domain = "predicted_domain", dist_method = "euclidean", slices_class = "slices_class"){
  fl = c("4","5","6","7","8","9","10")
  metric <- matrix(,,4)
  plot_list <- list()
  
  data_original = read_h5ad(original_file)
  metadata_original = data_original$obs
  
  for(f in fl){
    ex <- paste(paste(file,f,sep = "/"), paste(model,"h5ad",sep = "."), sep = "/")  
    data <- read_h5ad(ex)
    metadata<-data$obs
    
    inter_barcode <- intersect(metadata_original$barcode, metadata$barcode)
    metadata <- metadata[match(inter_barcode, metadata$barcode),]
    metadata$slices_class <- metadata_original[match(inter_barcode, metadata_original$barcode),]$slices_class
    
    slices_clustering_re <- slices_clustering(metadata,predicted_domain, dist_method, cluster_method, cluster_number, slices_class)
    metric <- rbind(metric, c(length(unique(metadata$predicted_domain)), slices_clustering_re$ari, slices_clustering_re$nmi, f)) 
    plot_list[[f]] <- slices_clustering_re$plot
    
    if(cluster_method == "hclust"){
      abundance_matrix <- slices_clustering_re$abundance
      write.table(abundance_matrix, file = paste(file, f, "abundance_matrix_normal.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
    }
  }
  
  metric <- metric[-1,]
  colnames(metric) <- c("domain_count", "ari", "nmi", "parameter_domain")
  metric <- as.data.frame(metric)
  metric$domain_count <- as.numeric(metric$domain_count)
  metric$ari <- as.numeric(metric$ari)
  metric$nmi <- as.numeric(metric$nmi)

  output_dir <- file.path(file, "Metric", cluster_method)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.table(metric, file = paste(output_dir, "model_metric.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
  custom_theme <- theme(
    axis.text.x = element_text(size = 4, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.key.size = unit(0.2, "cm")
  )
  
  plot_ari <- ggplot(metric, aes(x = domain_count)) +
    geom_bar(aes(y = ari), stat = "identity", fill = "#E89C9A", alpha = 0.5) +
    geom_line(aes(y = ari), color = "#E89C9A", size = 0.5) +
    geom_point(aes(y = ari), color = "#E89C9A", size = 1) +
    labs(title = "ARI vs Domain Count", x = "Domain Count", y = "ARI") +
    theme_minimal()+
    custom_theme
  
  plot_nmi <- ggplot(metric, aes(x = domain_count)) +
    geom_bar(aes(y = nmi), stat = "identity", fill = "#ACCBDE", alpha = 0.5) +
    geom_line(aes(y = nmi), color = "#ACCBDE", size = 0.5) +
    geom_point(aes(y = nmi), color = "#ACCBDE", size = 1) +
    labs(title = "NMI vs Domain Count", x = "Domain Count", y = "NMI") +
    theme_minimal()+
    custom_theme
  
  combined_plot <- plot_ari / plot_nmi
  ggsave(paste(output_dir, "metric.pdf", sep = "/"), combined_plot,height = 8,width = 6, units = "cm")
  
  domain_plot <- plot_list[[1]]
  for(d in 2:length(plot_list)){
    domain_plot<- domain_plot + plot_list[[d]]
  }
  ggsave(paste(output_dir, "slices_plot.pdf", sep = "/"), domain_plot,height = 40,width = 50, units = "cm")
  
}


args <- commandArgs(trailingOnly = TRUE)
original_file <- args[1]
file <- args[2]
model <- args[3]
cluster_number <- as.numeric(args[4])
cluster_method <- ifelse(length(args) >= 4, args[5], "leiden")
predicted_domain <- ifelse(length(args) >= 5, args[6], "predicted_domain")
dist_method <- ifelse(length(args) >= 6, args[7], "euclidean")
slices_class <- ifelse(length(args) >= 7, args[8], "slices_class")

slices_calculation(original_file, file, model, cluster_number,cluster_method, predicted_domain, dist_method, slices_class)

