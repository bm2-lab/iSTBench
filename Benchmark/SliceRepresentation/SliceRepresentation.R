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

use_python("/home/dongkj/anaconda3/envs/MultiSpatial/bin/python")

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
  }else{
    similarity_matrix <- cor(t(abundance_matrix))
    graph <- graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    optimal_resolution <- find_optimal_resolution(graph, cluster_number)
    
    leiden_result <- leiden.community(graph, resolution = optimal_resolution)
    clusters <- as.data.frame(leiden_result[["membership"]])
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
  
  # Step 6: Plot PCA and return results
  set.seed(1234)
  pca_result <- prcomp(abundance_matrix, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x)
  pca_df$slices <- rownames(pca_df)
  
  pca_df <- pca_df %>%
    left_join(clusters, by = "slices")
  
  pca_df$slices_class <- as.character(pca_df$slices_class)
  
  # Create plot title including ARI and NMI metrics
  plot_title <- paste("Domain count:", length(unique(metadata$predicted_domain)), "; ARI:", round(ari, 5), "; NMI:", round(nmi, 5), sep = "")
  
  # Plot PCA with clusters colored by slices_class
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = slices_class)) +
    geom_point(size = 3) +
    labs(title = plot_title, x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(plot.title = element_text(size = 6),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6))
  
  # Return results: ARI, NMI, plot, and abundance matrix
  return(list(ari = ari, nmi = nmi, plot = p, abundance = abundance_matrix_re))
}

#slices_calculation Function ----
# Function to perform clustering and calculate metrics across multiple datasets
# This function processes each dataset in a given directory, performs clustering, and stores results such as ARI, NMI, and PCA plots.
# It also generates and saves summary plots for ARI and NMI across all datasets.
slices_calculation <- function(file, model, cluster_number,cluster_method,predicted_domain = "predicted_domain", dist_method = "euclidean", slices_class = "slices_class"){
  
  # Step 1: Read all files from the specified directory
  fl <- list.files(file)
  fl <- fl[-grep("Metric", fl)]  # Exclude files related to metrics
  metric <- matrix(,,4)
  plot_list <- list()
  
  # Step 2: Process each dataset
  for(f in fl){
    ex <- paste(paste(file,f,sep = "/"), paste(model,"h5ad",sep = "."), sep = "/")  
    data <- read_h5ad(ex)
    metadata<-data$obs
    
    # Perform clustering and get results
    slices_clustering_re <- slices_clustering(metadata,predicted_domain, dist_method, cluster_method, cluster_number, slices_class)
    metric <- rbind(metric, c(length(unique(metadata$predicted_domain)), slices_clustering_re$ari, slices_clustering_re$nmi, f)) 
    plot_list[[f]] <- slices_clustering_re$plot
    
    # Save abundance matrix
    abundance_matrix <- slices_clustering_re$abundance
    write.table(abundance_matrix, file = paste(file, f, "abundance_matrix_normal.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  }
  
  # Step 3: Compile and save metrics
  metric <- metric[-1,]
  colnames(metric) <- c("domain_count", "ari", "nmi", "parameter_domain")
  metric <- as.data.frame(metric)
  metric$domain_count <- as.numeric(metric$domain_count)
  metric$ari <- as.numeric(metric$ari)
  metric$nmi <- as.numeric(metric$nmi)
  write.table(metric, file = paste(paste(file, "Metric", cluster_method,sep = "/"), "model_metric.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
  # Step 4: Generate and save combined ARI and NMI plots
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
  
  # Save the combined ARI and NMI plots
  combined_plot <- plot_ari / plot_nmi
  ggsave(paste(paste(file, "Metric",cluster_method, sep = "/"), "metric.pdf", sep = "/"), combined_plot,height = 8,width = 6, units = "cm")
  
  # Step 5: Combine and save all slice plots
  domain_plot <- plot_list[[1]]
  for(d in 2:length(plot_list)){
    domain_plot<- domain_plot + plot_list[[d]]
  }
  ggsave(paste(paste(file, "Metric",cluster_method, sep = "/"), "slices_plot.pdf", sep = "/"), domain_plot,height = 40,width = 50, units = "cm")
}



# Parse command-line arguments using getopt
opt = matrix(c("file", "f", 1, "character", "Input file path",
               "model", "m", 1, "character", "The model used for integrating",
               "cluster_number", "n", 1, "numeric", "Number of sample classes",
               "cluster_method", "M", 2, "character", "Cluster method, default is hclust",
               "predicted_domain", "p", 2, "character", "The name of metadata that stores domains",
               "dist_method", "d", 2, "character", "The method used for calculating distance between slices, default is euclidean",
               "slices_class", "s", 2, "character", "The name of metadata that stores slice classes"),
             byrow = TRUE, ncol = 5)

args = getopt(opt)

# Set default values for unspecified arguments
if (is.null(args$cluster_method)) {
  args$cluster_method = "hclust"
}

if (is.null(args$predicted_domain)) {
  args$predicted_domain = "predicted_domain"
}

if (is.null(args$dist_method)) {
  args$dist_method = "euclidean"
}

if (is.null(args$slices_class)) {
  args$slices_class = "slices_class"
}

# Run the slices_calculation function with parsed arguments
slices_calculation(args$file, args$model, args$cluster_number, args$cluster_method, args$predicted_domain, args$dist_method, args$slices_class)


