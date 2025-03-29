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
library(MASS)
# library(rlang)
# use_python("/home/dongkj/anaconda3/envs/MultiSpatial/bin/python")

#Domain abundance Function ----
# This function calculates the abundance matrix for each slice.
# It returns a matrix of domain proportions for each slice along with class information.
domainAbundance <- function(metadata, predicted_domain = "predicted_domain", slices_class = "slices_class"){

  # Step 1: Calculate abundance matrix for each slice
  # Summarizes the counts of predicted domains per slice and calculates the proportion for each domain
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
    dplyr::select(slices, !!ensym(predicted_domain), proportion)

  # Step 2: Prepare abundance matrix 
  # Converts the abundance data into a wide format and assigns slices as row names
  abundance_matrix <- abundance_df %>%
    pivot_wider(names_from = !!ensym(predicted_domain), values_from = proportion, values_fill = list(proportion = 0))

  # Step 3: Add slice class information 
  # Adds a slice class column from the metadata based on slice names
  abundance_matrix <- as.data.frame(abundance_matrix)
  rownames(abundance_matrix) <- abundance_matrix$slices
  abundance_matrix <- abundance_matrix[,-1]

  # Create an abundance matrix with slice classes
  abundance_matrix$slices_class <-  metadata[match(rownames(abundance_matrix), metadata$slices),which(colnames(metadata) == slices_class)]
  abundance_matrix$slices <- rownames(abundance_matrix)
  return(abundance_matrix)
}

#Domain cor Function ----
# Computes the correlation between domains for each slice, using 2D KDE density estimation
domainCor <- function(data, predicted_domain, step){
  data_metaa_all <- data$obs
  data_coor_all <- data$obsm$spatial
  slices <- unique(data_metaa_all$slices)
  
  # Loop through each slice and calculate domain correlations
  # For each slice, domain-specific 2D density distributions are computed
  domain_cor_all <- list()
  for(s in slices){
    s_index <- which(data_metaa_all$slices == s)
    data_meta <- data_metaa_all[s_index, ]
    data_coor <- data_coor_all[s_index, ]
    domains <- unique(data_meta[[predicted_domain]])
    domains <-sort(domains)
    min_x <- min(data_meta$X)
    min_y <- min(data_meta$Y)
    max_x <- max(data_meta$X)
    max_y <- max(data_meta$Y)
    
    # For each domain, estimate its 2D spatial density distribution
    # If there is only one point, assume zero density; otherwise, compute KDE
    domains_density<-data.frame()
    for(d in domains){
      domain_data <- data_meta[data_meta[[predicted_domain]] == d,]
      if(nrow(domain_data) == 1){
        domain_density <-  rep(0, step*step)
        domains_density <- rbind(domains_density,domain_density)
      }else{
        if(length(unique(domain_data$X)) == 1){
          domain_data$X[1] <- domain_data$X[1]+1
        }
        
        if(length(unique(domain_data$Y)) == 1){
          domain_data$Y[1] <- domain_data$Y[1]+1
        }
        
        domain_density <- kde2d(domain_data$X, domain_data$Y, n = step, lims=c(min_x, max_x, min_y, max_y))
        domains_density <- rbind(domains_density, as.vector(domain_density$z))
      }
    }
    
    rownames(domains_density) <- domains
    cor_matrix <- cor(t(domains_density), method = "pearson")
    cor_matrix[which(is.na(cor_matrix))] <- (-1)
    domain_cor_all[[s]] <- cor_matrix
  }
  
  # Harmonize domain correlations
  # Merge correlation matrices across all slices and ensure consistent row/column names
  common_row_names <- Reduce(union, lapply(domain_cor_all, rownames))
  common_col_names <- Reduce(union, lapply(domain_cor_all, colnames))
  
  filtered_domain_cor_all <- lapply(domain_cor_all, function(mat) {
    new <- matrix(-1,length(common_row_names),length(common_row_names))
    colnames(new) <- common_col_names
    rownames(new) <- common_row_names
    diag(new) <- 1
    
    inter_row <- intersect(rownames(mat), common_row_names)
    inter_col <- intersect(colnames(mat), common_col_names)
    for(r in inter_row){
      for(c in inter_col){
        new[r,c] <- mat[r,c]
      }
    }
    return(new)
  })
  
  # Convert correlation matrices to vectors
  # Flatten the upper triangular part of each matrix for further analysis
  vectorized_matrices <- lapply(filtered_domain_cor_all, function(mat) {
    return(mat[upper.tri(mat)])
  })
  
  cor_matrix <- do.call(rbind, vectorized_matrices)
  cor_matrix <- as.data.frame(cor_matrix)
  cor_matrix$slices <- slices

  cor_matrix <- as.data.frame(cor_matrix)
  rownames(cor_matrix) <- as.character(cor_matrix$slices)
  cor_matrix <- cor_matrix[,-match("slices", colnames(cor_matrix))]
  return(cor_matrix)
}

# Slice clustering function ----
# Performs clustering on slices based on domain abundance and correlation data
slices_clustering <- function(domain_abundance, domain_cor, cluster_number, cluster_method = "leiden", dist_method = "euclidean", seed = 123){
  set.seed(seed)
  
  # If hierarchical clustering method is selected
  # Calculate distance matrices for domain abundance and correlation, then combine them
  if(cluster_method == "hclust"){
    distance_matrix1 <- dist(domain_abundance, method = dist_method)
    distance_matrix1 <- as.matrix(distance_matrix1)
    distance_matrix1 <- distance_matrix1/ncol(domain_abundance)
    
    distance_matrix2 <- dist(domain_cor, method = dist_method)
    distance_matrix2 <- as.matrix(distance_matrix2)
    distance_matrix2 <- distance_matrix2/ncol(domain_cor)
    
    mean_matrix <- (distance_matrix1 + distance_matrix2) / 2
    final_distance_matrix <- as.dist(mean_matrix)
    
    hclust_result <- hclust(final_distance_matrix, method = "complete")
    clusters <- as.data.frame(cutree(hclust_result, k = cluster_number))
  }else{
    # If leiden clustering method is selected
    # Calculate similarity matrices from correlations and apply leiden algorithm
    similarity_matrix1 <- cor(t(domain_abundance))
    similarity_matrix2 <- cor(t(domain_cor))
    similarity_matrix1[which(is.na(similarity_matrix1))] <- -1
    similarity_matrix2[which(is.na(similarity_matrix2))] <- -1
    similarity_matrix <- (similarity_matrix1 + similarity_matrix2) / 2
    
    final_distance_matrix <- (1 - similarity_matrix)
    
    graph <- graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    optimal_resolution <- find_optimal_resolution(graph, cluster_number)
    
    leiden_result <- leiden.community(graph, resolution = optimal_resolution)
    clusters <- as.data.frame(leiden_result[["membership"]])
  }
  
  rownames(clusters) <- rownames(domain_abundance)
  colnames(clusters) <- "predicted_class"
  return(list(clusters = clusters, distance_matrix = final_distance_matrix))
}

# Optimal resolution finder for leiden clustering ----
# This function finds the best resolution parameter for leiden clustering based on the target number of clusters
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

metric_plot <- function(predicted_class){
  
  predicted_class$slices <- rownames(predicted_class)
  ari <- aricode::ARI(as.character(predicted_class$predicted_class), as.character(predicted_class$slices_class))
  nmi <- aricode::NMI(as.character(predicted_class$predicted_class), as.character(predicted_class$slices_class))
  return(list(ari = ari, nmi = nmi))
  
}

# slices_calculation function ----
# This function performs several steps for domain abundance and correlation analysis, clustering, and metric evaluation.
# It reads in the data, processes domain abundance, calculates correlations, performs clustering, and outputs the results.
# The function also handles domain filtering, embedding generation, and visualization of clustering metrics.
slices_calculation <- function(file, output, cluster_number, cluster_method, predicted_domain, slices_class, dist_method, step, sd_filter){
  # Step 1: Data loading and preprocessing
  # Read the input data from an H5AD file, extract metadata, and process domain abundance.
  data <- read_h5ad(file)
  metadata<-data$obs
  metadata[[predicted_domain]] <- as.character(metadata[[predicted_domain]])
  
  # Calculate domain abundance for each slice
  domain_abundance <- domainAbundance(metadata, predicted_domain, slices_class)
  
  # Step 2: Transform domain abundance values
  # Scale domain abundance values from [0,1] range to [-1,1].
  domain_count <- (ncol(domain_abundance)-2)
  domain_abundance_2 <- domain_abundance[,-(match(c("slices", slices_class), colnames(domain_abundance)))]
  domain_abundance_2 <- t(apply(domain_abundance_2, 1, function(x){
    2*x-1
  }))
  # Replace transformed values back into the domain_abundance dataframe
  domain_abundance[,-(match(c("slices", "slices_class"), colnames(domain_abundance)))] <- domain_abundance_2
  rownames(domain_abundance) <- domain_abundance$slices
  
  # Step 3: Calculate domain correlations
  # Compute pairwise domain correlations for the given step size.
  domain_cor <- domainCor(data, predicted_domain, step)
  domain_cor <- domain_cor[rownames(domain_abundance),]
  # Extract the slices_class column for later use
  slices_class <- domain_abundance$slices_class
  domain_abundance <- domain_abundance[,-match(c("slices", "slices_class"), colnames(domain_abundance))]
  
  # Step 4: Domain correlation filtering
  # Filter out columns from the domain correlation matrix with low standard deviation.
  cor_sd <- apply(domain_cor, 2, sd)
  domain_cor <- domain_cor[, which(cor_sd > sd_filter)]
  
  # Step 5: Clustering
  # Perform clustering using either hierarchical clustering or Leiden algorithm based on user input.
  predicted_class_leiden <- slices_clustering(domain_abundance, domain_cor, cluster_number, cluster_method = cluster_method,dist_method = dist_method)
  # Extract clustering results and distance matrix
  distance_matrix_leiden <- predicted_class_leiden$distance_matrix
  predicted_class_leiden <- predicted_class_leiden$clusters
  predicted_class_leiden$slices_class <- slices_class
  
  # Step 6: Plot clustering metrics 
  # Generate and plot clustering metrics (Adjusted Rand Index and Normalized Mutual Information).
  slices_clustering_re_leiden <- metric_plot(predicted_class_leiden)
  
  # Step 7: Domain embedding generation
  # Combine domain abundance and domain correlation into a single embedding matrix.
  domain_embedding <- cbind(domain_abundance, domain_cor)
  domain_embedding$slices <- rownames(domain_abundance)
  domain_embedding$slices_class <- slices_class
  colnames(domain_embedding) <- c(rep("domain_abundance", ncol(domain_abundance)), rep("domain_cor", ncol(domain_cor)), "slices", "slices_class")
  
  # write.table(domain_embedding, file = paste(output, paste(paste("DomainAbunCor", sd_filter, sep = "_"), "csv", sep = "."), sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
  colnames(distance_matrix_leiden) <- paste("S", colnames(distance_matrix_leiden), sep = "_")
  # write.table(distance_matrix_leiden, file = paste(output, paste(paste("DistanceMatrixLeiden", sd_filter, sep = "_"), "csv", sep = "."), sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
  write.table(predicted_class_leiden, file = paste(output, paste(paste("PredictedClass", sd_filter, sep = "_"), "csv", sep = "."), sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
  metric <- data_frame(metric = c("ari", "nmi"), value = c(slices_clustering_re_leiden$ari, slices_clustering_re_leiden$nmi))
  metric <- as.data.frame(metric)
  
  # Step 9: Plot final metric evaluation 
  # Create a bar plot for the clustering metrics and save as a PDF.
  plot_metric <- ggplot(metric, aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity") +
    labs(title = "Metric Values",
         x = "Metric",
         y = "Value") +
    geom_text(aes(label = round(value, 2)), vjust = -0.5, color = "black", size = 3)+
    theme_minimal()
  ggsave(paste(output, "metric.pdf", sep = "/"), plot_metric,height = 10,width = 8, units = "cm")
}


# Command line argument parsing ----
# This section handles the command-line arguments for the script, setting default values if necessary.
opt = matrix(c("file", "f", 1, "character", "Input file path",
               "output", "o", 1, "character", "Output file",
               "cluster_number", "c", 1, "numeric", "Number of sample classes",
               "cluster_method", "m", 2, "character", "hclust or leiden, default is leiden",
               "predicted_domain", "p", 2, "character", "The name of metadata that stores domains",
               "slices_class", "s", 2, "character", "The name of metadata that stores slice classes. If the slices class is not known, you need to add this column of information to the metadata of the data, which can be set to unknown.",
               "dist_method", "d", 2, "character", "If hclust clustering is carried out, the method of calculating distance is adopted. Default is euclidean",
               "sd_filter", "S", 2, "numeric", "Thresholds for filtering domains correlation embedding",
               "step", "P", 2, "numeric", "The step size is selected when the domain density distribution is collected. Default is 20"),
             byrow = TRUE, ncol = 5)

args = getopt(opt)

if (is.null(args$cluster_method)) {
  args$cluster_method = "leiden"
}

if (is.null(args$predicted_domain)) {
  args$predicted_domain = "predicted_domain"
}

if (is.null(args$slices_class)) {
  args$slices_class = "slices_class"
}

if (is.null(args$dist_method)) {
  args$dist_method = "euclidean"
}

if (is.null(args$sd_filter)) {
  args$sd_filter = 0
}

if (is.null(args$step)) {
  args$step = 20
}

file <- args$file
output <- args$output
cluster_number <- args$cluster_number
cluster_method <- args$cluster_method
predicted_domain <- args$predicted_domain
slices_class <- args$slices_class
dist_method <- args$dist_method
sd_filter <- args$sd_filter
step <- args$step

# Call the function
slices_calculation(file, output, cluster_number, cluster_method, predicted_domain, slices_class, dist_method, step, sd_filter)

