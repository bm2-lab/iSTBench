library(circlize)
library(reticulate)
library(Seurat)
library(SeuratDisk)
library(anndata)  
library(ggplot2)
library(MASS)
library(dplyr)
library(tidyr)
library(stats)
library(umap)
library(ggplot2)
library(patchwork)
library(aricode)
library(leidenAlg)
library(igraph)
library(getopt)

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

slices_clustering <- function(domain_abundance, domain_cor, cluster_number, cluster_method = "leiden", dist_method = "euclidean", seed = 123){
  set.seed(seed)
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

metric_plot <- function(predicted_class){
  
  predicted_class$slices <- rownames(predicted_class)
  ari <- aricode::ARI(as.character(predicted_class$predicted_class), as.character(predicted_class$slices_class))
  nmi <- aricode::NMI(as.character(predicted_class$predicted_class), as.character(predicted_class$slices_class))
  return(list(ari = ari, nmi = nmi))
  
}

slices_calculation <- function(file, cluster_number,abundance_data, predicted_domain, slices_class, sd_filter = 0){
  fl <- list.files(file)
  fl <- fl[-grep("Metric", fl)]
  metric_hclust <- matrix(,,4)
  metric_leiden <- matrix(,,4)
  
  for(f in fl){
    ###1.read data
    f1 <- paste(file, f, abundance_data, sep = "/") 
    f2 <- paste(file, f, "Slices_DomainCor_embedding.csv", sep = "/") 
    
    domain_abundance <- read.csv(f1)
    domain_count <- (ncol(domain_abundance)-2)
    domain_abundance_2 <- domain_abundance[,-(match(c("slices", "slices_class"), colnames(domain_abundance)))]
    domain_abundance_2 <- t(apply(domain_abundance_2, 1, function(x){
      2*x-1
    }))
    domain_abundance[,-(match(c("slices", "slices_class"), colnames(domain_abundance)))] <- domain_abundance_2
    rownames(domain_abundance) <- domain_abundance$slices
    
    domain_cor <- read.csv(f2)
    domain_cor <- domain_cor[,-which(colnames(domain_cor) == "slices_class")]
    rownames(domain_cor) <- domain_cor$slices
    domain_cor <- domain_cor[rownames(domain_abundance),]
    
    slices_class <- domain_abundance$slices_class
    domain_abundance <- domain_abundance[,-match(c("slices", "slices_class"), colnames(domain_abundance))]
    domain_cor <- domain_cor[,-match(c("slices"), colnames(domain_cor))]
    
    ###2. data filtering
    cor_sd <- apply(domain_cor, 2, sd)
    domain_cor <- domain_cor[, which(cor_sd > sd_filter)]
    
    ###3. clustring
    predicted_class_hclust <- slices_clustering(domain_abundance, domain_cor, cluster_number, cluster_method = "hclust",dist_method = "euclidean")
    distance_matrix_hclust <- as.matrix(predicted_class_hclust$distance_matrix)
    predicted_class_hclust <- predicted_class_hclust$clusters
    predicted_class_hclust$slices_class <- slices_class
    
    predicted_class_leiden <- slices_clustering(domain_abundance, domain_cor, cluster_number, cluster_method = "leiden",dist_method = "euclidean")
    distance_matrix_leiden <- predicted_class_leiden$distance_matrix
    predicted_class_leiden <- predicted_class_leiden$clusters
    predicted_class_leiden$slices_class <- slices_class
    
    slices_clustering_re_hclust <- metric_plot(predicted_class_hclust)
    slices_clustering_re_leiden <- metric_plot(predicted_class_leiden)
    
    domain_embedding <- cbind(domain_abundance, domain_cor)
    domain_embedding$slices <- rownames(domain_abundance)
    domain_embedding$slices_class <- slices_class
    colnames(domain_embedding) <- c(rep("domain_abundance", ncol(domain_abundance)), rep("domain_cor", ncol(domain_cor)), "slices", "slices_class")
    
    ab <- unlist(strsplit(unlist(strsplit(abundance_data, split = "[.]"))[1], split = "_"))[3]
    write.table(domain_embedding, file = paste(file, f, paste(paste("DomainAbunCor", ab, sd_filter, sep = "_"), "csv", sep = "."), sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
    
    colnames(distance_matrix_leiden) <- paste("S", colnames(distance_matrix_leiden), sep = "_")
    colnames(distance_matrix_hclust) <- paste("S", colnames(distance_matrix_hclust), sep = "_")
    write.table(distance_matrix_leiden, file = paste(file, f, paste(paste("DistanceMatrixLeiden", ab, sd_filter, sep = "_"), "csv", sep = "."), sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
    write.table(distance_matrix_hclust, file = paste(file, f, paste(paste("DistanceMatrixHclust", ab, sd_filter, sep = "_"), "csv", sep = "."), sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
    
    
    ####----
    metric_hclust <- rbind(metric_hclust, c(domain_count, slices_clustering_re_hclust$ari, slices_clustering_re_hclust$nmi, f)) 
    metric_leiden <- rbind(metric_leiden, c(domain_count, slices_clustering_re_leiden$ari, slices_clustering_re_leiden$nmi, f)) 
  }
  
  ###4.plot hclust
  metric <- metric_hclust
  metric <- metric[-1,]
  colnames(metric) <- c("domain_count", "ari", "nmi", "parameter_domain")
  metric <- as.data.frame(metric)
  
  ab <- unlist(strsplit(unlist(strsplit(abundance_data, split = "[.]"))[1], split = "_"))[3]
  output_file = paste(file, "Metric",paste("DomainAbunCor","hclust", ab, sd_filter, sep = "_"),sep = "/")
  if (!dir.exists(output_file)) {dir.create(output_file)}
  write.table(metric, file = paste(output_file , "model_metric.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
  metric$ari <- as.numeric(metric$ari)
  metric$nmi <- as.numeric(metric$nmi)
  metric$domain_count <- as.numeric(metric$domain_count)
  
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
  ggsave(paste(output_file, "metric.pdf", sep = "/"), combined_plot,height = 8,width = 6, units = "cm")
  
  ###plot leiden
  metric <- metric_leiden
  metric <- metric[-1,]
  colnames(metric) <- c("domain_count", "ari", "nmi", "parameter_domain")
  metric <- as.data.frame(metric)
  
  metric$ari <- as.numeric(metric$ari)
  metric$nmi <- as.numeric(metric$nmi)
  metric$domain_count <- as.numeric(metric$domain_count)
  
  ab <- unlist(strsplit(unlist(strsplit(abundance_data, split = "[.]"))[1], split = "_"))[3]
  output_file = paste(file, "Metric",paste("DomainAbunCor","leiden", ab, sd_filter, sep = "_"),sep = "/")
  if (!dir.exists(output_file)) {dir.create(output_file)}
  write.table(metric, file = paste(output_file , "model_metric.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
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
  ggsave(paste(output_file, "metric.pdf", sep = "/"), combined_plot,height = 8,width = 6, units = "cm")
  
}


# Get the corresponding value from the command line argument
opt = matrix(c("file", "f", 1, "character", "Input file path",
               "cluster_number", "cn", 1, "numeric", "Number of sample classes",
               "abundance_data", "ad", 1, "character", "Slice representation based domain abundance",
               "predicted_domain", "pd", 2, "character", "The name of metadata that stores domains",
               "slices_class", "sc", 2, "character", "The name of metadata that stores slice classes",
               "sd_filter", "sd", 2, "numeric", "Thresholds for filtering domains"),
             byrow = TRUE, ncol = 5)

args = getopt(opt)

if (is.null(args$sd_filter)) {
  args$sd_filter = 0
}

if (is.null(args$predicted_domain)) {
  args$predicted_domain = "predicted_domain"
}

if (is.null(args$slices_class)) {
  args$slices_class = "slices_class"
}

file <- args$file
cluster_number <- args$cluster_number
abundance_data <- args$abundance_data
predicted_domain <- as.numeric(args$predicted_domain)
slices_class <- args$slices_class
sd_filter <- args$sd_filter

# Call the function
slices_calculation(file, cluster_number,abundance_data, predicted_domain, slices_class, sd_filter)

