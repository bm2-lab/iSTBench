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

#use_python("/NFS_home/NFS_home_6/dongkj/.local/share/r-miniconda/envs/r-reticulate/bin/python")


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


slices_clustering <- function(metadata, predicted_domain, dist_method, cluster_method, cluster_number, slices_class, seed = 123){
  ###get slices-domain abundance matrix----
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
  
  ###sample_level
  # if(any(grepl("patients", colnames(metadata)))){
  #   slices2patients <- unique(metadata[,match(c("slices","patients"), colnames(metadata))])
  # 
  #   merged_data <- merge(abundance_matrix, slices2patients, by = "slices")
  #   merged_data <- merged_data[,-match("slices", colnames(merged_data))]
  # 
  #   abundance_matrix_patients <- merged_data %>%
  #     group_by(patients) %>%
  #     summarise(across(everything(), mean))
  # 
  #   abundance_matrix_patients <- as.data.frame(abundance_matrix_patients)
  #   rownames(abundance_matrix_patients) <- abundance_matrix_patients$patients
  #   abundance_matrix_patients <- abundance_matrix_patients[,-match("patients", colnames(abundance_matrix_patients)) ]
  #   abundance_matrix <- abundance_matrix_patients
  # }else{
  #   abundance_matrix <- as.data.frame(abundance_matrix)
  #   rownames(abundance_matrix) <- abundance_matrix$slices
  #   abundance_matrix <- abundance_matrix[,-1]
  # }
  
  ###slices level
  abundance_matrix <- as.data.frame(abundance_matrix)
  rownames(abundance_matrix) <- abundance_matrix$slices
  abundance_matrix <- abundance_matrix[,-1]
  
  ###clustering----
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
  
  
  ###calculation----
  ###sample level
  # colnames(clusters) <- "predicted_class"
  # if(any(grepl("patients", colnames(metadata)))){
  #   clusters$slices_class <- metadata[match(rownames(clusters), metadata$patients) ,match(slices_class,colnames(metadata))]
  # }else{
  #   clusters$slices_class <- metadata[match(rownames(clusters), metadata$slices) ,match(slices_class,colnames(metadata))]
  # }
  
  ###slices level
  colnames(clusters) <- "predicted_class"
  clusters$slices_class <- metadata[match(rownames(clusters), metadata$slices) ,match(slices_class,colnames(metadata))]
  
  clusters$slices <- rownames(clusters)
  ari <- aricode::ARI(as.character(clusters$predicted_class), as.character(clusters$slices_class))
  nmi <- aricode::NMI(as.character(clusters$predicted_class), as.character(clusters$slices_class))
  
  abundance_matrix_re <- as.data.frame(abundance_matrix)
  abundance_matrix_re$slices_class <-  clusters[match(rownames(abundance_matrix_re), rownames(clusters)),]$slices_class
  abundance_matrix_re$slices <- rownames(abundance_matrix_re)
  
  ###ploting----
  set.seed(1234)
  
  pca_result <- prcomp(abundance_matrix, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x)
  pca_df$slices <- rownames(pca_df)
  
  pca_df <- pca_df %>%
    left_join(clusters, by = "slices")
  
  pca_df$slices_class <- as.character(pca_df$slices_class)
  
  # 设置图表标题
  plot_title <- paste("Domain count:", length(unique(metadata$predicted_domain)), "; ARI:", round(ari, 5), "; NMI:", round(nmi, 5), sep = "")
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = slices_class)) +
    geom_point(size = 3) +
    labs(title = plot_title, x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(plot.title = element_text(size = 6),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6))
  
  return(list(ari = ari, nmi = nmi, plot = p, abundance = abundance_matrix_re))
}

slices_calculation <- function(file, model, cluster_number,cluster_method,predicted_domain = "predicted_domain", dist_method = "euclidean", slices_class = "slices_class"){
  fl <- list.files(file)
  fl <- fl[-grep("Metric", fl)]
  #fl <- fl[-which(fl == "2")]
  metric <- matrix(,,4)
  plot_list <- list()
  
  for(f in fl){
    ex <- paste(paste(file,f,sep = "/"), paste(model,"h5ad",sep = "."), sep = "/")  
    data <- read_h5ad(ex)
    metadata<-data$obs
    
    slices_clustering_re <- slices_clustering(metadata,predicted_domain, dist_method, cluster_method, cluster_number, slices_class)
    metric <- rbind(metric, c(length(unique(metadata$predicted_domain)), slices_clustering_re$ari, slices_clustering_re$nmi, f)) 
    plot_list[[f]] <- slices_clustering_re$plot
    
    abundance_matrix <- slices_clustering_re$abundance
    write.table(abundance_matrix, file = paste(file, f, "abundance_matrix_normal.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  }
  
  metric <- metric[-1,]
  colnames(metric) <- c("domain_count", "ari", "nmi", "parameter_domain")
  metric <- as.data.frame(metric)
  metric$domain_count <- as.numeric(metric$domain_count)
  metric$ari <- as.numeric(metric$ari)
  metric$nmi <- as.numeric(metric$nmi)
  write.table(metric, file = paste(paste(file, "Metric", cluster_method,sep = "/"), "model_metric.csv", sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
  
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
  ggsave(paste(paste(file, "Metric",cluster_method, sep = "/"), "metric.pdf", sep = "/"), combined_plot,height = 8,width = 6, units = "cm")
  
  domain_plot <- plot_list[[1]]
  for(d in 2:length(plot_list)){
    domain_plot<- domain_plot + plot_list[[d]]
  }
  ggsave(paste(paste(file, "Metric",cluster_method, sep = "/"), "slices_plot.pdf", sep = "/"), domain_plot,height = 40,width = 50, units = "cm")
  
}

###
opt = matrix(c("file", "f", 1, "character", "Input file path",
               "model", "m", 1, "character", "The model used for integrating",
               "cluster_number", "cn", 1, "numeric", "Number of sample classes",
               "cluster_method", "cm", 2, "character", "Cluster method, default is hclust",
               "predicted_domain", "pd", 2, "character", "The name of metadata that stores domains",
               "dist_method", "d", 2, "character", "The method used for calculating distance between slices, default is euclidean",
               "slices_class", "sc", 2, "character", "The name of metadata that stores slice classes"),
             byrow = TRUE, ncol = 5)

args = getopt(opt)

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

slices_calculation(file, model, cluster_number,cluster_method, predicted_domain, dist_method, slices_class)



