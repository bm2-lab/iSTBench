generate_spatial_plots <- function(data_file, output_file) {
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
  
  
  # 读取数据文件
  spatial_data_file <- read.table(data_file, header = FALSE)
  
  
  #DLPFC data color----
  #color <- c("L1" = "#ACCBDE", "L2" = "#3C76AE", "L3" = "#BAD691", "L4" = "#559D52", "L5" = "#E89C9A", "L6" = "#CF3732", "WM" ="#F2C07A" )
  color2 <- c("#E89C9A","#ACCBDE","#CF3732","#BAD691","#002FA7","#008C8C","#81D8CF","#B05923","#900021","#E60000","#FBD26A","#E85827","#432913")
  
  
  
  # 处理数据
  metadata <- matrix(,dim(read_h5ad(unlist(strsplit(spatial_data_file[1, 1], split = ":"))[2]))[1],)
  for(i in 1:nrow(spatial_data_file)){
    model <- unlist(strsplit(spatial_data_file[i, 1], split = ":"))[1]
    ex <- unlist(strsplit(spatial_data_file[i, 1], split = ":"))[2]
    data <- read_h5ad(ex)
    
    if(i == 1){
      data_meta<-data$obs
      metadata <- cbind(metadata, data_meta[,match(c("barcode","X","Y","slices","original_domain"), colnames(data_meta))]) 
      predicted_domain <- as.data.frame(as.character(data_meta$predicted_domain))
      colnames(predicted_domain) <- model
      metadata <- cbind(metadata,predicted_domain)
    } else {
      data_meta<-data$obs
      index <- match(metadata[,2], data_meta$barcode)
      predicted_domain <- as.data.frame(as.character(data_meta[index,]$predicted_domain))
      colnames(predicted_domain) <- model
      metadata <- cbind(metadata,predicted_domain)
      
    }
  }
  
  metadata <- metadata[, -1]
  metadata <- as.data.frame(metadata)
  
  ex <- unlist(strsplit(spatial_data_file[1, 1], split = ":"))[2]
  #paste(paste(unlist(strsplit(ex, split = "/"))[1:7], collapse = "/"), "sample_all_data/Slices_combind_data.h5ad", sep = "/")
  
  if(grepl("seed", ex)){
    Intergration_Benchmark
    spatial <- read_h5ad(paste(paste(unlist(strsplit(gsub("Intergration_Benchmark_seed=[0-9]+", "Intergration_Benchmark", ex), split = "/"))[1:8], collapse = "/"), "sample_all_data/Slices_combind_data.h5ad", sep = "/"))
  }else{
    spatial <- read_h5ad(paste(paste(unlist(strsplit(ex, split = "/"))[1:8], collapse = "/"), "sample_all_data/Slices_combind_data.h5ad", sep = "/"))
  }
  
  index <- match(metadata$barcode, spatial$obs$barcode)
  spatial <- spatial$obsm$spatial[index,]
  metadata$X <- as.numeric(spatial[,1])
  metadata$Y <- as.numeric(spatial[,2])
  #metadata$original_domain <- factor(metadata$original_domain)
  metadata<-na.omit(metadata)
  write.table(metadata,paste(output_file, "metaPredictedRe.csv", sep = "/"), col.names = T,row.names = F,sep = ",",quote = F)
  
  
  model_names <- colnames(metadata)[-c(1:5)]
  model_ari <- c()
  model_nmi <- c()
  for(m in model_names){
    o_domain <- metadata$original_domain
    p_domain <- metadata[,match(m,colnames(metadata))]
    
    ari <- ARI(as.character(o_domain), as.character(p_domain))
    nmi <- NMI(as.character(o_domain), as.character(p_domain))
    model_ari <- c(model_ari, ari)
    model_nmi <- c(model_nmi, nmi)
  }
  
  model_metric <- cbind(model_names,cbind(model_ari,model_nmi))
  model_metric <- as.data.frame(model_metric)
  model_metric$model_ari <- as.numeric(model_metric$model_ari)
  model_metric$model_nmi <- as.numeric(model_metric$model_nmi)
  write.table(model_metric,paste(output_file, "model_metric.csv", sep = "/"), col.names = T,row.names = F,sep = ",",quote = F)
  
  ari_plot <- ggplot(model_metric, aes(x = model_names, y = model_ari, fill = model_names)) +
    geom_bar(stat = "identity") +
    labs(title = "ARI Values", x = "Model", y = "ARI") +
    theme_minimal()+
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
    theme_minimal()+
    theme(plot.title = element_text(size = 3),
          legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 3),
          legend.title = element_text(size = 3),
          axis.title = element_text(size = 3),
          axis.text = element_text(size = 3),
          axis.text.x = element_text(angle = 30, hjust = 1)) 
  
  metric_plot <- ari_plot + nmi_plot
  ggsave(paste(output_file, "metric.pdf", sep = "/"), metric_plot, width = 10, height = 5, units = "cm",limitsize = FALSE)
  
  sample_names <- unique(metadata$slices)
  spatial_plots <- lapply(sample_names, function(snm) {
    sample_data <- metadata[metadata$slices == snm,]
    
    plot_list <- list()
    plot_list[[1]] <- ggplot(sample_data, aes(x = X, y = Y, col = original_domain)) +
      geom_point(size = 0.1) + 
      #scale_color_manual(values = color, name = "") +
      theme_classic() + 
      theme(
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
      labs(title = "Ground truth", y = paste("slices", as.character(snm))) +
      coord_equal()+
      theme(plot.title = element_text(size = 6),
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.1, 'cm'),
            axis.title.y = element_text(size = 6),
            legend.title = element_text(size = 6)) 
    
    for(z in 1:length(model_names)){
      m <- model_names[z]
      ari <- model_metric[z, 2]
      
      df <- cbind.data.frame(x_index = sample_data$X, y_index = sample_data$Y, predicted_domain = sample_data[, model_names[z]])
      df$predicted_domain <- as.character(df$predicted_domain)
      
      plot_list[[z+1]] <- ggplot(df, aes(x = x_index, y = y_index, col = predicted_domain )) +
        geom_point(size = 0.1) + 
        #scale_color_manual(values = color2[1: length(unique(df$predicted_domain))], name = "") +
        theme_classic() + 
        theme(
          legend.position = "right",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
        labs(title =  sprintf("%s - ARI: %s", m, round(ari, 3)) ) +
        coord_equal()+
        theme(plot.title = element_text(size = 6),
              legend.key.size = unit(0.1, 'cm'),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 6)) 
    }
    p <- wrap_plots(plot_list, ncol = length(plot_list))
  })
  
  p <- plot_grid(plotlist = spatial_plots, ncol = 1, byrow = TRUE)
  ggsave(paste(output_file, "domain.pdf", sep = "/"), p, width = 40, height = 25, units = "cm",limitsize = FALSE)
  
  
}

# 调用函数
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
output_file <- args[2]

generate_spatial_plots(data_file, output_file)
