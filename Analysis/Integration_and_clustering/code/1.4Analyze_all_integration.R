library(tidyverse)
library(tibble)
library(pheatmap)
library(anndata)
library(patchwork)
library(ggplot2)
library(Seurat)
library(tidyr)
library(funkyheatmap)
library(kableExtra)


###1.read data----
### Intergration mertic value
# The consolidated results from all datasets are read, and the results are merged into a final file for visualization. 
# The related result files, integration_re.csv and integration_re_all.csv, can be directly used in subsequent analyses 
# without needing to repeat this step.

datasets <- c("BaristaSeq","DLPFC_sample1", "DLPFC_sample2", "DLPFC_sample3", "MERFISH", "MERFISH_Brain_S2","MERFISH_Brain_S3",
              "MERFISH_Brain_S4", "MERFISH_Brain_S5", "MERFISH_Brain_S6", "MERFISH_Brain_S7", "MERFISH_Brain_S8",
              "MERFISH_Brain_S9", "MERFISH_Brain_S10", "MERFISH_Brain_S11", "MERFISH_Brain_S12", "STARMap", "Mouse")

models <- c("Banksy", "CellCharter", "CN", "GraphST", "GraphSTwithPASTE", "MENDER", "NicheCompass", "Spado")

intergration_re <- data.frame()
for(d in datasets){
  if(d != "Mouse"){
    data_re <- matrix(,8,11)
    path1 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/model_metric.csv", sep = "/")
    path2 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/CHAOS.csv", sep = "/")
    path3 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/PAS.csv", sep = "/")
    path4 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/Intergration_value.csv", sep = "/")
    
    d_data <- read_h5ad(paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "sample_all_data/Slices_combind_data.h5ad", sep = "/"))
    slices_count <- length(unique(d_data$obs$slices))
    
    metric1 <- read.csv(path1)
    chaos <- read.csv(path2)
    chaos_mean <- apply(chaos[,-which(colnames(chaos) == "Model")],1,mean)
    
    pas <- read.csv(path3)
    pas_mean <- apply(pas[,-which(colnames(pas) == "Model")],1,mean)
    
    intergration_value <- read.csv(path4)
    
    stat_fl <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe", sep = "/")
    stat_flies <- list.files(stat_fl, pattern = "stats")
    time <- data.frame()
    memory <- data.frame()
    for(sf in stat_flies){
      sf_stats <- read.csv(paste(stat_fl, sf, sep = "/"))
      time <- rbind(time, c(unlist(strsplit(sf, split = "_"))[1], sf_stats[1,2]))
      memory <- rbind(memory, c(unlist(strsplit(sf, split = "_"))[1], sf_stats[1,1]))
    }
    colnames(time) <- c("Model", "time")
    time$time <- as.numeric(time$time)
    colnames(memory) <- c("Model", "memory")
    memory$memory <- as.numeric(memory$memory)
    
    data_re[match(metric1$model_names, models),1] <- metric1$model_ari
    data_re[match(metric1$model_names, models),2] <- metric1$model_nmi
    data_re[match(chaos$Model, models),3] <- chaos_mean
    data_re[match(pas$Model, models),4] <- pas_mean
    data_re[match(time$Model, models),5] <- time$time
    data_re[match(memory$Model, models),6] <- memory$memory
    
    data_re[match(intergration_value$Model, models),7] <- intergration_value$bASW
    data_re[match(intergration_value$Model, models),8] <- intergration_value$kBET
    data_re[match(intergration_value$Model, models),9] <- intergration_value$dASW
    data_re[match(intergration_value$Model, models),10] <- intergration_value$ilF1
    data_re[,11] <- slices_count
    
    data_re <- as.data.frame(data_re)
    colnames(data_re) <- c("ARI", "NMI", "CHAOS", "PAS", "Time", "Memory", "bASW", "kBET", "dASW", "ilF1", "slices_count")
    
    ###添加data id
    data_re$model <- models
    data_re$datasets <- d
    data_re$CellCount <- nrow(d_data)
    data_re$GeneCount <- ncol(d_data)
    data_re$AllCount <- data_re$CellCount * data_re$GeneCount
    
    #data_re$name <- models
    data_re$Platform <- c("R", "Python", "R", "Python", "Python", "Python", "Python", "R")
    data_re$Modeltype <- c("Sta", "DL", "Sta", "DL", "DL", "DL", "DL", "Sta")
    data_re <- data_re[,c("model","datasets","Platform", "Modeltype" , "ARI", "NMI", "CHAOS", "PAS", "Time", "Memory", "CellCount", "GeneCount", "AllCount", "bASW", "kBET", "dASW", "ilF1", "slices_count")]
    
    intergration_re <- rbind(intergration_re, data_re)
  }else{
    data_re <- matrix(,8,11)
    path1 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/model_metric.csv", sep = "/")
    path2 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/CHAOS.csv", sep = "/")
    path3 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/PAS.csv", sep = "/")
    d_data <- read_h5ad(paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "sample_all_data/Slices_combind_data.h5ad", sep = "/"))
    slices_count <- length(unique(d_data$obs$slices))
    
    metric1 <- read.csv(path1)
    chaos <- read.csv(path2)
    chaos_mean <- apply(chaos[,-which(colnames(chaos) == "Model")],1,mean)
    
    pas <- read.csv(path3)
    pas_mean <- apply(pas[,-which(colnames(pas) == "Model")],1,mean)
    
    stat_fl <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe", sep = "/")
    stat_flies <- list.files(stat_fl, pattern = "stats")
    time <- data.frame()
    memory <- data.frame()
    for(sf in stat_flies){
      sf_stats <- read.csv(paste(stat_fl, sf, sep = "/"))
      time <- rbind(time, c(unlist(strsplit(sf, split = "_"))[1], sf_stats[1,2]))
      memory <- rbind(memory, c(unlist(strsplit(sf, split = "_"))[1], sf_stats[1,1]))
    }
    colnames(time) <- c("Model", "time")
    time$time <- as.numeric(time$time)
    colnames(memory) <- c("Model", "memory")
    memory$memory <- as.numeric(memory$memory)
    
    data_re[match(metric1$model_names, models),1] <- metric1$model_ari
    data_re[match(metric1$model_names, models),2] <- metric1$model_nmi
    data_re[match(chaos$Model, models),3] <- chaos_mean
    data_re[match(pas$Model, models),4] <- pas_mean
    data_re[match(time$Model, models),5] <- time$time
    data_re[match(memory$Model, models),6] <- memory$memory
    
    data_re[ ,7] <- NA
    data_re[ ,8] <- NA
    data_re[ ,9] <- NA
    data_re[ ,10] <- NA
    data_re[,11] <- slices_count
    
    data_re <- as.data.frame(data_re)
    colnames(data_re) <- c("ARI", "NMI", "CHAOS", "PAS", "Time", "Memory", "bASW", "kBET", "dASW", "ilF1", "slices_count")
    
    ###添加data id
    data_re$model <- models
    data_re$datasets <- d
    data_re$CellCount <- nrow(d_data)
    data_re$GeneCount <- ncol(d_data)
    data_re$AllCount <- data_re$CellCount * data_re$GeneCount
    
    #data_re$name <- models
    data_re$Platform <- c("R", "Python", "R", "Python", "Python", "Python", "Python", "R")
    data_re$Modeltype <- c("Sta", "DL", "Sta", "DL", "DL", "DL", "DL", "Sta")
    data_re <- data_re[,c("model","datasets","Platform", "Modeltype" , "ARI", "NMI", "CHAOS", "PAS", "Time", "Memory", "CellCount", "GeneCount", "AllCount", "bASW", "kBET", "dASW", "ilF1", "slices_count")]
    
    intergration_re <- rbind(intergration_re, data_re)
  }
}

write.table(intergration_re, "/home/dongkj/home_dkj/FD_yzy/Result/1.Intergration/intergration_re_all.csv", col.names = T, row.names = F, sep = ",", quote = F)

###2. integration result plotting----
intergration_re1 <- read.csv("/home/dongkj/home_dkj/FD_yzy/Result/1.Intergration/intergration_re_all.csv")
intergration_re1[which(intergration_re1$model == "GraphSTwithPASTE"),]$model  <- "GraphST-PASTE"
intergration_re1 <- intergration_re1[intergration_re1$model != "CN",]
intergration_re1$model <- factor(intergration_re1$model, levels = c("Banksy", "CellCharter", "GraphST", "GraphST-PASTE", "MENDER", "NicheCompass", "Spado"))

#### 2.1分别绘制指标----
datasets <- unique(intergration_re1$datasets)
metric <- c("bASW", "kBET", "dASW", "ilF1")
for(dataset in datasets){
  plotData <- intergration_re1[intergration_re1$datasets == dataset,]
  plotData <- plotData[,match(c(metric,"model"), colnames(plotData))]
  plotData <- na.omit(plotData)
  p1 <- ggplot(plotData, aes(x = model, y = bASW, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = bASW), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.3, color = "#74A8D2") +
    labs(y = "bASW", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p2 <- ggplot(plotData, aes(x = model, y = kBET, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = kBET), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.3, color = "#C74B3E") +
    labs(y = "kBET", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p3 <- ggplot(plotData, aes(x = model, y = dASW, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = dASW), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.3, color = "#92c057") +
    labs(y = "dASW", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p4 <- ggplot(plotData, aes(x = model, y = ilF1, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = ilF1), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.3, color = "#d7abb6") +
    labs(y = "ilF1", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_all <- wrap_plots(list(p1,p2,p3,p4), ncol = 4)
  savePath <- paste("/NFS_home/NFS_home_6/dongkj/home_dkj/FD_yzy/Result/1.Intergration_umap/MetricPlot", paste(dataset, "pdf", sep = "."), sep = "/")
  ggsave(savePath, plot = p_all, width = 22, height = 4, units = "cm")
  ### 
}

#### 2.2 plotting function for extended fig 2----
plotFun <- function(plotData){
  # color <- c("Banksy" = "#C74B3E", "CellCharter" = "#539B4F", "GraphST" = "#ECBD7C", "GraphSTwithPASTE" = "#DC813B", "MENDER" = "#BCD3E6", "NicheCompass" ="#8481B4", "Spado" = "#B8D394" )
  color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "Spado" = "#D57BBE" )
  
  # p_bASW_bar <- ggplot(plotData, aes(x = model, y = bASW, fill = model)) +
  #   stat_summary(fun = mean, geom = "bar", alpha = 0.7, size = 0.3, width = 0.5) +  
  #   stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, size = 0.3) +
  #   geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.3), size = 0.2) +
  #   scale_fill_manual(values = color) +
  #   scale_color_manual(values = color) +
  #   theme_classic() +  
  #   labs(x = NULL, y = "bASW") +  
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
  #         axis.text.y = element_text(size = 6, color = "black"),
  #         axis.title.y = element_text(size = 6, color = "black"),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank(),
  #         axis.line = element_line(color = "black", size = 0.3)) +
  #   NoLegend() 
  p_bASW_box <- ggplot(plotData, aes(x = model, y = bASW, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "bASW") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_kBET_box <- ggplot(plotData, aes(x = model, y = kBET, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "kBET") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_dASW_box <- ggplot(plotData, aes(x = model, y = dASW, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "dASW") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_ilF1_box <- ggplot(plotData, aes(x = model, y = ilF1, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "isF1") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_box <- wrap_plots(list(p_bASW_box, p_kBET_box ,p_dASW_box, p_ilF1_box), ncol = 4)
  #p_bar <- wrap_plots(list(p_bASW_bar, p_kBET_bar ,p_dASW_bar, p_ilF1_bar), ncol = 4)
  
  return(list(boxplot =  p_box))
}

#### 2.3 DLPFC----
plotData <- intergration_re1[grep("DLPFC", intergration_re1$datasets),]
plotData <- plotData[plotData$model != "CN",]
DLPFC_plot <- plotFun(plotData)
ggsave("~/FD_yzy/Result/1.Intergration_umap/MetricPlot_with_errobar/DLPFC_box.pdf", DLPFC_plot$boxplot, width = 21, height = 5, units = "cm")

#### 2.4 MERFISH Brain----
plotData <- intergration_re1[grep("MERFISH_Brain", intergration_re1$datasets),]
plotData <- plotData[,c("model", "datasets","bASW", "kBET", "dASW", "ilF1")]
plotData <- na.omit(plotData)
MERFISH_Brain_plot <- plotFun(plotData)
ggsave("/home/dongkj/home_dkj/FD_yzy/Result/1.Intergration_umap/MetricPlot_with_errobar/MERFISH_Brain_box.pdf", MERFISH_Brain_plot$boxplot, width = 21, height = 5, units = "cm")

### 3.funck heatmap----
#### 3.1load color data-----
data("dynbenchmark_data")
palettes <- dynbenchmark_data$palettes

overall <- palettes$colours[[1]]
intergrationColor <- palettes$colours[[2]]
simulationColor <- palettes$colours[[4]]
landmarkColor <- palettes$colours[[5]]
scalabilityColor <- palettes$colours[[3]]
gradient_colors <- colorRampPalette(c("#400764", "#9F1F7A", "#E46F9F", "#F2C7C1", "#FDF5F3"))
all <- gradient_colors(100)

colours <- list(overall, intergrationColor, simulationColor, landmarkColor, scalabilityColor, all)
color <- tibble(palette = c("overall", "intergrationColor", "simulationColor", "landmarkColor", "scalabilityColor", "all"),
                colours = colours)

#### 3.2 plot heatmap for extended fig 2----
heatmapData <- intergration_re1[intergration_re1$datasets != "Mouse",match(c("model","datasets", "bASW", "kBET", "dASW", "ilF1") ,colnames(intergration_re1))]
heatmapData <- heatmapData %>%
  pivot_longer(cols = c("bASW", "kBET", "dASW", "ilF1"), 
               names_to = "Metric", 
               values_to = "Value") %>%
  unite("Metric_Dataset", datasets, Metric, sep = "_") %>%
  pivot_wider(names_from = "Metric_Dataset", values_from = "Value")
colnames(heatmapData)[1] <- "id"

heatmapData[,2:ncol(heatmapData)] <- apply(heatmapData[,2:ncol(heatmapData)], 2, function(x){
  x[which(is.na(x))] <- 0
  rank(x)
})

datasets <- c("DLPFC", "MERFISH", "BaristaSeq", "STARMap")
metric <- c("_bASW", "_kBET", "_dASW", "_ilF1")
integration_wide <- data.frame(id = heatmapData$id)
for(d in datasets){
  d_data <- heatmapData[,grep(d, colnames(heatmapData))]
  for(m in metric){
    index <- grep(m, colnames(d_data))
    d_m_data <- d_data[,index]
    d_m_data_mean <- apply(d_m_data, 1, mean)
    c_name <- paste(d,m,sep = "")
    integration_wide[[c_name]] <- d_m_data_mean
  }
}

integration_wide$id <- as.character(integration_wide$id)
#heatmapData[which(heatmapData$model == "GraphSTwithPASTE"),]$model <- "GraphST-PASTE"

averagePerformance <- data.frame(id = integration_wide$id)
for(d in datasets){
  d_data <-integration_wide[,grep(d, colnames(integration_wide))]
  d_average <- as.data.frame(apply(d_data, 1, mean))
  colnames(d_average) <- paste(d, "overall", sep = "_")
  averagePerformance <- cbind(averagePerformance, d_average)
}
heatmapData <- merge(integration_wide, averagePerformance, by = "id")

colInfo <- tribble(
  ~id,             ~group,                ~name,       ~geom,          ~palette,              ~size,   ~options,
  "id",            "Method",              "",          "text",         NA,                    1,       list(hjust = 0, width = 2),
  #"Platform",      "Meta information",    "Platform",  "text",         NA,                    1,       list(width = 1),
  #"Modeltype",     "Meta information",    "Modeltype", "text",         NA,                    1,       list(width = 1),
  "DLPFC_overall",           "10X Visium",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "DLPFC_bASW",           "10X Visium",         "bASW",       "circle",          "intergrationColor",   NA,      list(width = 1),
  "DLPFC_kBET",           "10X Visium",         "kBET",       "circle",          "intergrationColor",     NA,      list(width = 1),
  "DLPFC_dASW",         "10X Visium",         "dASW",     "circle",          "intergrationColor",       NA,      list(width = 1),
  "DLPFC_ilF1",           "10X Visium",         "ilF1",       "circle",          "intergrationColor",    NA,      list(width = 1),
  "MERFISH_overall",        "MERFISH",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "MERFISH_bASW",           "MERFISH",         "bASW",       "circle",          "simulationColor",   NA,      list(width = 1),
  "MERFISH_kBET",           "MERFISH",         "kBET",       "circle",          "simulationColor",     NA,      list(width = 1),
  "MERFISH_dASW",         "MERFISH",         "dASW",     "circle",          "simulationColor",       NA,      list(width = 1),
  "MERFISH_ilF1",           "MERFISH",         "ilF1",       "circle",          "simulationColor",    NA,      list(width = 1),
  "BaristaSeq_overall",        "BaristaSeq",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "BaristaSeq_bASW",           "BaristaSeq",         "bASW",       "circle",          "landmarkColor",   NA,      list(width = 1),
  "BaristaSeq_kBET",           "BaristaSeq",         "kBET",       "circle",          "landmarkColor",     NA,      list(width = 1),
  "BaristaSeq_dASW",         "BaristaSeq",         "dASW",     "circle",          "landmarkColor",       NA,      list(width = 1),
  "BaristaSeq_ilF1",           "BaristaSeq",         "ilF1",       "circle",          "landmarkColor",    NA,      list(width = 1),
  "STARMap_overall",           "STARMap",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "STARMap_bASW",           "STARMap",         "bASW",       "circle",          "scalabilityColor",   NA,      list(width = 1),
  "STARMap_kBET",           "STARMap",         "kBET",       "circle",          "scalabilityColor",     NA,      list(width = 1),
  "STARMap_dASW",         "STARMap",         "dASW",     "circle",          "scalabilityColor",       NA,      list(width = 1),
  "STARMap_ilF1",           "STARMap",         "ilF1",       "circle",          "scalabilityColor",    NA,      list(width = 1)
)


colGroup <- data.frame(group = c("Method", "10X Visium", "MERFISH", "BaristaSeq", "STARMap"),
                       Experiment = c("Method", "10X Visium", "MERFISH", "BaristaSeq", "STARMap"),
                       palette = rep("overall", 5)
)


p <- funky_heatmap(
  data = heatmapData,
  column_info = colInfo,
  column_groups = colGroup,
  palettes = color,
  position_args = position_arguments(col_annot_offset = 2, col_annot_angle = 30),
  scale_column = T
)
p
ggsave("/home/dongkj/home_dkj/FD_yzy/Result/1.Intergration_umap/MetricPlot_with_errobar/heatmap_all.pdf", plot = p, device = "pdf", width = 210, height = 100, units = "mm")


###4. clustering result plotting----
intergration_re1 <- read.csv("/home/dongkj/home_dkj/FD_yzy/Result/1.Intergration/intergration_re_all.csv")
intergration_re1[which(intergration_re1$model == "GraphSTwithPASTE"),]$model  <- "GraphST-PASTE"
intergration_re1$model <- factor(intergration_re1$model, levels = c("Banksy", "CellCharter", "CN", "GraphST", "GraphST-PASTE", "MENDER", "NicheCompass", "Spado"))

#### 4.1 plotting function for extended fig 3----
plotFun <- function(plotData, plotData2){
  #color <- c("Banksy" = "#C74B3E", "CellCharter" = "#539B4F", "CN" = "#5E93BF", "GraphST" = "#ECBD7C", "GraphSTwithPASTE" = "#DC813B", "MENDER" = "#BCD3E6", "NicheCompass" ="#8481B4", "Spado" = "#B8D394" )
  color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "Spado" = "#D57BBE" )
  
  plotData$Time <- log2(plotData$Time)
  plotData$Memory <- log2(plotData$Memory)
  CHAOSData <- plotData2[plotData2$Type == "CHAOS",]
  PASData <- plotData2[plotData2$Type == "PAS",]
  
  p_ARI_bar <- ggplot(plotData, aes(x = model, y = ARI, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "ARI") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_NMI_bar <- ggplot(plotData, aes(x = model, y = NMI, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "NMI") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_Time_bar <- ggplot(plotData, aes(x = model, y = Time, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "Log2(time)") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_Memory_bar <- ggplot(plotData, aes(x = model, y = Memory, fill = model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "Log2(memory)") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_CHAOS_bar <- ggplot(CHAOSData, aes(x = Model, y = Value, fill = Model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = Model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "CHAOS") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_PAS_bar <- ggplot(PASData, aes(x = Model, y = Value, fill = Model)) +
    #stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, size = 0.3)+
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = Model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "PAS") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  

  p_all <- wrap_plots(list(p_ARI_bar,p_NMI_bar,p_CHAOS_bar,p_PAS_bar), ncol = 4)
  return(p_all)
}

#### 4.2 DLPFC----
plotData <- intergration_re1[grep("DLPFC", intergration_re1$datasets),]
plotData2 <- data.frame()
datasets <- c("DLPFC_sample1", "DLPFC_sample3", "DLPFC_sample3")
for(d in datasets){
  path2 <- paste("~/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/CHAOS.csv", sep = "/")
  path3 <- paste("~/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/PAS.csv", sep = "/")
  
  chaos <- read.csv(path2)
  chaos <- chaos %>%
    pivot_longer(
      cols = setdiff(colnames(chaos), "Model"),
      names_to = "Slices",
      values_to = "Value"
    )
  chaos$Type <- "CHAOS"
  
  pas <- read.csv(path3)
  pas <- pas %>%
    pivot_longer(
      cols = setdiff(colnames(pas), "Model"),
      names_to = "Slices",
      values_to = "Value"
    )
  pas$Type <- "PAS"
  
  plotData2 <- rbind(plotData2, rbind(chaos, pas))
}
plotData2[which(plotData2$Model == "GraphSTwithPASTE"),]$Model  <- "GraphST-PASTE"
DLPFC_plot <- plotFun(plotData, plotData2)
ggsave("~/FD_yzy/Result/1.Intergration/Metric_with_errobar/DLPFC.pdf", DLPFC_plot, width = 22, height = 5, units = "cm")

#### 4.3 MERFISH_Brain----
plotData <- intergration_re1[grep("MERFISH_Brain", intergration_re1$datasets),]
plotData2 <- data.frame()
datasets <- c("MERFISH_Brain_S2", "MERFISH_Brain_S3", "MERFISH_Brain_S4", "MERFISH_Brain_S5",
              "MERFISH_Brain_S6", "MERFISH_Brain_S7", "MERFISH_Brain_S8", "MERFISH_Brain_S9",
              "MERFISH_Brain_S10", "MERFISH_Brain_S11", "MERFISH_Brain_S12")

for(d in datasets){
  path2 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/CHAOS.csv", sep = "/")
  path3 <- paste("/home/dongkj/home_dkj/FD_yzy/Dataset/Intergration_Benchmark", d, "IntergrationRe/Metric/PAS.csv", sep = "/")
  
  chaos <- read.csv(path2)
  chaos <- chaos %>%
    pivot_longer(
      cols = setdiff(colnames(chaos), "Model"),
      names_to = "Slices",
      values_to = "Value"
    )
  chaos$Type <- "CHAOS"
  
  pas <- read.csv(path3)
  pas <- pas %>%
    pivot_longer(
      cols = setdiff(colnames(pas), "Model"),
      names_to = "Slices",
      values_to = "Value"
    )
  pas$Type <- "PAS"
  
  plotData2 <- rbind(plotData2, rbind(chaos, pas))
}
plotData2[which(plotData2$Model == "GraphSTwithPASTE"),]$Model  <- "GraphST-PASTE"
MERFISH_Brain_plot <- plotFun(plotData, plotData2)
ggsave("/home/dongkj/home_dkj/FD_yzy/Result/1.Intergration/Metric_with_errobar/MERFISH_Brain.pdf", MERFISH_Brain_plot, width = 22, height = 5, units = "cm")


### 5.clustering funck heatmap----
#### 5.1 load color data-----
library(funkyheatmap)
library(kableExtra)
data("dynbenchmark_data")
palettes <- dynbenchmark_data$palettes

overall <- palettes$colours[[1]]
intergrationColor <- palettes$colours[[2]]
simulationColor <- palettes$colours[[4]]
landmarkColor <- palettes$colours[[5]]
scalabilityColor <- palettes$colours[[3]]
largeColor <- colorRampPalette(c("#370578", "#64509c", "#9c99c3", "#d9d9ea", "#fafafd"))
largeColor <- largeColor(100)
gradient_colors <- colorRampPalette(c("#400764", "#9F1F7A", "#E46F9F", "#F2C7C1", "#FDF5F3"))
all <- gradient_colors(100)

colours <- list(overall, intergrationColor, simulationColor, landmarkColor, scalabilityColor, largeColor, all)
color <- tibble(palette = c("overall", "intergrationColor", "simulationColor", "landmarkColor", "scalabilityColor", "largeColor", "all"),
                colours = colours)

#### 5.2 plot heatmap for extendex fig 3----
heatmapData <- intergration_re1[,match(c("model","datasets", "ARI", "NMI", "CHAOS", "PAS") ,colnames(intergration_re1))]
heatmapData$CHAOS <- 1-heatmapData$CHAOS
heatmapData$PAS <- 1-heatmapData$PAS

heatmapData <- heatmapData %>%
  pivot_longer(cols = c("ARI", "NMI", "CHAOS", "PAS"), 
               names_to = "Metric", 
               values_to = "Value") %>%
  unite("Metric_Dataset", datasets, Metric, sep = "_") %>%
  pivot_wider(names_from = "Metric_Dataset", values_from = "Value")
colnames(heatmapData)[1] <- "id"

heatmapData[,2:ncol(heatmapData)] <- apply(heatmapData[,2:ncol(heatmapData)], 2, function(x){
  x[which(is.na(x))] <- 0
  rank(x)
})

datasets <- c("DLPFC", "MERFISH", "BaristaSeq", "STARMap", "Mouse")
metric <- c("_ARI", "_NMI", "_CHAOS", "_PAS")
integration_wide <- data.frame(id = heatmapData$id)
for(d in datasets){
  d_data <- heatmapData[,grep(d, colnames(heatmapData))]
  for(m in metric){
    index <- grep(m, colnames(d_data))
    d_m_data <- d_data[,index]
    d_m_data_mean <- apply(d_m_data, 1, mean)
    c_name <- paste(d,m,sep = "")
    integration_wide[[c_name]] <- d_m_data_mean
  }
}

integration_wide$id <- as.character(integration_wide$id)

averagePerformance <- data.frame(id = integration_wide$id)
for(d in datasets){
  d_data <-integration_wide[,grep(d, colnames(integration_wide))]
  d_average <- as.data.frame(apply(d_data, 1, mean))
  colnames(d_average) <- paste(d, "overall", sep = "_")
  averagePerformance <- cbind(averagePerformance, d_average)
}
heatmapData <- merge(integration_wide, averagePerformance, by = "id")


colInfo <- tribble(
  ~id,             ~group,                ~name,       ~geom,          ~palette,              ~size,   ~options,
  "id",            "Method",              "",          "text",         NA,                    1,       list(hjust = 0, width = 2),
  #"Platform",      "Meta information",    "Platform",  "text",         NA,                    1,       list(width = 1),
  #"Modeltype",     "Meta information",    "Modeltype", "text",         NA,                    1,       list(width = 1),
  "DLPFC_overall",           "10X Visium",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "DLPFC_ARI",           "10X Visium",         "ARI",       "circle",          "intergrationColor",   NA,      list(width = 1),
  "DLPFC_NMI",           "10X Visium",         "NMI",       "circle",          "intergrationColor",     NA,      list(width = 1),
  "DLPFC_CHAOS",         "10X Visium",         "CHAOS",     "circle",          "intergrationColor",       NA,      list(width = 1),
  "DLPFC_PAS",           "10X Visium",         "PAS",       "circle",          "intergrationColor",    NA,      list(width = 1),
  "MERFISH_overall",        "MERFISH",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "MERFISH_ARI",           "MERFISH",         "ARI",       "circle",          "simulationColor",   NA,      list(width = 1),
  "MERFISH_NMI",           "MERFISH",         "NMI",       "circle",          "simulationColor",     NA,      list(width = 1),
  "MERFISH_CHAOS",         "MERFISH",         "CHAOS",     "circle",          "simulationColor",       NA,      list(width = 1),
  "MERFISH_PAS",           "MERFISH",         "PAS",       "circle",          "simulationColor",    NA,      list(width = 1),
  "BaristaSeq_overall",        "BaristaSeq",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "BaristaSeq_ARI",           "BaristaSeq",         "ARI",       "circle",          "landmarkColor",   NA,      list(width = 1),
  "BaristaSeq_NMI",           "BaristaSeq",         "NMI",       "circle",          "landmarkColor",     NA,      list(width = 1),
  "BaristaSeq_CHAOS",         "BaristaSeq",         "CHAOS",     "circle",          "landmarkColor",       NA,      list(width = 1),
  "BaristaSeq_PAS",           "BaristaSeq",         "PAS",       "circle",          "landmarkColor",    NA,      list(width = 1),
  "STARMap_overall",           "STARMap",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "STARMap_ARI",           "STARMap",         "ARI",       "circle",          "scalabilityColor",   NA,      list(width = 1),
  "STARMap_NMI",           "STARMap",         "NMI",       "circle",          "scalabilityColor",     NA,      list(width = 1),
  "STARMap_CHAOS",         "STARMap",         "CHAOS",     "circle",          "scalabilityColor",       NA,      list(width = 1),
  "STARMap_PAS",           "STARMap",         "PAS",       "circle",          "scalabilityColor",    NA,      list(width = 1),
  
  "Mouse_overall",           "Large data",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "Mouse_ARI",           "Large data",         "ARI",       "circle",          "largeColor",   NA,      list(width = 1),
  "Mouse_NMI",           "Large data",         "NMI",       "circle",          "largeColor",     NA,      list(width = 1),
  "Mouse_CHAOS",         "Large data",         "CHAOS",     "circle",          "largeColor",       NA,      list(width = 1),
  "Mouse_PAS",           "Large data",         "PAS",       "circle",          "largeColor",    NA,      list(width = 1)
)


colGroup <- data.frame(group = c("Method", "10X Visium", "MERFISH", "BaristaSeq", "STARMap", "Large data"),
                       Experiment = c("Method", "10X Visium", "MERFISH", "BaristaSeq", "STARMap", "Large data"),
                       palette = rep("overall", 6)
)


p <- funky_heatmap(
  data = heatmapData,
  column_info = colInfo,
  column_groups = colGroup,
  palettes = color,
  position_args = position_arguments(col_annot_offset = 2, col_annot_angle = 30),
  scale_column = T
)
p
ggsave("/home/dongkj/home_dkj/FD_yzy/Result/1.Intergration/Metric_with_errobar/heatmap_all_new.pdf", plot = p, device = "pdf", width = 210, height = 100, units = "mm")


###6. time plot----
intergration_re <- read.csv("/NFS_home/NFS_home_6/home_dkj/FD_yzy/Result/1.Intergration/intergration_re_all.csv")
intergration_re[which(intergration_re$model == "GraphSTwithPASTE"),]$model  <- "GraphST-PASTE"
intergration_re$model <- factor(intergration_re$model, levels = c("Banksy", "CellCharter", "CN", "GraphST", "GraphST-PASTE", "MENDER", "NicheCompass", "Spado"))

color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "Spado" = "#D57BBE" )

plotData <- intergration_re[,c(1,2,9,10,11,12,13)]
plotData <- na.omit(plotData)
plotData$Time <- log2(plotData$Time)
plotData$Memory <- log2(plotData$Memory)

p_time_1 <- ggplot(plotData, aes(x = CellCount, y = Time, color = model, group = model)) +
  geom_line(size = 0.5) + 
  geom_point(size = 1) +  
  labs(x = "log10(Cell Counts)", 
       y = "Log2(Time(s))") +
  scale_color_manual(values = color)+
  scale_x_log10() +
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        legend.position = "none") 

p_time_2 <- ggplot(plotData, aes(x = GeneCount, y = Time, color = model, group = model)) +
  geom_line(size = 0.5) + 
  geom_point(size = 1) + 
  labs(x = "log10(Gene Counts)", 
       y = "Log2(Time(s))") + 
  scale_color_manual(values = color)+
  scale_x_log10() + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        legend.position = "none") 

p_time_3 <- ggplot(plotData, aes(x = AllCount, y = Time, color = model, group = model)) +
  geom_line(size = 0.5) + 
  geom_point(size = 1) +  
  labs(x = "log10(All Counts)", 
       y = "Log2(Time(s))") +
  scale_color_manual(values = color)+
  scale_x_log10() + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        legend.position = "none") 

p_memory_1 <- ggplot(plotData, aes(x = CellCount, y = Memory, color = model, group = model)) +
  geom_line(size = 0.5) + 
  geom_point(size = 1) +  
  labs(x = "log10(Cell Counts)", 
       y = "Log(Memory(Mb))") +  
  scale_color_manual(values = color)+
  scale_x_log10() + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        legend.position = "none") 

p_memory_2 <- ggplot(plotData, aes(x = GeneCount, y = Memory, color = model, group = model)) +
  geom_line(size = 0.5) + 
  geom_point(size = 1) + 
  labs(x = "log10(Gene Counts)", 
       y = "Log(Memory(Mb))") + 
  scale_color_manual(values = color)+
  scale_x_log10() + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        legend.position = "none")

p_memory_3 <- ggplot(plotData, aes(x = AllCount, y = Memory, color = model, group = model)) +
  geom_line(size = 0.5) + 
  geom_point(size = 1) + 
  labs(x = "log10(All Counts)", 
       y = "Log(Memory(Mb))") +
  scale_color_manual(values = color)+
  scale_x_log10() + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        legend.position = "none")

p_all <- wrap_plots(list(p_time_1, p_memory_1,
                         p_time_2, p_memory_2,
                         p_time_3, p_memory_3), ncol = 2)
ggsave("~/FD_yzy/Result/1.Intergration/Metric_with_errobar/time_memory.pdf", plot = p_all, width = 20, height = 15, units = "cm")


