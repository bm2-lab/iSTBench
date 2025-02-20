library(tidyverse)
library(tibble)
library(pheatmap)
library(anndata)
library(dplyr)
library(ggplot2)

setwd("~/home_dkj/FD_yzy/Result/GitHub_test/iSTBench")

###1.read data----
###1.1 Intergration value
# The consolidated results from all datasets are read, and the results are merged into a final file for visualization. 
# The related result files, intergration_re.csv and intergration_sd.csv, can be directly used in subsequent analyses 
# without needing to repeat this step.
datasets <- c("DLPFC_sample3", "MERFISH", "MERFISH_Brain_S3","BaristaSeq", "STARMap")
models <- c("Banksy", "CellCharter", "CN", "GraphST", "GraphSTwithPASTE", "MENDER", "NicheCompass", "Spado", "unCorrect")

intergration_re <- data.frame()
for(d in datasets){
  data_re <- matrix(,length(models),2)

  # Read accuracy and ratio metrics from CSV files for each dataset
  path1 <- paste("/home/dongkj/home_dkj/FD_yzy/Landmark_based_intergration/Dataset", d, "Metric/domain_match_rate.csv", sep = "/")
  path2 <- paste("/home/dongkj/home_dkj/FD_yzy/Landmark_based_intergration/Dataset", d, "Metric/domain_match_ratio.csv", sep = "/")
  # path1 <- paste("Benchmark/Alignment/Result", d, "Metric/Accuracy.csv", sep = "/")
  # path2 <- paste("Benchmark/Alignment/Result", d, "Metric/Ratio.csv", sep = "/")
  
  if(d == "MERFISH_Brain_S3"){
    Accuracy <-read.csv(path1)
    data_re[match(Accuracy$model, models),1] <- Accuracy[,1]
    Ratio <-read.csv(path2)
    Ratio[,1] <- abs(log2(Ratio[,1]))
    data_re[match(Ratio$model, models),2] <- Ratio[,1]
    
  }else{
    Accuracy <-read.csv(path1)
    rownames(Accuracy) <- Accuracy$model
    Accuracy <- Accuracy[,-match("model", colnames(Accuracy))]
    Accuracy <- Accuracy[models, ]
    Accuracy2 <- apply(Accuracy, 1, mean)
    
    Ratio <- read.csv(path2)
    rownames(Ratio) <- Ratio$model
    Ratio <- Ratio[,-match("model", colnames(Ratio))]
    Ratio <- Ratio[models, ]
    Ratio <- abs(log2(Ratio))
    Ratio2 <- apply(Ratio, 1, mean)
    
    data_re[,1] <- Accuracy2
    data_re[,2] <- Ratio2
  }
  
  data_re <- as.data.frame(data_re)
  colnames(data_re) <- c("Accuracy", "Ratio")
  data_re$model <- models
  data_re$datasets <- d
  intergration_re <- rbind(intergration_re, data_re)
}

# Similar data processing steps are applied for the standard deviation calculations
intergration_sd <- data.frame()
for(d in datasets){
  data_re <- matrix(,length(models),4)
  path1 <- paste("/home/dongkj/home_dkj/FD_yzy/Landmark_based_intergration/Dataset", d, "Metric/domain_match_rate.csv", sep = "/")
  path2 <- paste("/home/dongkj/home_dkj/FD_yzy/Landmark_based_intergration/Dataset", d, "Metric/domain_match_ratio.csv", sep = "/")
  # path1 <- paste("Benchmark/Alignment/Result", d, "Metric/Accuracy.csv", sep = "/")
  # path2 <- paste("Benchmark/Alignment/Result", d, "Metric/Ratio.csv", sep = "/")
  
  if(d == "MERFISH_Brain_S3"){
    Accuracy <-read.csv(path1)
    data_re[,2] <- 0
    data_re[match(Accuracy$model, models),1] <- Accuracy$s1.s2
    
    Ratio <-read.csv(path2)
    Ratio[,1] <- abs(log2(Ratio[,1]))
    data_re[match(Ratio$model, models),3] <- Ratio$s1.s2
    data_re[,4] <- 0
    
  }else{
    Accuracy <-read.csv(path1)
    rownames(Accuracy) <- Accuracy$model
    Accuracy <- Accuracy[,-match("model", colnames(Accuracy))]
    Accuracy <- Accuracy[models, ]
    Accuracy_mean <- apply(Accuracy, 1, mean)
    Accuracy_sd <- apply(Accuracy, 1, sd)
    
    Ratio <- read.csv(path2)
    rownames(Ratio) <- Ratio$model
    Ratio <- Ratio[,-match("model", colnames(Ratio))]
    Ratio <- Ratio[models, ]
    Ratio <- abs(log2(Ratio))
    Ratio_mean <- apply(Ratio, 1, mean)
    Ratio_sd <- apply(Ratio, 1, sd)
    
    data_re[,1] <- Accuracy_mean
    data_re[,2] <- Accuracy_sd
    data_re[,3] <- Ratio_mean
    data_re[,4] <- Ratio_sd
  }
  
  data_re <- as.data.frame(data_re)
  colnames(data_re) <- c("Accuracy", "Accuracy_sd", "Ratio", "Ratio_sd")
  data_re$model <- models
  data_re$datasets <- d
  intergration_sd <- rbind(intergration_sd, data_re)
}

# Save the integration result and standard deviation to CSV files
write.table(intergration_re, "Analysis/Spatial_alignment/result/intergration_re.csv", col.names = T, row.names = F, sep = ",", quote = F)
write.table(intergration_sd, "Analysis/Spatial_alignment/result/intergration_sd.csv", col.names = T, row.names = F, sep = ",", quote = F)

#intergration_re <- read.csv("Analysis/Spatial_alignment/result/intergration_re.csv")
#intergration_sd <- read.csv("Analysis/Spatial_alignment/result/intergration_sd.csv")


### 2. Point with Error Bar ----
# This section generates scatter plots with error bars for each dataset. 
# The plot shows the relationship between accuracy and ratio, 
# with error bars representing the standard deviation for each model and dataset.

color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "Spado" = "#D57BBE" )
datasets <- c("DLPFC_sample3", "MERFISH", "MERFISH_Brain_S3", "BaristaSeq", "STARMap")
plotlist <- list()
for(d in datasets){
  plotData <- intergration_sd[intergration_sd$datasets == d,]
  plotData <- plotData[!(plotData$model == "unCorrect" | plotData$model == "GraphSTwithPASTE"),]
  plotData$model <- factor(plotData$model, levels = c("Banksy", "CellCharter", "CN", "GraphST", 
                                                      "MENDER", "NicheCompass", "Spado"))
  plotData <- na.omit(plotData)
  p <- ggplot(plotData, aes(x = Accuracy, y = Ratio, label = model, color = model)) +
    geom_point(size = 0.7) +  # 绘制点
    geom_errorbar(aes(ymin = Ratio - Ratio_sd, ymax = Ratio + Ratio_sd), width = 0.05, size = 0.2) +
    geom_errorbarh(aes(xmin = Accuracy - Accuracy_sd, xmax = Accuracy + Accuracy_sd), height = 0.05, size = 0.2) +
    scale_color_manual(values = color) +  
    labs(x = "Accuracy", y = "Ratio") + 
    theme_classic() +  
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.2),
          legend.position = "none")
  plotlist[[d]] <- p
}
p_all <- wrap_plots(plotlist, ncol = 5)
ggsave("Analysis/Spatial_alignment/result/Metric.pdf", p_all, width = 22, height = 4.5, units = "cm")

### 3. Funky Heatmap ----
# This section creates a heatmap of the integration results across all datasets. 
# The heatmap visualizes the accuracy and ratio metrics, with color scales 
# customized for different datasets. The results are displayed in a visually 
# informative format, highlighting the performance of different models across datasets.

library(funkyheatmap)
library(kableExtra)
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

heatmapData <- intergration_re[,match(c("model", "datasets", "Accuracy", "Ratio") ,colnames(intergration_re))]
heatmapData <- heatmapData[heatmapData$datasets %in% c("DLPFC_sample3", "MERFISH", "MERFISH_Brain_S3", "BaristaSeq", "STARMap"),]
heatmapData$model <- as.character(heatmapData$model)
heatmapData$Ratio <- (100-heatmapData$Ratio)

heatmapData <- heatmapData %>%
  pivot_longer(cols = c("Accuracy", "Ratio"), 
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
metric <- c("_Accuracy", "_Ratio")
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
  "DLPFC_Accuracy",           "10X Visium",         "Accuracy",       "circle",          "intergrationColor",   NA,      list(width = 2),
  "DLPFC_Ratio",           "10X Visium",         "Ratio",       "circle",          "intergrationColor",     NA,      list(width = 2),
  
  "MERFISH_overall",        "MERFISH",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "MERFISH_Accuracy",           "MERFISH",         "Accuracy",       "circle",          "simulationColor",   NA,      list(width = 2),
  "MERFISH_Ratio",           "MERFISH",         "Ratio",       "circle",          "simulationColor",     NA,      list(width = 2),
  
  "BaristaSeq_overall",        "BaristaSeq",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "BaristaSeq_Accuracy",           "BaristaSeq",         "Accuracy",       "circle",          "landmarkColor",   NA,      list(width = 2),
  "BaristaSeq_Ratio",           "BaristaSeq",         "Ratio",       "circle",          "landmarkColor",     NA,      list(width = 2),
  
  "STARMap_overall",           "STARMap",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "STARMap_Accuracy",           "STARMap",         "Accuracy",       "circle",          "scalabilityColor",   NA,      list(width = 2),
  "STARMap_Ratio",           "STARMap",         "Ratio",       "circle",          "scalabilityColor",     NA,      list(width = 2)
)


colGroup <- data.frame(group = c("Method", "10X Visium", "MERFISH", "BaristaSeq", "STARMap"),
                       Experiment = c("Method", "10X Visium", "MERFISH", "BaristaSeq", "STARMap" ),
                       palette = rep("overall", 5)
)


p <- funky_heatmap(
  data = heatmapData,
  column_info = colInfo,
  column_groups = colGroup,
  palettes = color,
  position_args = position_arguments(col_annot_offset = 2, col_annot_angle = 0),
  scale_column = T
)
p
ggsave("Analysis/Spatial_alignment/result/heatmap_all_new.pdf", plot = p, device = "pdf", width = 210, height = 100, units = "mm")

