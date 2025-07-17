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


setwd("iSTBench/")

###1.read data----
### Intergration mertic value
# The consolidated results from all datasets are read, and the results are merged into a final file for visualization. 
# The related result files, integration_re.csv and integration_re_all.csv, can be directly used in subsequent analyses 
# without needing to repeat this step.
models <- c("Banksy", "CellCharter", "CN", "GraphST", "GraphSTwithPASTE", "MENDER", "NicheCompass", "PRECAST", "Spado", "SPIRAL", "STAIG", "STAligner")

intergration_re <- data.frame()
for(d in datasets){
  data_re <- matrix(,12,13)
  path1 <- paste("Data", d, "IntergrationRe/Metric/model_metric.csv", sep = "/")
  path2 <- paste("Data", d, "IntergrationRe/Metric/CHAOS.csv", sep = "/")
  path3 <- paste("Data", d, "IntergrationRe/Metric/PAS.csv", sep = "/")
  path4 <- paste("Data", d, "IntergrationRe/Metric/Intergration_value.csv", sep = "/")
  
  d_data <- read_h5ad(paste("Data", d, "sample_all_data/Slices_combind_data.h5ad", sep = "/"))
  slices_count <- length(unique(d_data$obs$slices))
  
  metric1 <- read.csv(path1)
  chaos <- read.csv(path2)
  chaos_mean <- apply(chaos[,-which(colnames(chaos) == "Model")],1,mean)
  
  pas <- read.csv(path3)
  pas_mean <- apply(pas[,-which(colnames(pas) == "Model")],1,mean)
  
  intergration_value <- read.csv(path4)
  
  stat_fl <- paste("Data", d, "IntergrationRe", sep = "/")
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
  
  data_re[match(intergration_value$Model, models),7] <- intergration_value$dASW
  data_re[match(intergration_value$Model, models),8] <- intergration_value$dLISI
  data_re[match(intergration_value$Model, models),9] <- intergration_value$ILL
  data_re[match(intergration_value$Model, models),10] <- intergration_value$bASW
  data_re[match(intergration_value$Model, models),11] <- intergration_value$iLISI
  data_re[match(intergration_value$Model, models),12] <- intergration_value$GC
  data_re[,13] <- slices_count
  
  data_re <- as.data.frame(data_re)
  colnames(data_re) <- c("ARI", "NMI", "CHAOS", "PAS", "Time", "Memory", "dASW", "dLISI", "ILL", "bASW", "iLISI", "GC", "slices_count")
  
  data_re$model <- models
  data_re$datasets <- d
  data_re$CellCount <- nrow(d_data)
  data_re$GeneCount <- ncol(d_data)
  data_re$AllCount <- data_re$CellCount * data_re$GeneCount
  
  data_re$Platform <- c("R", "Python", "R", "Python", "Python", "Python", "Python", "R", "R", "Python", "Python", "Python")
  data_re$Modeltype <- c("Sta", "DL", "Sta", "DL", "DL", "DL", "DL", "Sta", "Sta", "DL", "DL", "DL")
  data_re <- data_re[,c("model","datasets","Platform", "Modeltype" , "ARI", "NMI", "CHAOS", "PAS", "Time", "Memory", "CellCount", "GeneCount", "AllCount", "dASW", "dLISI", "ILL", "bASW","iLISI","GC", "slices_count")]
  
  intergration_re <- rbind(intergration_re, data_re)
}
write.table(intergration_re, "Analysis/Integration_and_clustering/result/intergration_re_all.csv", col.names = T, row.names = F, sep = ",", quote = F)


setwd("/iSTBench")
###2. integration result plotting----
intergration_re1 <- read.csv("Analysis/Integration_and_clustering/result/intergration_re_all.csv")
intergration_re1[which(intergration_re1$model == "GraphSTwithPASTE"),]$model  <- "GraphST-PASTE"
intergration_re1 <- intergration_re1[intergration_re1$model != "CN",]
intergration_re1$model <- factor(intergration_re1$model, levels = c("Banksy", "CellCharter", "GraphST", "GraphST-PASTE", "MENDER", "NicheCompass", "Spado"))

#### 2.1 Draw indicators respectively----
datasets <- unique(intergration_re1$datasets)
metric <- c("dASW", "dLISI", "ILL", "bASW", "iLISI", "GC")
for(dataset in datasets){
  plotData <- intergration_re1[intergration_re1$datasets == dataset,]
  plotData <- plotData[,match(c(metric,"model"), colnames(plotData))]
  plotData <- na.omit(plotData)
  
  p1 <- ggplot(plotData, aes(x = model, y = dASW, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = dASW), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.1, color = "#74A8D2") +
    labs(y = "dASW", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p2 <- ggplot(plotData, aes(x = model, y = dLISI, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = dLISI), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.1, color = "#92c057") +
    labs(y = "dLISI", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p3 <- ggplot(plotData, aes(x = model, y = ILL, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = ILL), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.1, color = "#58BACC") +
    labs(y = "ILL", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p4 <- ggplot(plotData, aes(x = model, y = bASW, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = bASW), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.1, color = "#C74B3E") +
    labs(y = "bASW", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p5 <- ggplot(plotData, aes(x = model, y = iLISI, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = iLISI), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.1, color = "#8D68B8") +
    labs(y = "iLISI", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p6 <- ggplot(plotData, aes(x = model, y = GC, col = model)) +
    geom_segment(aes(x = model, xend = model, y = 0, yend = GC), 
                 arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.1, color = "#d7abb6") +
    labs(y = "GC", x  = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_all <- wrap_plots(list(p1,p2,p3,p4,p5,p6), ncol = 6)
  savePath <- paste("Analysis/Integration_and_clustering/result", paste(dataset, "pdf", sep = "."), sep = "/")
  ggsave(savePath, plot = p_all, width = 22, height = 4, units = "cm")
  ### 
}

#### 2.2 plotting function for fig 2----
plotFun <- function(plotData){
  color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8",
             "NicheCompass" ="#BABD45", "PRECAST" = "#82564c", "Spado" = "#D57BBE", "SPIRAL" = "#b3c6e3", "STAIG" = "#f19c97", "STAligner" = "#a7dc92")
  
  p_dASW_box <- ggplot(plotData, aes(x = model, y = dASW, fill = model)) +
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
  
  p_dLISI_box <- ggplot(plotData, aes(x = model, y = dLISI, fill = model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "dLISI") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_ILL_box <- ggplot(plotData, aes(x = model, y = ILL, fill = model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "ILL") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_bASW_box <- ggplot(plotData, aes(x = model, y = bASW, fill = model)) +
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
  
  p_iLISI_box <- ggplot(plotData, aes(x = model, y = iLISI, fill = model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "iLISI") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_GC_box <- ggplot(plotData, aes(x = model, y = GC, fill = model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.3, width = 0.5) +  
    geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
                size = 0.2) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = color)+
    theme_classic() +  
    labs(x = NULL, y = "GC") +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.3))+
    NoLegend()
  
  p_box <- wrap_plots(list(p_dASW_box, p_dLISI_box ,p_ILL_box, p_bASW_box, p_iLISI_box, p_GC_box), ncol = 6)
  
  return(list(boxplot =  p_box))
}

#### 2.3 DLPFC----
plotData <- intergration_re1[grep("DLPFC", intergration_re1$datasets),]
plotData <- plotData[plotData$model != "NicheCompass",]
DLPFC_plot <- plotFun(plotData)
ggsave("Analysis/Integration_and_clustering/result/DLPFC_box.pdf", DLPFC_plot$boxplot, width = 21, height = 5, units = "cm")

### 3.funck heatmap----
#### 3.1load color data-----
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

#### 3.2 plot heatmap for fig 2----
heatmapData <- intergration_re1[,match(c("model","datasets", "dASW", "dLISI", "ILL", "bASW", "iLISI", "GC") ,colnames(intergration_re1))]
heatmapData <- heatmapData %>%
  pivot_longer(cols = c("dASW", "dLISI", "ILL", "bASW", "iLISI", "GC"), 
               names_to = "Metric", 
               values_to = "Value") %>%
  unite("Metric_Dataset", datasets, Metric, sep = "_") %>%
  pivot_wider(names_from = "Metric_Dataset", values_from = "Value")
colnames(heatmapData)[1] <- "id"

heatmapData[,2:ncol(heatmapData)] <- apply(heatmapData[,2:ncol(heatmapData)], 2, function(x){
  x[which(is.na(x))] <- 0
  rank_custom_jump(x)
})

datasets <- c("DLPFC", "MERFISH", "BaristaSeq", "STARMap", "Mouse")
metric <- c("_dASW", "_dLISI", "_ILL", "_bASW", "_iLISI", "_GC")

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
  "DLPFC_overall",           "10X Visium",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "DLPFC_dASW",           "10X Visium",         "dASW",       "circle",          "intergrationColor",   NA,      list(width = 1),
  "DLPFC_dLISI",           "10X Visium",         "dLISI",       "circle",          "intergrationColor",     NA,      list(width = 1),
  "DLPFC_ILL",         "10X Visium",         "ILL",     "circle",          "intergrationColor",       NA,      list(width = 1),
  "DLPFC_bASW",           "10X Visium",         "bASW",       "circle",          "intergrationColor",    NA,      list(width = 1),
  "DLPFC_iLISI",           "10X Visium",         "iLISI",       "circle",          "intergrationColor",    NA,      list(width = 1),
  "DLPFC_GC",           "10X Visium",         "GC",       "circle",          "intergrationColor",    NA,      list(width = 1),
  "MERFISH_overall",           "MERFISH",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "MERFISH_dASW",           "MERFISH",         "dASW",       "circle",          "simulationColor",   NA,      list(width = 1),
  "MERFISH_dLISI",           "MERFISH",         "dLISI",       "circle",          "simulationColor",     NA,      list(width = 1),
  "MERFISH_ILL",         "MERFISH",         "ILL",     "circle",          "simulationColor",       NA,      list(width = 1),
  "MERFISH_bASW",           "MERFISH",         "bASW",       "circle",          "simulationColor",    NA,      list(width = 1),
  "MERFISH_iLISI",           "MERFISH",         "iLISI",       "circle",          "simulationColor",    NA,      list(width = 1),
  "MERFISH_GC",           "MERFISH",         "GC",       "circle",          "simulationColor",    NA,      list(width = 1),
  "BaristaSeq_overall",           "BaristaSeq",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "BaristaSeq_dASW",           "BaristaSeq",         "dASW",       "circle",          "landmarkColor",   NA,      list(width = 1),
  "BaristaSeq_dLISI",           "BaristaSeq",         "dLISI",       "circle",          "landmarkColor",     NA,      list(width = 1),
  "BaristaSeq_ILL",         "BaristaSeq",         "ILL",     "circle",          "landmarkColor",       NA,      list(width = 1),
  "BaristaSeq_bASW",           "BaristaSeq",         "bASW",       "circle",          "landmarkColor",    NA,      list(width = 1),
  "BaristaSeq_iLISI",           "BaristaSeq",         "iLISI",       "circle",          "landmarkColor",    NA,      list(width = 1),
  "BaristaSeq_GC",           "BaristaSeq",         "GC",       "circle",          "landmarkColor",    NA,      list(width = 1),
  "STARMap_overall",           "STARMap",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "STARMap_dASW",           "STARMap",         "dASW",       "circle",          "scalabilityColor",   NA,      list(width = 1),
  "STARMap_dLISI",           "STARMap",         "dLISI",       "circle",          "scalabilityColor",     NA,      list(width = 1),
  "STARMap_ILL",         "STARMap",         "ILL",     "circle",          "scalabilityColor",       NA,      list(width = 1),
  "STARMap_bASW",           "STARMap",         "bASW",       "circle",          "scalabilityColor",    NA,      list(width = 1),
  "STARMap_iLISI",           "STARMap",         "iLISI",       "circle",          "scalabilityColor",    NA,      list(width = 1),
  "STARMap_GC",           "STARMap",         "GC",       "circle",          "scalabilityColor",    NA,      list(width = 1),
  "Mouse_overall",           "Large data",         "Overall",       "bar",          "all",   NA,      list(width = 3),
  "Mouse_dASW",           "Large data",         "dASW",       "circle",          "largeColor",   NA,      list(width = 1),
  "Mouse_dLISI",           "Large data",         "dLISI",       "circle",          "largeColor",     NA,      list(width = 1),
  "Mouse_ILL",         "Large data",         "ILL",     "circle",          "largeColor",       NA,      list(width = 1),
  "Mouse_bASW",           "Large data",         "bASW",       "circle",          "largeColor",    NA,      list(width = 1),
  "Mouse_iLISI",           "Large data",         "iLISI",       "circle",          "largeColor",    NA,      list(width = 1),
  "Mouse_GC",           "Large data",         "GC",       "circle",          "largeColor",    NA,      list(width = 1)
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
ggsave("Analysis/Integration_and_clustering/result/Integration_heatmap.pdf", plot = p, device = "pdf", width = 210, height = 100, units = "mm")


###4. clustering result plotting----
intergration_re1 <- read.csv("Analysis/Integration_and_clustering/result/intergration_re_all.csv")
intergration_re1[which(intergration_re1$model == "GraphSTwithPASTE"),]$model  <- "GraphST-PASTE"
intergration_re1$model <- factor(intergration_re1$model, levels = c("Banksy", "CellCharter", "CN", "GraphST", "GraphST-PASTE", "MENDER", "NicheCompass", "PRECAST",
                                                                    "Spado", "SPIRAL", "STAIG", "STAligner"))
#### 4.1 plotting function for fig 3----
plotFun <- function(plotData, plotData2){
  color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8",
             "NicheCompass" ="#BABD45", "PRECAST" = "#82564c", "Spado" = "#D57BBE", "SPIRAL" = "#b3c6e3", "STAIG" = "#f19c97", "STAligner" = "#a7dc92")
  
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
# It is necessary to complete the corresponding integration analysis and 
# preliminary result analysis of the data first
plotData <- intergration_re1[grep("DLPFC", intergration_re1$datasets),]
plotData <- plotData[plotData$model != "NicheCompass",]
plotData2 <- data.frame()
datasets <- c("DLPFC_sample1", "DLPFC_sample3", "DLPFC_sample3")
for(d in datasets){
  path2 <- paste("Data", d, "IntergrationRe/Metric/CHAOS.csv", sep = "/")
  path3 <- paste("Data", d, "IntergrationRe/Metric/PAS.csv", sep = "/")
  
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
ggsave("Analysis/Integration_and_clustering/result/DLPFC_clustering.pdf", DLPFC_plot, width = 22, height = 5, units = "cm")

### 5.clustering funck heatmap----
#### 5.1 load color data-----
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

#### 5.2 plot heatmap for fig 3----
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
  x[which(is.na(x))] <- -10
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
  d_average <- as.data.frame(apply(d_data, 1, function(x){
    return((0.8*(x[1] + x[2]) + 0.2*(x[3] + x[4])) / 2)
  }))
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
ggsave("Analysis/Integration_and_clustering/result/Clustering_heatmap.pdf", plot = p, device = "pdf", width = 210, height = 100, units = "mm")


###6. time plot----
intergration_re <- read.csv("Analysis/Integration_and_clustering/result/intergration_re_all.csv")
intergration_re[which(intergration_re$model == "GraphSTwithPASTE"),]$model  <- "GraphST-PASTE"
intergration_re$model <- factor(intergration_re$model, levels = c("Banksy", "CellCharter", "CN", "GraphST", "GraphST-PASTE", "MENDER", "NicheCompass", "PRECAST","Spado", "SPIRAL", "STAIG", "STAligner"))

###time plot----
color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8",
           "NicheCompass" ="#BABD45", "PRECAST" = "#82564c", "Spado" = "#D57BBE", "SPIRAL" = "#b3c6e3", "STAIG" = "#f19c97", "STAligner" = "#a7dc92")

plotData <- intergration_re[,c(1,2,9,10,11,12,13)]
plotData <- na.omit(plotData)
plotData$Time <- log2(plotData$Time)
plotData$Memory <- log2(plotData$Memory)

p_time_1 <- ggplot(plotData, aes(x = CellCount, y = Time, color = model, group = model)) +
  geom_line(size = 0.5) + 
  geom_point(size = 1) +  
  labs(x = "Cell Counts", 
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
  labs(x = "Gene Counts", 
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
  labs(x = "All Counts", 
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
  labs(x = "Cell Counts", 
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
  labs(x = "Gene Counts", 
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
  labs(x = "All Counts", 
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
ggsave("Analysis/Integration_and_clustering/result/time_memory.pdf", plot = p_all, width = 20, height = 15, units = "cm")

