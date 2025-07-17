# Load required libraries
library(Seurat)
library(magrittr)
library(SeuratDisk)
library(anndata)
library(getopt)
library(reshape2) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(scales) 
library(patchwork)
# Set working directory for storing results
setwd("/NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub_test/iSTBench")

## hclust----
### 1. Read data and prepare classification metrics ----
# This section reads metric data for different models and combines them into a single data frame.
# It calculates ARI (Adjusted Rand Index) and NMI (Normalized Mutual Information) for each model.
models <- c("Banksy", "CellCharter", "CN", "MENDER", "NicheCompass", "PRECAST")

# Initialize an empty dataframe to store classification metrics for all models
classiMetricRe <- data.frame()
for(m in models){
  file <- paste("Data/TNBC/SlicesEmbedding/", m, "/Metric/hclust/model_metric.csv", sep = "")

  re <- read.csv(file)
  colnames(re) <- c("ncluster","ARI", "NMI", "ncluster2")
  re$Model <- m
  classiMetricRe <- rbind(classiMetricRe, re)
}
write.table(classiMetricRe, "Analysis/Slices_representation/result/representation_re.csv", col.names = T, row.names = F, quote = F, sep = ",")

classiMetricRe2 <- classiMetricRe
classiMetricRe <- classiMetricRe[,-1]

# Calculate mean, standard deviation, and coefficient of variation for ARI and NMI per model
classiMetricMean <- classiMetricRe %>%
  group_by(Model) %>%
  summarise(
    ARI_mean = mean(ARI),
    ARI_sd = sd(ARI),
    NMI_mean = mean(NMI),
    NMI_sd = sd(NMI)
  )

# Calculate the coefficient of variation for ARI and NMI
classiMetricMean$ARI_CV <- classiMetricMean$ARI_sd/classiMetricMean$ARI_mean
classiMetricMean$NMI_CV <- classiMetricMean$NMI_sd/classiMetricMean$NMI_mean
classiMetricMean <- classiMetricMean[match(models, classiMetricMean$Model),]


### 2. Data preparation for plotting ----
# This section reshapes the data for plotting and calculates the metrics for Mean, SD, and CV.
# The data is then written to a CSV file for further use in visualization.
plotData <- classiMetricRe %>%
  pivot_wider(
    names_from = Model,
    values_from = c(ARI, NMI),
    names_sep = "_"
  )

# Convert to matrix and order by the ARI values
plotData <- as.matrix(plotData)
plotData <- plotData[order(plotData[,1], decreasing = F),]
rownames(plotData) <- paste("ncluster:",plotData[,1], sep = "")
plotData <- plotData[,-1]

# Add mean, SD, and CV values for each model
plotData2 <- rbind(rbind(c(classiMetricMean$ARI_mean, classiMetricMean$NMI_mean),c(classiMetricMean$ARI_sd, classiMetricMean$NMI_sd)),
                   c(classiMetricMean$ARI_CV, classiMetricMean$NMI_CV))
rownames(plotData2) <- c("Mean", "SD", "CV")
plotData <- rbind(plotData, plotData2)

write.table(plotData, "Analysis/Slices_representation/result/MetricRe.csv", col.names = T, row.names = T, sep = ",", quote = F)

### 3. Heatmap visualization ----
# This section generates a heatmap for the ARI and NMI values across different models.
# The heatmap is customized with color schemes and annotations for easy interpretation.
col_group <- data.frame(Model = factor(c(models, models)), Metric =  factor(c(rep("ARI", length(models)), rep("NMI", length(models)))))
rownames(col_group) <- colnames(plotData)

# Define color palettes for models and metrics
modelColor = c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC",  "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "PRECAST" = "#82564c")
metricColor = c("ARI" = "#74A8D2", "NMI" = "#C74B3E")
ann_colors<-list(Model = modelColor, Metric = metricColor)

# Custom color scale for the heatmap
custom_colors <- colorRampPalette(c("#5F83b7", "#b6d1e4", "#fefec7", "#f4c486", "#cf533e"))(100)
breaks <- c(seq(0, 0.4, length.out = 80),
            seq(0.41, 1, length.out = 20))

pdf("Analysis/Slices_representation/result/heatmap.pdf",height = 6,width  = 5)
pheatmap(plotData[-which(rownames(plotData) == "SD" ),], scale = "none",
         show_colnames = F, show_rownames = TRUE, cluster_cols = F, cluster_rows = F,
         annotation_col = col_group, annotation_colors = ann_colors,
         border_color = "black", fontsize = 6, display_numbers = TRUE,
         number_color = "black", number_format = "%.2f",
         gaps_col = 5, gaps_row = nrow(plotData)-3,
         color = custom_colors,
         breaks = breaks)
dev.off()

## kmeans----
### 1. Read data and prepare classification metrics ----
models <- c("Banksy", "CellCharter", "CN", "MENDER", "NicheCompass", "PRECAST")

classiMetricRe <- data.frame()
for(m in models){
  file <- paste("Data/TNBC/SlicesEmbedding/", m, "/Metric/kmeans/model_metric.csv", sep = "")
  re <- read.csv(file)
  colnames(re) <- c("ncluster","ARI", "NMI", "ncluster2")
  re$Model <- m
  classiMetricRe <- rbind(classiMetricRe, re)
}
write.table(classiMetricRe, "Analysis/Slices_representation/result/representation_re_kmeans.csv", col.names = T, row.names = F, quote = F, sep = ",")

classiMetricRe2 <- classiMetricRe
classiMetricRe <- classiMetricRe[,-1]

classiMetricMean <- classiMetricRe %>%
  group_by(Model) %>%
  summarise(
    ARI_mean = mean(ARI),
    ARI_sd = sd(ARI),
    NMI_mean = mean(NMI),
    NMI_sd = sd(NMI)
  )

classiMetricMean$ARI_CV <- classiMetricMean$ARI_sd/classiMetricMean$ARI_mean
classiMetricMean$NMI_CV <- classiMetricMean$NMI_sd/classiMetricMean$NMI_mean

classiMetricMean <- classiMetricMean[match(models, classiMetricMean$Model),]

### 2. Data preparation for plotting ----
plotData <- classiMetricRe %>%
  pivot_wider(
    names_from = Model,
    values_from = c(ARI, NMI),
    names_sep = "_"
  )
plotData <- as.matrix(plotData)
plotData <- plotData[order(plotData[,1], decreasing = F),]
rownames(plotData) <- paste("ncluster:",plotData[,1], sep = "")
plotData <- plotData[,-1]

plotData2 <- rbind(rbind(c(classiMetricMean$ARI_mean, classiMetricMean$NMI_mean),c(classiMetricMean$ARI_sd, classiMetricMean$NMI_sd)),
                   c(classiMetricMean$ARI_CV, classiMetricMean$NMI_CV))
rownames(plotData2) <- c("Mean", "SD", "CV")
plotData <- rbind(plotData, plotData2)

### 3. Heatmap visualization ----
col_group <- data.frame(Model = factor(c(models, models)), Metric =  factor(c(rep("ARI", length(models)), rep("NMI", length(models)))))
rownames(col_group) <- colnames(plotData)

modelColor = c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC",  "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "PRECAST" = "#82564c")
metricColor = c("ARI" = "#74A8D2", "NMI" = "#C74B3E")
ann_colors<-list(Model = modelColor, Metric = metricColor)

custom_colors <- colorRampPalette(c("#5F83b7", "#b6d1e4", "#fefec7", "#f4c486", "#cf533e"))(100)
breaks <- c(seq(0, 0.4, length.out = 80),
            seq(0.41, 1, length.out = 20))

pdf("Analysis/Slices_representation/result/heatmap_kmeans.pdf",height = 6,width  = 5)
pheatmap(plotData[-which(rownames(plotData) == "SD" ),], scale = "none",
         show_colnames = F, show_rownames = TRUE, cluster_cols = F, cluster_rows = F,
         annotation_col = col_group, annotation_colors = ann_colors,
         border_color = "black", fontsize = 6, display_numbers = TRUE,
         number_color = "black", number_format = "%.2f",
         gaps_col = 6, gaps_row = nrow(plotData)-3,
         color = custom_colors,
         breaks = breaks)
dev.off()

## 1. clustering method compaire----
### 1.1 read data----
#### hclust 数据读取----
metric <- read.csv("Analysis/Slices_representation/result/representation_re.csv")
metric <- metric %>%
  rownames_to_column(var = "nclust")  # 将行名转为一列，命名为 "Sample"

metric <- metric %>%
  pivot_longer(
    cols = -nclust,
    names_to = "Metric",
    values_to = "Value"
  )

metric$Model <- sapply(metric$Metric, function(x){
  unlist(strsplit(x, split = "_"))[2]
})

metric$MetricType <- sapply(metric$Metric, function(x){
  unlist(strsplit(x, split = "_"))[1]
})
metric$Type <- "hclust"

#### kmeans数据读取----
metric_kmeans <- read.csv("Analysis/Slices_representation/result/representation_re_kmeans.csv")
metric_kmeans <- metric_kmeans %>%
  rownames_to_column(var = "nclust")  # 将行名转为一列，命名为 "Sample"

metric_kmeans <- metric_kmeans %>%
  pivot_longer(
    cols = -nclust,
    names_to = "Metric",
    values_to = "Value"
  )

metric_kmeans$Model <- sapply(metric_kmeans$Metric, function(x){
  unlist(strsplit(x, split = "_"))[2]
})

metric_kmeans$MetricType <- sapply(metric_kmeans$Metric, function(x){
  unlist(strsplit(x, split = "_"))[1]
})
metric_kmeans$Type <- "Kmeans"

### 1.2.ARI box plot----
plotData <- rbind(metric[metric$MetricType == "ARI", ], metric_kmeans[metric_kmeans$MetricType == "ARI",])
plotData <- plotData[plotData$nclust != "Mean" &
                       plotData$nclust != "SD" & 
                       plotData$nclust != "CV",]

plotData <- as.data.frame(plotData)
plotData <- plotData[,match(c("nclust", "Model", "Value", "Type"), colnames(plotData))]
plotData$Type <- factor(plotData$Type, levels = c("hclust", "Kmeans"))

p <- ggplot(plotData, aes(x = Model, y = Value, fill = Type)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.7), 
               alpha = 0.7, size = 0.5, width = 0.5) +  
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", 
               position = position_dodge(0.7), width = 0.2, size = 0.5) +
  geom_jitter(aes(color = Type), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), 
              size = 0.4) +
  stat_compare_means(aes(group = Type), method = "t.test", 
                     label = "p.signif", size = 2.5) +
  scale_fill_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  scale_color_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  theme_classic() +  
  labs(x = NULL, y = "ARI") +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3)) +
  NoLegend() 
ggsave("Analysis/Slices_representation/result/ARI.pdf", p, width = 80, height = 50, units = "mm")


###3.NMI box plot----
plotData <- rbind(metric[metric$MetricType == "NMI", ],  metric_kmeans[metric_kmeans$MetricType == "NMI",])
plotData <- plotData[plotData$nclust != "Mean" &
                       plotData$nclust != "SD" & 
                       plotData$nclust != "CV",]

plotData <- as.data.frame(plotData)
plotData <- plotData[,match(c("nclust", "Model", "Value", "Type"), colnames(plotData))]
plotData$Type <- factor(plotData$Type, levels = c("hclust", "Kmeans"))

p <- ggplot(plotData, aes(x = Model, y = Value, fill = Type)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.7), 
               alpha = 0.7, size = 0.5, width = 0.5) +  
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", 
               position = position_dodge(0.7), width = 0.2, size = 0.5) +
  geom_jitter(aes(color = Type), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), 
              size = 0.4) +
  stat_compare_means(aes(group = Type), method = "t.test", 
                     label = "p.signif", size = 2.5) +
  scale_fill_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  scale_color_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  theme_classic() +  
  labs(x = NULL, y = "NMI") +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3)) +
  NoLegend() 
ggsave("Analysis/Slices_representation/result/NMI.pdf", p, width = 80, height = 50, units = "mm")


###4.ARI CV plot----
plotData <- rbind(metric[metric$MetricType == "ARI", ],  metric_kmeans[metric_kmeans$MetricType == "ARI",])
plotData <- plotData[plotData$nclust == "CV",]
plotData$Type <- factor(plotData$Type, levels = c("hclust", "Kmeans"))

p <- ggplot(plotData, aes(x = Type, y = Value, fill = Type)) +
  stat_summary(fun = mean, geom = "bar", alpha = 0.7, size = 0.5, width = 0.5) +  
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, size = 0.5) +
  geom_jitter(aes(color = Type), position = position_jitterdodge(jitter.width = 0.15), size = 0.4) +
  scale_fill_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  scale_color_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  theme_classic() +  
  labs(x = NULL, y = "ARI CV") +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3)) +
  NoLegend() 
ggsave("Analysis/Slices_representation/result/ARI_CV_barplot.pdf", p, width = 30, height = 50, units = "mm")


###5.NMI CV plot----
plotData <- rbind(metric[metric$MetricType == "NMI", ], metric_kmeans[metric_kmeans$MetricType == "NMI",])
plotData <- plotData[plotData$nclust == "CV",]
plotData$Type <- factor(plotData$Type, levels = c("hclust", "Kmeans"))

p <- ggplot(plotData, aes(x = Type, y = Value, fill = Type)) +
  stat_summary(fun = mean, geom = "bar", alpha = 0.7, size = 0.5, width = 0.5) +  
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, size = 0.5) +
  geom_jitter(aes(color = Type), position = position_jitterdodge(jitter.width = 0.15), size = 0.4) +
  scale_fill_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  scale_color_manual(values = c("#5d80b5","#509C3D", "#e66231")) +
  theme_classic() +  
  labs(x = NULL, y = "NMI CV") +  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3)) +
  NoLegend() 
ggsave("Analysis/Slices_representation/result/NMI_CV_barplot.pdf", p, width = 30, height = 50, units = "mm")



## 2.domain-celltype analysis----
domain_plot <- function(file, ex, saveFile){
  
  ### Step 1: Load input data ----
  # Read domain abundance matrix (output from domain decomposition)
  abundance_matrix <- read.csv(file, check.names = F)
  
  # Read spatial transcriptomics dataset (AnnData format)
  data <- read_h5ad(ex)
  metadata <- data$obs  # Extract metadata including predicted domains and cell types
  
  ### Step 2: Plot domain abundance per slice ----
  # Define color palette for domains
  color <- c("#509C3D", "#C33931", "#58BACC", "#3A73AE", "#EE8435", "#8D68B8", "#BABD45",
             "#a7dc92", "#D57BBE", "#b3c6e3", "#f19c97", "#82564c", "#FFCC00", "#6D9EC1",
             "#E6739F", "#76C16E", "#A469FF")
  
  # Reshape the abundance matrix to long format for ggplot2
  abundance_long <- reshape2::melt(abundance_matrix,
                                   id.vars = c("slices", "slices_class"),
                                   variable.name = "domain",
                                   value.name = "abundance")
  
  # Ensure 'slices_class' is a factor with defined order
  abundance_long$slices_class <- factor(abundance_long$slices_class, 
                                        levels = c("Cold", "Mixed", "Compartmentalized"))
  
  # Create stacked bar plot showing domain abundance per slice, grouped by slice class
  p1 <- ggplot(abundance_long, aes(x = factor(slices), y = abundance, fill = domain)) +
    geom_bar(stat = "identity", width = 0.8) + 
    facet_grid(. ~ slices_class, scales = "free_x", space = "free_x") +  # Facet by slice class
    labs(x = "Sample (Slice)", y = "Abundance", fill = "Domain") +
    scale_fill_manual(values = color) +
    theme_classic() + 
    theme(
      legend.key.size = unit(0.1, 'cm'),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
      axis.text.x = element_text(angle = 45, size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      axis.title.y = element_text(size = 6, color = "black"),
      axis.title.x = element_text(size = 6, color = "black"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.2)
    )
  
  ### Step 3: Plot cell type composition per domain ----
  # Summarize cell type proportions within each predicted domain
  plot_data <- metadata %>%
    group_by(predicted_domain, celltype) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(predicted_domain) %>%
    mutate(proportion = count / sum(count))  # Calculate cell type proportions
  
  # Set factor levels for consistent ordering in plots
  plot_data$predicted_domain <- factor(plot_data$predicted_domain)
  plot_data$celltype <- factor(plot_data$celltype, 
                               levels = c("B", "CD3 T", "CD4 T", "CD8 T", "DC", "DC/Mono", 
                                          "Endothelial", "Keratin+ tumor", "Macrophages",
                                          "Mesenchymal like", "Mono/Neu", "NK", 
                                          "Neutrophils", "Other immune", "Tregs", 
                                          "Tumor", "Unidentified"))
  
  # Define color palette for cell types
  color2 <- c("#509C3D", "#3A73AE", "#58BACC", "#b3c6e3", "#a7dc92", "#BABD45", "#FFCC00", "#E6739F",
              "#76C16E", "#EE8435", "#82564c", "#6D9EC1", "#8D68B8", "#69B4DF", "#25C18E", "#C33931", "grey70")
  
  # Create stacked bar plot showing cell type composition within each domain
  p2 <- ggplot(plot_data, aes(x = predicted_domain, y = proportion, fill = celltype)) +
    geom_bar(stat = "identity", width = 0.8) +
    labs(x = "Domain", y = "Celltype proportion", fill = "Celltype") +
    scale_fill_manual(values = color2) +
    theme_classic() + 
    theme(
      legend.key.size = unit(0.1, 'cm'),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
      axis.text.x = element_text(size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      axis.title.y = element_text(size = 6, color = "black"),
      axis.title.x = element_text(size = 6, color = "black"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.2)
    )
  
  ### Step 4: Combine and save plots ----
  # Combine domain abundance plot (p1) and cell type composition plot (p2) side-by-side
  p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(2, 1))  # p1 is wider than p2
  
  # Save the combined figure to file
  ggsave(saveFile, p, width = 15, height = 6, units = "cm")
}

###MENDER 5
file <- "Data/TNBC/SlicesEmbedding/MENDER/5/abundance_matrix_normal.csv"
ex <- "Data/TNBC/SlicesEmbedding/MENDER/5/MENDER.h5ad"
saveFile <- "Data/TNBC/SlicesEmbedding/MENDER/5/MENDER.pdf"
domain_plot(file, ex, saveFile)

