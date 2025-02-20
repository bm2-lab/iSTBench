# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(scales) 

# Set working directory for storing results
setwd("/NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub_test/iSTBench")

### 1. Read data and prepare classification metrics ----
# This section reads metric data for different models and combines them into a single data frame.
# It calculates ARI (Adjusted Rand Index) and NMI (Normalized Mutual Information) for each model.
models <- c("Banksy", "CellCharter", "CN", "MENDER", "NicheCompass")

# Initialize an empty dataframe to store classification metrics for all models
classiMetricRe <- data.frame()
for(m in models){
  file <- paste("~/home_dkj/FD_yzy/Dataset/Slices_embedding/6_TNBC/SlicesEmbedding/", m, "/Metric/hclust/model_metric.csv", sep = "")
  # file <- paste("Data/TNBC/SlicesEmbedding/", m, "/Metric/hclust/model_metric.csv", sep = "")
  
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
modelColor = c("Banksy" = "#C74B3E", "CellCharter" = "#539B4F", "CN" = "#5E93BF",  "MENDER" = "#BAD691", "NicheCompass" ="#BCD3E6")
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

### 4. PCA plot function ----
# This function generates PCA plots for a given dataset, showing the distribution of slices based on PCA components.
# It also generates boxplots for the first and second principal components to show differences across classes.
PCA_plot_fun <- function(file, ex, saveFile, saveFile2){
  # Read input data: abundance matrix and metadata
  data <- read_h5ad(ex)
  metadata<-data$obs
  real <- unique(data.frame(slices = metadata$slices, 
                            real_class = metadata$slices_class))
  real$slices <- as.integer(as.character(real$slices))
  
  # Read and preprocess abundance matrix
  abundance_matrix <- read.csv(file)
  rownames(abundance_matrix) <- abundance_matrix$slices
  abundance_matrix <- abundance_matrix[,-match(c("slices_class", "slices"), colnames(abundance_matrix))]
  
  # Perform PCA on the abundance matrix
  set.seed(1234)
  pca_result <- prcomp(abundance_matrix, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x)
  pca_df$slices <- rownames(pca_df)
  
  # Merge PCA results with real class information
  pca_df <- merge(pca_df, real, by = "slices")
  pca_df$real_class <- as.character(pca_df$real_class)
  if(length(which(pca_df$real_class == 0)) != 0){
    pca_df[which(pca_df$real_class == 0),]$real_class <- "Cold"
    pca_df[which(pca_df$real_class == 1),]$real_class <- "Compartmentalized"
    pca_df[which(pca_df$real_class == 2),]$real_class <- "Mixed"
  }
  
  # Set factor levels for class labels and define custom colors for plotting
  pca_df$real_class <- factor(pca_df$real_class, levels = c("Cold", "Mixed", "Compartmentalized"))
  custom_colors <- c("Cold" = "#93C356", "Mixed" = "#69B4DF", "Compartmentalized" = "#D44D40") 
  
  # Create the main PCA plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = real_class)) +
    geom_point(size = 0.5) +
    labs(title = NULL, x = "PC1", y = "PC2") +
    scale_color_manual(values = custom_colors) +
    scale_x_continuous(labels = number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    #scale_y_continuous(labels = scientific_format(accuracy = 0.1)) +
    theme_classic() + 
    theme(legend.key.size = unit(0.1, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          axis.text.x = element_text(size = 6, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 6, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.2))+
    NoLegend()
  
  # Add marginal boxplots for PC1 and PC2
  p <- ggMarginal(p, type = "boxplot", margins = "both", size = 6, colour = c("#93C356", "#69B4DF", "#D44D40"), groupFill = T,
                  xparams = list(size = 0.3, outlier.shape = NA), yparams = list(size = 0.3, outlier.shape = NA))
  ggsave(saveFile, p, width = 5, height = 5, units = "cm")
  
  # Create boxplots for PC1 and PC2 comparisons across real classes
  comparisons <- list(c("Cold", "Mixed"), c("Cold", "Compartmentalized"), c("Mixed", "Compartmentalized"))
  p_box_PC1 <-  ggplot(pca_df, aes(x = real_class, y = PC1, fill = real_class)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif", size = 2) +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    scale_fill_manual(values = custom_colors) + 
    theme_classic() +
    theme(axis.text.x = element_text(size = 0, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 6, color = "black"))+
    NoLegend()
  
  p_box_PC2 <-  ggplot(pca_df, aes(x = real_class, y = PC2, fill = real_class)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", size = 2) +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    scale_fill_manual(values = custom_colors) + 
    theme_classic() +
    theme(axis.text.x = element_text(size = 0, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 6, color = "black"))+
    NoLegend()
  
  p2 <- grid.arrange(p_box_PC1, p_box_PC2, ncol = 2)
  ggsave(saveFile2, p2, width = 4, height = 5, units = "cm")
  
}

### 5. Run PCA analysis for CN model (example) ----
# Example usage of the PCA plotting function for the "CN" model, slice 46
file <- "Data/TNBC/SlicesEmbedding/CN/46/abundance_matrix_normal.csv"
ex <- "Data/TNBC/SlicesEmbedding/CN/46/CN.h5ad"
saveFile <- "Analysis/Slices_representation/result/CN_46.pdf"
saveFile2 <- "Analysis/Slices_representation/result/CN_46_box.pdf"
PCA_plot_fun(file, ex, saveFile, saveFile2)

