library(anndata)
library(ggplot2)
library(patchwork)
library(plotly)
library(cowplot)
library(reshape2)
library(scatterplot3d)
library(Seurat)
setwd("iSTBench")

### 1. Read Data ----
# Load metadata and prepare necessary variables
metadata <- read.csv("Data/BaristaSeq/IntergrationRe/Metric/metaPredictedRe.csv")
domain = unique(metadata$original_domain)
model <- colnames(metadata)[6:ncol(metadata)]
slices <- sort(unique(metadata$slices), decreasing = T)

# Load spatial coordinate data from the rotated files
fl <- list.files("Data/BaristaSeq/sample_data_rotate")
rotateSpatial <- data.frame()
for(f in fl){
  f_data <- read_h5ad(paste("Data/BaristaSeq/sample_data_rotate", f, sep = "/"))
  f_metadata <- f_data$obs
  f_spatial <- f_data$obsm$spatial
  f_rotate <- cbind(f_metadata$barcode, f_spatial)
  rotateSpatial <- rbind(rotateSpatial, f_rotate)
}
colnames(rotateSpatial) <- c("barcode", "X", "Y")
rotateSpatial <- rotateSpatial[match(metadata$barcode, rotateSpatial$barcode) ,]

# Add spatial coordinates to metadata
metadata$X <- rotateSpatial$X
metadata$Y <- rotateSpatial$Y


### 2. Domain Match ----
# Loop through each model and calculate domain matching for plotting
for(m in model){
  model_doamin <- metadata[,c("original_domain", m)]
  model_domain_count <- as.data.frame(table(model_doamin[,1], model_doamin[,2]))
  pre_domain <- names(sort(table(model_doamin[,2]),decreasing = T))
  
  pd2domain <- data.frame()
  for(pd in pre_domain){
    a <- model_domain_count[model_domain_count[,2] == pd,][order(model_domain_count[model_domain_count[,2] == pd,]$Freq,decreasing = T)[1],1] 
    pd2domain <- rbind(pd2domain, c(as.character(pd), as.character(a)))
    model_domain_count <- model_domain_count[-which(model_domain_count[,1] == a),]
  }
  colnames(pd2domain) <- c(m, paste("New",m,sep = "_"))
  
  metadata <- merge(metadata, pd2domain, by = m)
  
}

# Clean up metadata by removing old model columns
metadata <- metadata[, -match(model, colnames(metadata))]
model <- sapply(colnames(metadata)[6:ncol(metadata)], function(x){unlist(strsplit(x, split = "_"))[2] })
colnames(metadata) <- c(colnames(metadata)[1:4],"Ground_truth", model)
metadata$Ground_truth <- as.character(metadata$Ground_truth)

### 3. Add Corrected Spatial ----
# Loop through each model and add corrected spatial data
for(m in model){
  m_data <- read_h5ad(paste("Benchmark/Alignment/Result/BaristaSeq", m, "Spatial_correct_data.h5ad", sep = "/"))
  m_metadata <- m_data$obs
  index <- match(metadata$barcode, m_metadata$barcode )
  
  spatial <- m_data$obsm$spatial[index,]
  colnames(spatial) <- paste(m, c("X", "Y"), sep = "_")
  metadata <- cbind(metadata, spatial)
}


### 4. Plotting ----
# 4.1: 2D plot of color-coded slices ----
# Define colors for slices
slices
color <- c("1" = "#FFFF54", "2" = "#6CE3FB", "3" = "#EB49F7" )
plot_list <- list()

# Prepare data for plotting
plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", "X", "Y"), colnames(metadata)))]
plot_data$X <- as.numeric(plot_data$X)
plot_data$Y <- as.numeric(plot_data$Y)
plot_data$slices <- factor(as.character(plot_data$slices), levels = c("1", "2", "3"))

# Generate plot for original data
p <- ggplot(plot_data, aes(x = X, y = Y, col = slices )) +
  geom_point(size = 0.2, alpha= 0.75) +
  scale_color_manual(values = color) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        plot.title = element_blank()) +
  coord_fixed(ratio = 1)
plot_list[["original"]] <- p

# Generate plots for each model
for(m in model){
  plot_data <- metadata[ ,c(match(c("barcode", "slices", m, paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
  colnames(plot_data)[3:5] <- c("model", "X", "Y")
  plot_data$X <- as.numeric(plot_data$X)
  plot_data$Y <- as.numeric(plot_data$Y)
  plot_data$slices <- factor(as.character(plot_data$slices), levels = c("1", "2", "3"))
  
  
  p <- ggplot(plot_data, aes(x = X, y = Y, col = slices )) +
    geom_point(size = 0.2, alpha= 0.75) +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          plot.title = element_blank()) +
    coord_fixed(ratio = 1)
  
  plot_list[[m]] <- p
  
}

# Combine and save all 2D plots
p_all <- wrap_plots(plot_list, ncol = 9, heights = 1, widths = 1)
ggsave("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_2D_plot_slices.pdf", plot = p_all, width = 20, height = 6.75, units = "in")


# 4.2: 2D plot color-coded by domain ----
slices
color <- c("VISp_I" = "#2F6DA1", "VISp_II/III" = "#EB49F7", "VISp_IV" = "#6CE3FB", "VISp_V" = "#EB594F", "VISp_VI" = "#FFFF54", "VISp_wm" ="#3B8749" )

plot_list <- list()

plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", "X", "Y"), colnames(metadata)))]
plot_data$X <- as.numeric(plot_data$X)
plot_data$Y <- as.numeric(plot_data$Y)

p <- ggplot(plot_data, aes(x = X, y = Y, col = Ground_truth )) +
  geom_point(size = 0.2, alpha= 0.75) +
  scale_color_manual(values = color) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        plot.title = element_blank()) +
  coord_fixed(ratio = 1)
plot_list[["original"]] <- p

for(m in model){
  plot_data <- metadata[ ,c(match(c("barcode", "slices", m, paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
  colnames(plot_data)[3:5] <- c("Ground_truth", "X", "Y")
  plot_data$X <- as.numeric(plot_data$X)
  plot_data$Y <- as.numeric(plot_data$Y)
  
  p <- ggplot(plot_data, aes(x = X, y = Y, col = Ground_truth )) +
    geom_point(size = 0.2, alpha= 0.75) +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          plot.title = element_blank()) +
    coord_fixed(ratio = 1)
  
  plot_list[[m]] <- p
  
}

p_all <- wrap_plots(plot_list, ncol = 9, heights = 1, widths = 1)
ggsave("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_2D_plot_domain.pdf", plot = p_all, width = 20, height = 6.75, units = "in")

# 4.3: 3D plot color-coded by domain ----
color <- c("VISp_I" = "#2F6DA1", "VISp_II/III" = "#EB49F7", "VISp_IV" = "#6CE3FB", "VISp_V" = "#EB594F", "VISp_VI" = "#FFFF54", "VISp_wm" ="#3B8749" )

plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", "X", "Y"), colnames(metadata)))]
Z <- rep(1, nrow(plot_data))
Z[which(plot_data$slices == 2)] <-11
Z[which(plot_data$slices == 1)] <- 21
plot_data$Z <- Z
colnames(plot_data)[3] <- c("model")
plot_data$color <- color[plot_data$model]
plot_data$X <- as.numeric(plot_data$X)
plot_data$Y <- as.numeric(plot_data$Y)

pdf("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_3D_plot.pdf", width = 60, height = 6.75)
par(mfrow = c(1, (length(model)+1)))
p <- scatterplot3d(plot_data$X, plot_data$Y, plot_data$Z, color = plot_data$color, pch = 20, 
                   xlab = "", ylab = "", zlab = "", 
                   ylim = c(min(plot_data$Y), max(plot_data$Y)), xlim = c(min(plot_data$X), max(plot_data$X)), zlim = c(min(plot_data$Z), max(plot_data$Z)),
                   grid = F, box = F, cex.symbols = 0.5, tick.marks = F,angle = 45, axis = F)


for(m in model){
  plot_data <- metadata[ ,c(match(c("barcode", "slices", m, paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
  Z <- rep(1, nrow(plot_data))
  Z[which(plot_data$slices == 2)] <- 11
  Z[which(plot_data$slices == 1)] <- 21
  plot_data$Z <- Z
  colnames(plot_data)[3:5] <- c("model", "X", "Y")
  plot_data$color <- color[plot_data$model]
  plot_data$X <- as.numeric(plot_data$X)
  plot_data$Y <- as.numeric(plot_data$Y)
  
  scatterplot3d(plot_data$X, plot_data$Y, plot_data$Z, color = plot_data$color, pch = 20, 
                xlab = "", ylab = "", zlab = "", 
                ylim = c(min(plot_data$Y), max(plot_data$Y)), xlim = c(min(plot_data$X), max(plot_data$X)), zlim = c(min(plot_data$Z), max(plot_data$Z)),
                grid = F, box = F, cex.symbols = 0.5, tick.marks = F, angle = 45, axis = F)
  
}
dev.off()
###

### 5. Metric plot ----
accuracy <- read.csv("Benchmark/Alignment/Result/BaristaSeq/Metric/Accuracy.csv")
ratio <- read.csv("Benchmark/Alignment/Result/BaristaSeq/Metric/Ratio.csv")

model <- c("Banksy", "CellCharter", "CN", "GraphST", "MENDER", "NicheCompass", "Spado")

accuracy <- accuracy[match(model, accuracy$model),]
ratio <- ratio[match(model, ratio$model),]

accuracy_long <- melt(accuracy, id.vars = "model", variable.name = "slices", value.name = "value")
ratio_long <- melt(ratio, id.vars = "model", variable.name = "slices", value.name = "value")
color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "Spado" = "#D57BBE" )

p1 <- ggplot(accuracy_long, aes(x = model, y = value, fill = model)) +
  geom_boxplot(lwd = 0.3, alpha = 1,outlier.size = 0.3, outlier.shape = 16) +
  geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
              size = 0.2) +
  scale_color_manual(values = color)+
  labs(x = NULL, y = "Accuracy") +
  scale_fill_manual(values = color) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3))+
  NoLegend()

p2 <- ggplot(ratio_long, aes(x = model, y = value, fill = model)) +
  geom_boxplot(lwd = 0.3, alpha = 1,outlier.size = 0.3, outlier.shape = 16) +  # 去除离群点
  geom_jitter(aes(color = model), position = position_jitterdodge(jitter.width = 0.15), 
              size = 0.2) +
  scale_color_manual(values = color)+
  labs(x = NULL, y = "Ratio") +
  scale_fill_manual(values = color) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3))+
  NoLegend()

p <- p1 + p2
ggsave("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_metric.pdf", plot = p, width = 10, height = 5, units = "cm")

