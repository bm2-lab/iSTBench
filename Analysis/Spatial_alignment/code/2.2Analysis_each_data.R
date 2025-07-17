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
model <- model[-which(model == "GraphSTwithPASTE")]
model <- c("PASTE", "STalign", model, paste("SPACEL", model, sep = "_"))
for(m in model){
  m_data <- read_h5ad(paste("/home/dongkj/home_dkj/FD_yzy/Landmark_based_intergration/Dataset/BaristaSeq", m, "Spatial_correct_data.h5ad", sep = "/"))
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
ggsave("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_2D_plot_slices_original.pdf", plot = p, width = 20/9, height = 6.75, units = "in")


# Generate plots for each model
for(m in model){
  if(m %in% c("PASTE", "STalign")){
    plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
    title <- m
  }else if(m %in% c("Banksy", "CellCharter", "CN", "GraphST", "MENDER","NicheCompass", "PRECAST", "Spado", "SPIRAL", "STAIG", "STAligner")){
    plot_data <- metadata[ ,c(match(c("barcode", "slices", m, paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
    title <- paste("STAligner", m, sep = "-")
  }else{
    m1 <- unlist(strsplit(m, split = "_"))[2]
    plot_data <- metadata[ ,c(match(c("barcode", "slices", m1, paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
    title <- paste("SPACEL", m1, sep = "-")
  }
  
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
          plot.title = element_text(hjust = 0.5, size = 6)) +
    coord_fixed(ratio = 1)+
    ggtitle(title) 
  
  plot_list[[m]] <- p
  
}

# Combine and save all 2D plots
p_all <- wrap_plots(plot_list, ncol = 8, heights = 1, widths = 1)
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
ggsave("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_2D_plot_domain_original.pdf", plot = p, width = 20/9, height = 6.75, units = "in")

for(m in model){
  if(m %in% c("PASTE", "STalign")){
    plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
    title <- m
  }else if(m %in% c("Banksy", "CellCharter", "CN", "GraphST", "MENDER","NicheCompass", "PRECAST", "Spado", "SPIRAL", "STAIG", "STAligner")){
    plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
    title <- paste("STAligner", m, sep = "-")
  }else{
    m1 <- unlist(strsplit(m, split = "_"))[2]
    plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
    title <- paste("SPACEL", m1, sep = "-")
  }
  
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
          plot.title = element_text(hjust = 0.5, size = 6)) +
    coord_fixed(ratio = 1)+
    ggtitle(title) 
  
  plot_list[[m]] <- p
  
}

p_all <- wrap_plots(plot_list, ncol = 8, heights = 1, widths = 1)
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

pdf("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_3D_plot_original.pdf", width = 60/9, height = 6.75)
p <- scatterplot3d(plot_data$X, plot_data$Y, plot_data$Z, color = plot_data$color, pch = 20, 
                   xlab = "", ylab = "", zlab = "", 
                   ylim = c(min(plot_data$Y), max(plot_data$Y)), xlim = c(min(plot_data$X), max(plot_data$X)), zlim = c(min(plot_data$Z), max(plot_data$Z)),
                   grid = F, box = F, cex.symbols = 0.5, tick.marks = F,angle = 45, axis = F)
dev.off()


pdf("Benchmark/Alignment/Result/BaristaSeq/Metric/BaristaSeq_3D_plot.pdf", width = 60/9*8, height = 6.75)
par(mfrow = c(3,8))
for(m in model){
  if(m %in% c("PASTE", "STalign")){
    plot_data <- metadata[ ,c(match(c("barcode", "slices", "Ground_truth", paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
  }else if(m %in% c("Banksy", "CellCharter", "CN", "GraphST", "MENDER","NicheCompass", "PRECAST", "Spado", "SPIRAL", "STAIG", "STAligner")){
    plot_data <- metadata[ ,c(match(c("barcode", "slices", m, paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
  }else{
    m1 <- unlist(strsplit(m, split = "_"))[2]
    plot_data <- metadata[ ,c(match(c("barcode", "slices", m1, paste(m,"X", sep = "_"), paste(m,"Y", sep = "_") ), colnames(metadata)))]
  }
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
accuracy <- accuracy[match(model, accuracy$model),]
ratio <- ratio[match(model, ratio$model),]

accuracy$model <- sapply(accuracy$model, function(x){
  if(x %in% c("Banksy", "CellCharter", "CN", "GraphST", "MENDER", "NicheCompass", "PRECAST", "Spado", "SPIRAL", "STAIG", "STAligner")){
    x <- paste("STAligner", x, sep = "_")
  }
  return(x)
})

ratio$model <- sapply(ratio$model, function(x){
  if(x %in% c("Banksy", "CellCharter", "CN", "GraphST", "MENDER", "NicheCompass", "PRECAST", "Spado", "SPIRAL", "STAIG", "STAligner")){
    x <- paste("STAligner", x, sep = "_")
  }
  return(x)
})

model <- sapply(model, function(x){
  if(x %in% c("Banksy", "CellCharter", "CN", "GraphST", "MENDER", "NicheCompass", "PRECAST", "Spado", "SPIRAL", "STAIG", "STAligner")){
    x <- paste("STAligner", x, sep = "_")
  }
  return(x)
})

accuracy_long <- melt(accuracy, id.vars = "model", variable.name = "slices", value.name = "value")
ratio_long <- melt(ratio, id.vars = "model", variable.name = "slices", value.name = "value")
accuracy_long$model <- factor(accuracy_long$model, levels = model)
ratio_long$model <- factor(ratio_long$model, levels = model)

color <- c("PASTE" = "#6CE3FB", "STalign" = "#EB49F7", 
           "STAligner_Banksy" = "#509C3D", "STAligner_CellCharter" = "#C33931", "STAligner_CN" = "#58BACC", "STAligner_GraphST" = "#3A73AE", "STAligner_GraphST-PASTE" = "#EE8435", "STAligner_MENDER" = "#8D68B8",
           "STAligner_NicheCompass" ="#BABD45", "STAligner_PRECAST" = "#82564c", "STAligner_Spado" = "#D57BBE", "STAligner_SPIRAL" = "#b3c6e3", "STAligner_STAIG" = "#f19c97", "STAligner_STAligner" = "#a7dc92",
           "SPACEL_Banksy" = "#78B963", "SPACEL_CellCharter" = "#D65A52", "SPACEL_CN" = "#81CCE0", "SPACEL_GraphST" = "#5F8FC3","SPACEL_GraphST-PASTE" = "#F3A765","SPACEL_MENDER" = "#A88ACC",
           "SPACEL_NicheCompass" = "#D0D358","SPACEL_PRECAST" = "#A0766A", "SPACEL_Spado" = "#E39ED2","SPACEL_SPIRAL" = "#c7d6eb","SPACEL_STAIG" = "#f5b3af","SPACEL_STAligner" = "#bde7a8" )



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
