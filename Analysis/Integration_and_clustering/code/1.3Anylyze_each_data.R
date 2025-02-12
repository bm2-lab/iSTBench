#Take BaristaSeq as an example
library(ggplot2)
library(patchwork)
library(pdist)
require(parallel)
library(ggbreak)
library(reshape)
library(Seurat)
library(gridExtra)
library(anndata)
library(Seurat)

fx_1NN = function(i,location_in){
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  return(min(line_i))
}

fx_kNN = function(i,location_in,k,cluster_in){
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  ind = order(line_i)[1:k]
  cluster_use = cluster_in[-i]
  if(sum(cluster_use[ind] != cluster_in[i])>(k/2)){
    return(1)
  }else{
    return(0)
  }
  
}

fx_CHAOS = function(clusterlabel, location){
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  matched_location = scale(matched_location)
  dist_val = rep(0,length(unique(clusterlabel)))
  count = 0
  for(k in unique(clusterlabel)){
    count = count + 1
    location_cluster = matched_location[which(clusterlabel == k),]
    if(length(location_cluster)==2){next}
    results = mclapply(1:dim(location_cluster)[1], fx_1NN, location_in=location_cluster,mc.cores = 5)
    dist_val[count] = sum(unlist(results))
  }
  dist_val = na.omit(dist_val)
  return(sum(dist_val)/length(clusterlabel))
  
}

fx_PAS = function(clusterlabel, location){
  
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  
  results = mclapply(1:dim(matched_location)[1], fx_kNN, location_in=matched_location,k=10,cluster_in=clusterlabel, mc.cores = 5)
  return(sum(unlist(results))/length(clusterlabel))
}

metric1_cal <- function(metadata, m, domain){
  match_domain <- names(sort(table(metadata[metadata$Ground_truth == domain, which(colnames(metadata) == m)]), decreasing = T))[1]
  real <- metadata$Ground_truth
  real[which(real != domain)] <- "No"
  pre <- metadata[,which(colnames(metadata) == m)]
  pre[which(pre != match_domain)] <- "No"
  pre[which(pre == match_domain)] <- domain
  
  real <- factor(real, levels = c("No", domain))
  pre <- factor(pre, levels = c("No", domain))
  
  confusion_matrix <- table(real, pre)
  
  TP <- confusion_matrix[2, 2]
  TN <- confusion_matrix[1, 1]
  FP <- confusion_matrix[1, 2]
  FN <- confusion_matrix[2, 1] 
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  return(c(sensitivity, specificity))
  ###
}
metric2_cal <- function(metadata, m, domain){
  
  match_domain <- names(sort(table(metadata[metadata$Ground_truth == domain,which(colnames(metadata) == m)]), decreasing = T))[1]
  slices <- unique(metadata$slices)
  
  Slices_Sen <- c()
  Slices_Spe <- c()
  for(s in slices){
    s_metadata <- metadata[metadata$slices == s,]
    
    real <- s_metadata$Ground_truth
    real[which(real != domain)] <- "No"
    pre <- s_metadata[,which(colnames(s_metadata) == m)]
    pre[which(pre != match_domain)] <- "No"
    pre[which(pre == match_domain)] <- domain
    
    real <- factor(real, levels = c("No", domain))
    pre <- factor(pre, levels = c("No", domain))
    
    confusion_matrix <- table(real, pre)
    
    TP <- confusion_matrix[2, 2]
    TN <- confusion_matrix[1, 1]
    FP <- confusion_matrix[1, 2]
    FN <- confusion_matrix[2, 1] 
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    
    Slices_Sen <- c(Slices_Sen, sensitivity)
    Slices_Spe <- c(Slices_Spe, specificity)
  }
  
  Slices_Sen <- as.data.frame(t(as.matrix(Slices_Sen)))
  colnames(Slices_Sen) <- slices
  
  Slices_Spe <- as.data.frame(t(as.matrix(Slices_Spe)))
  colnames(Slices_Spe) <- slices
  return(list(Slices_Spe = Slices_Spe,  Slices_Sen = Slices_Sen))
  
}

###set path to iSTBench
setwd("/iSTBench")

metadata <- read.csv("Data/BaristaSeq/IntergrationRe/Metric/metaPredictedRe.csv")
metadata_original <- read.csv("Data/BaristaSeq/IntergrationRe/Metric/metaPredictedRe.csv")

###1. domain matching----
domain = unique(metadata$original_domain)
model <- colnames(metadata)[6:ncol(metadata)]
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
metadata <- metadata[, -match(model, colnames(metadata))]
model <- sapply(colnames(metadata)[6:ncol(metadata)], function(x){unlist(strsplit(x, split = "_"))[2] })
colnames(metadata) <- c(colnames(metadata)[1:4],"Ground_truth", model)
colnames(metadata_original) <- c(colnames(metadata)[1:4],"Ground_truth", model)

###2.clustering metric calculation----
model_ari <- c()
model_nmi <- c()
for(m in model){
  o_domain <- metadata_original$Ground_truth
  p_domain <- metadata_original[,match(m,colnames(metadata_original))]
  
  ari <- aricode::ARI(as.character(o_domain), as.character(p_domain))
  nmi <- aricode::NMI(as.character(o_domain), as.character(p_domain))
  model_ari <- c(model_ari, ari)
  model_nmi <- c(model_nmi, nmi)
}

model_metric <- cbind(model,cbind(model_ari,model_nmi))
model_metric <- as.data.frame(model_metric)
model_metric$model_ari <- as.numeric(model_metric$model_ari)
model_metric$model_nmi <- as.numeric(model_metric$model_nmi)

slices <- unique(metadata_original$slices)
slices <- sort(slices)

CHAOS <- matrix(,length(model),)
PAS <- matrix(,length(model),)
for(s in slices){
  slices_metadata <- metadata_original[metadata_original$slices == s,]
  coor <- cbind(slices_metadata$X, slices_metadata$Y)
  
  model_slices_chaos <- c()
  model_slices_pas <- c()
  for(m in model){
    p_domain <- slices_metadata[,match(m,colnames(slices_metadata))]
    model_slices_chaos <- c(model_slices_chaos, fx_CHAOS(p_domain, coor))
    model_slices_pas <- c(model_slices_pas, fx_PAS(p_domain, coor))
  }
  
  CHAOS <- cbind(CHAOS, model_slices_chaos)
  PAS <- cbind(PAS, model_slices_pas)
}

CHAOS <- as.data.frame(CHAOS)
PAS <- as.data.frame(PAS)

CHAOS <- CHAOS[,-1]
PAS <- PAS[,-1]

CHAOS$Model <- model
rownames(CHAOS) <- model
colnames(CHAOS) <- c(slices, "Model")

PAS$Model <- model
rownames(PAS) <- model
colnames(PAS) <- c(slices, "Model")

###The Settings are based on the actual running model
CHAOS$Model <- factor(CHAOS$Model, levels = c("Banksy", "CellCharter", "CN", "GraphST", "GraphSTwithPASTE", "MENDER", "NicheCompass", "Spado"))
PAS$Model <- factor(PAS$Model, levels = c("Banksy", "CellCharter", "CN", "GraphST", "GraphSTwithPASTE", "MENDER", "NicheCompass", "Spado"))

write.table(CHAOS, "Data/BaristaSeq/IntergrationRe/Metric/CHAOS.csv", col.names = T,row.names = F,sep = ",",quote = F)
write.table(PAS, "Data/BaristaSeq/IntergrationRe/Metric/PAS.csv", col.names = T,row.names = F,sep = ",",quote = F)


###3.clustering metric plot----
CHAOS[which(CHAOS$Model == "GraphSTwithPASTE"), ]$Model <- "GraphST-PASTE"
PAS[which(PAS$Model == "GraphSTwithPASTE"), ]$Model <- "GraphST-PASTE"
model_metric[which(model_metric$model == "GraphSTwithPASTE"), ]$model <- "GraphST-PASTE"

color <- c("Banksy" = "#509C3D", "CellCharter" = "#C33931", "CN" = "#58BACC", "GraphST" = "#3A73AE", "GraphST-PASTE" = "#EE8435", "MENDER" = "#8D68B8", "NicheCompass" ="#BABD45", "Spado" = "#D57BBE" )

CHAOS_long <- melt(CHAOS, id.vars = "Model", variable.name = "slices", value.name = "value")
PAS_long <- melt(PAS, id.vars = "Model", variable.name = "slices", value.name = "value")

p1 <- ggplot(model_metric, aes(x = model, y = model_ari, col = model)) +
  geom_segment(aes(x = model, xend = model, y = 0, yend = model_ari), 
               arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.2, color = "#74A8D2") +
  labs(y = "ARI", x  = NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3))+
  NoLegend()

p2 <- ggplot(model_metric, aes(x = model, y = model_nmi, col = model)) +
  geom_segment(aes(x = model, xend = model, y = 0, yend = model_nmi), 
               arrow = arrow(length = unit(2, "mm"), type = "closed"), size = 0.2, color = "#C74B3E") +
  labs(y = "NMI", x  = NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3))+
  NoLegend()


p3 <- ggplot(CHAOS_long, aes(x = Model, y = value, fill = Model)) +
  geom_boxplot(lwd = 0.3, alpha = 1,outlier.size = 0.3, outlier.shape = 16) + 
  geom_jitter(aes(color = Model), position = position_jitterdodge(jitter.width = 0.15), 
              size = 0.2) +
  scale_color_manual(values = color)+
  scale_fill_manual(values = color) +
  labs(x = NULL, y = "CHAOS") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3))+
  NoLegend()


p4 <- ggplot(PAS_long, aes(x = Model, y = value, fill = Model)) +
  geom_boxplot(lwd = 0.3, alpha = 1,outlier.size = 0.3, outlier.shape = 16) + 
  geom_jitter(aes(color = Model), position = position_jitterdodge(jitter.width = 0.15), 
              size = 0.2) +
  scale_color_manual(values = color)+
  scale_fill_manual(values = color) +
  labs(x = NULL, y = "PAS") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.3))+
  NoLegend()

p_all <- wrap_plots(list(p1,p2,p3,p4), ncol = 4)
ggsave("Data/BaristaSeq/IntergrationRe/Metric/BaristaSeq_metric.pdf", plot = p_all, width = 22, height = 5, units = "cm")


###4. domain plotting---- 
slices <- unique(metadata$slices)
slices <- sort(slices)
model <- colnames(metadata)[5:ncol(metadata)]
domain <- unique(metadata$Ground_truth)
color <- c("VISp_I" = "#2F6DA1", "VISp_II/III" = "#EB49F7", "VISp_IV" = "#6CE3FB", "VISp_V" = "#EB594F", "VISp_VI" = "#FFFF54", "VISp_wm" ="#3B8749" )

plot_list <- list()
for(m in model){
  domain_data <- metadata[,match(c("X", "Y", "slices", m), colnames(metadata))]
  domain_plot_list <- list()
  for(s in slices){
    slices_data <- domain_data[domain_data$slices == s,]
    df <- cbind.data.frame(x_index = slices_data$X, y_index = slices_data$Y, predicted_domain = slices_data[, m])
    p <- ggplot(df, aes(x = x_index, y = y_index, col = predicted_domain )) +
      geom_point(size = 0.2) +
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
            #plot.title = element_blank(),
            plot.title = element_text(size = 24, hjust = 0.5)) +
      ggtitle(s) +
      coord_equal() 
    
    domain_plot_list[[paste(s,m,sep = "_")]] <- p
  }
  p_1 <- grid.arrange(grobs = domain_plot_list, ncol = 2)
  plot_list[[m]] <- p_1
}
p_all <- grid.arrange(grobs = plot_list, ncol = 9)
ggsave("Data/BaristaSeq/IntergrationRe/Metric/BaristaSeq_Result_domain.pdf", plot = p_all, width = 40, height = 6.75, units = "in")


###5. Umap plotting----
umapPlot <- function(originalEmb, metadata, col_slices, col_domain, isOriginal = FALSE){
  ###originalEmb 行是细胞 列是基因
  SeuratObj <- CreateSeuratObject(t(originalEmb),meta.data = metadata)
  SeuratObj <- NormalizeData(SeuratObj)
  all.genes <- rownames(SeuratObj)
  SeuratObj <- ScaleData(SeuratObj, features = all.genes)
  SeuratObj@assays[["RNA"]]@layers[["data"]] <- t(originalEmb)
  SeuratObj@assays[["RNA"]]@layers[["scale.data"]] <- t(originalEmb)
  SeuratObj <- RunPCA(SeuratObj, features = all.genes)
  if(dim(SeuratObj@reductions[["pca"]])[2] < 30){
    SeuratObj <- RunUMAP(SeuratObj, dims = 1:dim(SeuratObj@reductions[["pca"]])[2])
  }else{
    SeuratObj <- RunUMAP(SeuratObj, dims = 1:30)
  }
  
  
  SeuratObj$slices <- as.character(SeuratObj$slices)
  SeuratObj$original_domain <- as.character(SeuratObj$original_domain)
  
  if(isOriginal){
    p1 <- DimPlot(SeuratObj, reduction = "umap", group.by = "slices", pt.size = 0.01, cols = col_slices)+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, 'in'),
        legend.title = element_text(size = 8)
      )
    
    
    p2 <- DimPlot(SeuratObj, reduction = "umap", group.by = "original_domain", pt.size = 0.01, cols = col_domain)+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, 'in'),
        legend.title = element_text(size = 8)
      ) 
    
  }else{
    p1 <- DimPlot(SeuratObj, reduction = "umap", group.by = "slices", pt.size = 0.01, cols = col_slices)+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, 'in'),
        legend.title = element_text(size = 8)
      ) 
    
    p2 <- DimPlot(SeuratObj, reduction = "umap", group.by = "original_domain", pt.size = 0.01, cols = col_domain)+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, 'in'),
        legend.title = element_text(size = 8)
      ) 
    
  }
  
  return(list(SlicesUmap = p1, DomainUmap = p2))
}
###设置颜色----
col_slices <- c("1" = "#FFFF54", "2" = "#6CE3FB", "3" = "#EB49F7" )
col_domain <- c("VISp_I" = "#2F6DA1", "VISp_II/III" = "#EB49F7", "VISp_IV" = "#6CE3FB", "VISp_V" = "#EB594F", "VISp_VI" = "#FFFF54", "VISp_wm" ="#3B8749" )

###读取路径文件----
dataPath <- read.table("Data/BaristaSeq/IntergrationRe/Metric/domain_plot.txt")
dataPath$Model <- apply(dataPath, 1, function(x){
  unlist(strsplit(x, split = ":"))[1]
})
dataPath$Path <- apply(dataPath, 1, function(x){
  unlist(strsplit(x, split = ":"))[2]
})

dataPath <- dataPath[match(c("Banksy", "CellCharter", "GraphST", "GraphSTwithPASTE", "MENDER", "NicheCompass", "Spado") ,dataPath$Model),]


###original数据读取和绘图----
SlicesUmapList <- list()
DomainUmapList <- list()
originalData <- read_h5ad("Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad")
metadata <- originalData$obs

###model embedding数据读取和绘图----
for(i in 1:nrow(dataPath)){
  if(dataPath[i,]$Model == "GraphST" | dataPath[i,]$Model == "GraphSTwithPASTE"){
    preData <- read_h5ad(dataPath[i,]$Path) 
    preEmb <- preData$obsm$emb
  }else{
    preData <- read_h5ad(dataPath[i,]$Path)
    preEmb <- preData$obsm$X_embedding
  }
  
  prePlot <- umapPlot(preEmb, metadata, col_slices, col_domain, isOriginal = F)
  SlicesUmapList[[dataPath[i,]$Model]] <- prePlot$SlicesUmap
  DomainUmapList[[dataPath[i,]$Model]] <- prePlot$DomainUmap
}

SlicesUmap <- wrap_plots(SlicesUmapList, ncol = 7, heights = 1, widths = 1)
ggsave("Data/BaristaSeq/IntergrationRe/Metric/BaristaSeq_SlicesUmap.pdf", plot = SlicesUmap, width = 20, height = 2.5, units = "in")

DomainUmap <- wrap_plots(DomainUmapList, ncol = 7, heights = 1, widths = 1)
ggsave("Data/BaristaSeq/IntergrationRe/Metric/BaristaSeq_DomainUmap.pdf", plot = DomainUmap, width = 20, height = 2.5, units = "in")


