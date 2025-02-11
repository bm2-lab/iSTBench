# Load necessary libraries
library(Banksy)
library(harmony)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Seurat)
library(SeuratDisk)
library(scran)
library(data.table)
library(scater)
library(cowplot)
library(ggplot2)
library(pryr)
library(getopt)

find_optimal_resolution <- function(data, target_clusters, tol = 1e-3, max_iter = 100) {
  lower_res <- 0.0
  upper_res <- 5.0
  optimal_res <- (lower_res + upper_res) / 2
  iter <- 0
  
  while (iter < max_iter) {
    iter <- iter + 1
    data <- clusterBanksy(data, dimred = "Harmony_BANKSY", resolution = optimal_res, seed = 1000)
    data_seurat = as.Seurat(data, data = NULL)
    i <- paste("clust_Harmony_BANKSY_k50_res", optimal_res, sep = "")
    num_clusters =length(unique(data_seurat@meta.data[, i]))
    
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

run_banksy <- function(input_file, output_file, sample, nclust, dataType, hvgs, runNormalization, k_geom, lambda, SEED) {
  # Load input data
  if (dataType == "RData") {
    load(input_file)
  } else if (dataType == "RDS") {
    data <- readRDS(input_file)
  } else {
    print("You should input data as RData or RDS format")
  }
  
  # Stagger spatial coordinates
  locs <- spatialCoords(data)
  locs <- cbind(locs, sample_id = factor(data$slices))
  locs_dt <- data.table(locs)
  colnames(locs_dt) <- c("sdimx", "sdimy", "group")
  locs_dt[, sdimx := sdimx - min(sdimx), by = group]
  global_max <- max(locs_dt$sdimx) * 1.5
  locs_dt[, sdimx := sdimx + group * global_max]
  locs <- as.matrix(locs_dt[, 1:2])
  rownames(locs) <- colnames(data)
  spatialCoords(data) <- locs
  
  if (runNormalization) {
    seu <- as.Seurat(data, data = NULL)
    seu <- FindVariableFeatures(seu, nfeatures = hvgs)
    
    # Normalize data
    scale_factor <- median(colSums(assay(data, "counts")))
    seu <- NormalizeData(seu, scale.factor = scale_factor, normalization.method = "RC")
    
    # Add data to SpatialExperiment and subset to HVGs
    aname <- "normcounts"
    assay(data, aname,withDimnames=FALSE) <- GetAssayData(seu)
    data <- data[VariableFeatures(seu), ]
  } else {
    seu <- as.Seurat(data, data = NULL)
    aname <- "normcounts"
    assay(data, aname, withDimnames = FALSE) <- GetAssayData(seu)
  }
  
  # The time and memory usage of computeBanksy were recorded
  start_time <- Sys.time()
  start_mem <- pryr::mem_used()
  
  # Running BANKSY
  compute_agf <- TRUE
  #k_geom <- 18
  data <- computeBanksy(data, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)
  
  #lambda <- 0.2
  npcs <- 20
  use_agf <- TRUE
  data <- runBanksyPCA(data, use_agf = use_agf, lambda = lambda, npcs = npcs, seed = SEED)
  
  
  end_time <- Sys.time()
  end_mem <- pryr::mem_used()
  
  time_taken <- end_time - start_time
  mem_used <- end_mem - start_mem
  
  log_df <- data.frame(
    Memory_Usage_MiB = as.numeric(mem_used) / (1024 ^ 2),  # 转换为MB
    Execution_Time_s = as.numeric(time_taken, units = "secs")
  )
  
  
  cat("Time taken for computeBanksy: ", time_taken, "\n")
  cat("Memory used for computeBanksy: ", mem_used, "\n")
  
  set.seed(SEED)
  harmony_embedding <- HarmonyMatrix(
    data_mat = reducedDim(data, paste("PCA_M1_lam", lambda, sep = "")),
    meta_data = colData(data),
    vars_use = c("slices"),
    do_pca = FALSE,
    max.iter.harmony = 20,
    verbose = FALSE
  )
  reducedDim(data, "Harmony_BANKSY") <- harmony_embedding
  
  data <- runBanksyUMAP(data, use_agf = TRUE, lambda = lambda, npcs = npcs)
  data <- runBanksyUMAP(data, dimred = "Harmony_BANKSY")
  
  # Cluster the Harmony corrected PCA embedding
  optimal_resolution = find_optimal_resolution(data,nclust)
  data <- clusterBanksy(data, dimred = "Harmony_BANKSY", resolution = optimal_resolution, seed = SEED)
  
  ###Add spatial information
  data@colData@listData[["X"]] <- data@int_colData@listData[["spatialCoords"]][,1]
  data@colData@listData[["Y"]] <- data@int_colData@listData[["spatialCoords"]][,2]
  
  # Transform to h5ad
  data_seurat = as.Seurat(data, data = NULL)
  i <- paste("clust_Harmony_BANKSY_k50_res", optimal_resolution, sep = "")
  data_seurat$predicted_domain = data_seurat@meta.data[, i]
  #seu$predicted_domain <- data_seurat$predicted_domain
  
  embedding_exp <- data@assays@data@listData[["H1"]]
  embedding_exp <- t(embedding_exp)
  colnames(embedding_exp) <- paste0("embedding_", 1:ncol(embedding_exp))
  embedding <- CreateDimReducObject(embeddings = embedding_exp, key = "embedding_", assay = "originalexp")
  data_seurat[["embedding"]] <- embedding
  
  data_seurat$original_domain <- as.character(data_seurat$original_domain)
  data_seurat$predicted_domain <- as.character(data_seurat$predicted_domain)
  data_seurat$slices <- as.character(data_seurat$slices)
  
  h5Seurat_file_1 <- paste(output_file, paste(sample, "h5Seurat", sep = "."), sep = "/")
  SaveH5Seurat(data_seurat, filename = h5Seurat_file_1)
  Convert(h5Seurat_file_1, dest = "h5ad")
  
  write.table(log_df, paste(output_file, paste(sample, "_stats.csv", sep = ""), sep = "/"), col.names = T, row.names = F,sep = ",",quote = F)
}



# Pass arguments from the command line
opt = matrix(c("input_file", "i", 1, "character", "Input file path",
               "output_file", "o", 1, "character", "Output file path",
               "sample", "s", 1, "character", "Sample name",
               "nclust", "c", 1, "numeric", "Cluster number",
               "dataType", "d", 2, "character", "the type of data storage(RData or RDS)",
               "hvgs", "h", 2, "numeric", "Number of HVGs (default: 2000)",
               "runNormalization", "n", 2, "logical", "Run data normalization (default: TRUE)",
               "k_geom", "k", 2, "numeric", "hiperparameter k",
               "lambda", "l", 2, "numeric", "hiperparameter lambda",
               "SEED", "S", 2, "numeric", "Random seed. default is 1000"),
             byrow = TRUE, ncol = 5)

args = getopt(opt)

if (is.null(args$dataType)) {
  args$dataType = "RData"
}

if (is.null(args$hvgs)) {
  args$hvgs = 2000
}

if (is.null(args$runNormalization)) {
  args$runNormalization = TRUE
}

if (is.null(args$k_geom)) {
  args$k_geom = c(15,30)
}

if (is.null(args$lambda)) {
  args$lambda = 0.8
}

if (is.null(args$SEED)) {
  args$SEED = 1000
}

input_file <- args$input_file
output_file <- args$output_file
sample <- args$sample
nclust <- as.numeric(args$nclust)
dataType <- args$dataType
hvgs <- args$hvgs
runNormalization <- args$runNormalization
k_geom <- args$k_geom
lambda <- args$lambda
SEED <- args$SEED

run_banksy(input_file, output_file, sample, nclust, dataType, hvgs, runNormalization, k_geom, lambda, SEED)
