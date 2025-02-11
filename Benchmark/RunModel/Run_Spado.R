library(SpaDo)
library(Seurat)
library(magrittr)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(tidyr)
library(fastcluster)
library(pryr)  
library(SeuratData)


run_spatial_analysis <- function(input_file, output_file, sample, hvgs, n_pcs , nclust, runNormalization ) {
  
  data <- read_h5ad(input_file)
  exp <- data$X
  metadata <- data$obs
  rownames(exp) <- metadata$barcode
  data <- CreateSeuratObject(counts = t(exp), min.cells = 0, min.features = 0, meta.data = metadata)
  data$slices <- as.character(data$slices)
  
  if (runNormalization) {
    data <- FindVariableFeatures(data, nfeatures = hvgs)
    data <- NormalizeData(data)
    data <- ScaleData(data)
    data <- RunPCA(data, npcs = n_pcs)
    data <- FindNeighbors(data, dims = 1:n_pcs)
    data <- FindClusters(data, resolution = 2, random.seed = 666)
  } else {
    data <- FindVariableFeatures(data, nfeatures = ncol(exp))
    data <- NormalizeData(data)
    data <- ScaleData(data)
    data@assays[["RNA"]]@layers[["scale.data"]] <- t(exp)
    data <- RunPCA(data, npcs = n_pcs)
    data <- FindNeighbors(data, dims = 1:n_pcs)
    data <- FindClusters(data, resolution = 2, random.seed = 666)
  }
  
  coor_list <- list()
  cellAn_list <- list()
  slices <- unique(data$slices)
  
  for (s in slices) {
    slices_data_seurat <- subset(data, slices == s)
    coor_data <- cbind(slices_data_seurat$X, slices_data_seurat$Y)
    rownames(coor_data) <- slices_data_seurat$barcode
    colnames(coor_data) <- c("X", "Y")
    
    if ("celltype" %in% colnames(slices_data_seurat@meta.data)) {
      celltype <- as.character(slices_data_seurat$celltype)
      names(celltype) <- slices_data_seurat$barcode
    } else {
      celltype <- as.character(slices_data_seurat$seurat_clusters)
      names(celltype) <- slices_data_seurat$barcode
    }
    
    coor_list[[s]] <- as.data.frame(coor_data)
    cellAn_list[[s]] <- celltype
  }
  
  start_time <- Sys.time()
  start_mem <- pryr::mem_used()
  
  cell_type_distribution_multiple <- SpatialCellTypeDistribution_multiple(
    sample_information_coordinate_list = coor_list,
    sequence_resolution = 'single_cell',
    sample_information_cellType_list = cellAn_list
  )
  
  cell_type_distribution = cell_type_distribution_multiple$cell_type_distribution_combine
  DD <- DistributionDistance(
    cell_type_distribution,
    distance = "JSD", no_cores = 10
  )
  
  domain_hclust <- DomainHclust(distribution_distance = DD, autoselection = FALSE, domain_num = nclust)
  
  end_time <- Sys.time()
  end_mem <- pryr::mem_used()
  
  time_taken <- end_time - start_time
  mem_used <- end_mem - start_mem
  
  re <- domain_hclust[["hclust_result_df"]]
  rownames(re) <- sapply(rownames(re), function(x) {
    unlist(strsplit(x, split = "_"))[2]
  })
  
  index <- match(as.character(data$barcode), rownames(re))
  re <- re[index,]
  predicted_domain <- re[, ncol(re)]
  
  data_seurat <- data
  data_seurat$predicted_domain <- predicted_domain
  data_seurat$original_domain <- as.character(data_seurat$original_domain)
  data_seurat$slices <- as.character(data_seurat$slices)
  
  data_seurat[["RNA3"]] <- as(object = data_seurat[["RNA"]], Class = "Assay")
  DefaultAssay(data_seurat) <- "RNA3"
  data_seurat[["RNA"]] <- NULL
  data_seurat <- RenameAssays(object = data_seurat, RNA3 = 'RNA')
  
  embedding_exp <- cell_type_distribution
  rowname_embedding_exp<- sapply(rownames(embedding_exp), function(x){
    unlist(strsplit(x, split = "_"))[2]
  })
  rownames(embedding_exp) <- rowname_embedding_exp
  colnames(embedding_exp) <- paste0("embedding_", 1:ncol(embedding_exp))
  embedding <- CreateDimReducObject(embeddings = embedding_exp, key = "embedding_", assay = "RNA")
  data_seurat[["embedding"]] <- embedding
  
  h5Seurat_file_1 <- paste(output_file, paste(sample, "h5Seurat", sep = "."), sep = "/")
  SaveH5Seurat(data_seurat, filename = h5Seurat_file_1, overwrite = TRUE)
  Convert(h5Seurat_file_1, dest = "h5ad", overwrite = TRUE)
  
  cat("Time taken: ", time_taken, "\n")
  cat("Memory used: ", mem_used, "\n")
  
  log_df <- data.frame(
    Memory_Usage_MiB = as.numeric(mem_used) / (1024 ^ 2), 
    Execution_Time_s = as.numeric(time_taken, units = "secs")
  )
  
  write.table(log_df, paste(output_file, paste(sample, "_stats.csv", sep = ""), sep = "/"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
}


library(getopt)
spec <- matrix(
  c(
    'input_file',     'i', 1, 'character',  'Path to input file',
    'output_file',    'o', 1, 'character',  'Path to output directory',
    'sample',         's', 1, 'character',  'Sample name',
    'hvgs',           'v', 2, 'integer',    'Number of highly variable genes (default: 2000)',
    'n_pcs',          'p', 2, 'integer',    'Number of principal components (default: 50)',
    'nclust',         'c', 1, 'integer',    'Number of clusters',
    'runNormalization', 'n', 2, 'logical',  'Flag to run normalization (default: TRUE)'
  ),
  byrow = TRUE,
  ncol = 5
)

opt <- getopt(spec)


hvgs <- ifelse(is.null(opt$hvgs), 2000, opt$hvgs)
n_pcs <- ifelse(is.null(opt$n_pcs), 50, opt$n_pcs)
runNormalization <- ifelse(is.null(opt$runNormalization), TRUE, opt$runNormalization)

run_spatial_analysis(opt$input_file, opt$output_file, opt$sample, hvgs, n_pcs, opt$nclust, runNormalization)
