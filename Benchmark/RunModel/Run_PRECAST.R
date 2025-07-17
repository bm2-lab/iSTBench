# Load necessary libraries
library(PRECAST)
library(Seurat)
library(reticulate)
library(SeuratDisk)
library(anndata)
library(SpatialExperiment)
library(pryr) 

# Main function to run the PRECAST pipeline
run_PRECAST <- function(input_file, output_file, sample, nclust, hvgs, runNormalization, tech, species, seed) {
  # input_file (abbreviation: i): Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.
  # 
  # output_file (abbreviation: o): Path to the directory where the output results will be saved. This is a required parameter.
  # 
  # sample (abbreviation: s): Name for the output result files. It is recommended to use the model name, for example, setting it as "PRECAST" will result in output files named "PRECAST.h5ad". This is a required parameter.
  # 
  # nclust (abbreviation: c): Number of identified domains. This is a required parameter.
  # 
  # hvgs (abbreviation: v): The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 2000.
  # 
  # runNormalization (abbreviation: n): A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.
  # 
  # tech (abbreviation: t): Sequencing technology that generated the data. Can be set to "ST", "Visium" or "Other_SRT". This is a required parameter.
  # 
  # species (abbreviation: p): The detected species type can be set to "Human", "Mouse" or "Unknown". This is an optional parameter, with the default set to "human".
  # 
  # seed (abbreviation: S): Random seed for model execution. This is an optional parameter, with the default set to 1234.
  
  # Step 1: Read and convert input .h5ad files into Seurat objects
  files <- list.files(input_file)
  seuList <- list()
  for (i in 1:length(files)) {
    file <- files[i]
    d_data <- read_h5ad(paste(input_file, file, sep = "/"))
    count <- t(d_data$X)
    metadata <- d_data$obs
    metadata$row <- d_data$obsm$spatial[,1]
    metadata$col <- d_data$obsm$spatial[,2]
    
    seu1 <- CreateSeuratObject(counts = count, meta.data = metadata, min.cells = 0, min.features = 0)
    colnames(seu1) <- seu1$barcode
    seuList[[i]] <- seu1
  }
  
  # Step 2: Merge multiple Seurat objects and join spatial layers
  data.seurat <- merge(seuList[[1]], y = seuList[-1], add.cell.ids = names(seuList))
  data.seurat <- JoinLayers(data.seurat)
  
  # Step 3: Create PRECAST object with or without normalization
  set.seed(seed)
  if(runNormalization){
    preobj <- CreatePRECASTObject(seuList = seuList, selectGenesMethod = "HVGs", gene.number = hvgs, numCores_sparkx = 10,
                                  premin.spots = 0, premin.features=0, postmin.spots=0, postmin.features=0)  
  }else{ 
    preobj <- CreatePRECASTObject(seuList = seuList, selectGenesMethod = "HVGs", gene.number = hvgs, numCores_sparkx = 10,
                                  premin.spots = 0, premin.features=0, postmin.spots=0, postmin.features=0)  
    
    for(i in 1:length(preobj@seulist)){
      preobj@seulist[[i]]@assays[["RNA"]]@layers[["data"]] <- preobj@seulist[[i]]@assays[["RNA"]]@layers[["counts"]]
    }
  }
  
  # Step 4: Record memory and time before integration
  start_time <- Sys.time()
  start_mem <- pryr::mem_used()
  
  # Step 5: Run PRECAST method
  PRECASTObj <- AddAdjList(preobj, platform = tech)
  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = TRUE, coreNum = 10, verbose = TRUE)
  PRECASTObj <- PRECAST(PRECASTObj, K = nclust)
  
  # Step 6: Model selection and integration
  PRECASTObj <- SelectModel(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species = species)
  
  end_time <- Sys.time()
  end_mem <- pryr::mem_used()
  mem_used <- end_mem - start_mem
  time_taken <- end_time - start_time
  
  # Save memory and time usage
  log_df <- data.frame(
    Memory_Usage_MiB = as.numeric(mem_used) / (1024 ^ 2),  # 转换为MB
    Execution_Time_s = as.numeric(time_taken, units = "secs")
  )
  cat("Time taken for PRECAST: ", time_taken, "\n")
  cat("Memory used for PRECAST: ", mem_used, "\n")
  
  # Step 8: Prepare Seurat object for saving
  interCell <- intersect(colnames(seuInt), colnames(data.seurat))
  data.seurat <- data.seurat[,interCell]
  
  embedding <- as.matrix(seuInt@reductions[["PRECAST"]]@cell.embeddings)
  rownames(embedding) <- colnames(data.seurat)
  colnames(embedding) <- paste0("embedding_", 1:ncol(embedding))
  
  data.seurat[["RNA3"]] <- as(object = data.seurat[["RNA"]], Class = "Assay")
  DefaultAssay(data.seurat) <- "RNA3"
  data.seurat[["RNA"]] <- NULL
  data.seurat <- RenameAssays(object = data.seurat, RNA3 = 'RNA')
  
  data.seurat$predicted_domain = as.character(seuInt$cluster)
  embedding <- CreateDimReducObject(embeddings = embedding, key = "embedding_", assay = "RNA")
  data.seurat[["embedding"]] <- embedding
  
  # Step 9: Save Seurat object and performance log
  h5Seurat_file_1 <- paste(output_file, paste(sample, "h5Seurat", sep = "."), sep = "/")
  SaveH5Seurat(data.seurat, filename = h5Seurat_file_1, overwrite = T)
  Convert(h5Seurat_file_1, dest = "h5ad", overwrite = T)
  
  write.table(log_df, paste(output_file, paste(sample, "_stats.csv", sep = ""), sep = "/"), col.names = T, row.names = F,sep = ",",quote = F)
}



# ---------------------- #
# Command-line interface
# ---------------------- #

# Define command-line arguments
library(getopt)
opt = matrix(c("input_file", "i", 1, "character", "Input file path",
               "output_file", "o", 1, "character", "Output file path",
               "sample", "s", 1, "character", "Model name",
               "nclust", "c", 1, "numeric", "Cluster number",
               "hvgs", "h", 2, "numeric", "Number of HVGs (default: 2000)",
               "runNormalization", "n", 2, "logical", "Run data normalization (default: TRUE)",
               "tech", "t", 1, "character", "Visuim, ST and Other_SRT",
               "species", "p", 1, "character", "Human, Mouse and Unknown",
               "seed", "S", 2, "numeric", "Seed (default: 1234)"),
             byrow = TRUE, ncol = 5)

# Parse arguments
args = getopt(opt)

if ( is.null(args$hvgs) ) {
  args$hvgs = 2000
}

if ( is.null(args$runNormalization) ) {
  args$runNormalization = TRUE
}

if (is.null(args$seed)) {
  args$seed = 1234
}

input_file <- args$input_file
output_file <- args$output_file
sample <- args$sample
nclust <- as.numeric(args$nclust)
hvgs <- args$hvgs
runNormalization <- args$runNormalization
tech <- args$tech
species <- args$species
seed <- args$seed

# Run the pipeline
run_PRECAST(input_file, output_file, sample, nclust, hvgs, runNormalization, tech, species, seed)
