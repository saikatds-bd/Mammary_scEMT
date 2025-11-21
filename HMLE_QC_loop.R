library(Seurat)
library(tidyverse)
library(scCustomize)
library(scDblFinder)
library(HGNChelper)

# Input files for HMLE_Epi and HMLE_Mes
input_files <- list(
  HMLE_Epi = "~/filtered_feature_bc_matrix",
  HMLE_Mes = "~/filtered_feature_bc_matrix"
)

# Output directories
data.dir <- "~/EMT Model Cell lines/"
plot_dir <- "~/Plots/"
save_dir <- "~/EMT Model Cell lines/"

# Initialize an empty list to store QC metrics
qc_metrics <- list()

# Initialize an empty environment to store filtered Seurat objects
filtered_objects <- new.env()

# Loop through the input files
for (sample_name in names(input_files)) {
  cat("Processing:", sample_name, "\n")
  
  # Read the data
  data <- Read10X(input_files[[sample_name]])
  
  # Filter genes expressed in at least 1% of cells
  o <- order(rowSums(data), decreasing = TRUE)
  data <- data[o, ]

  # Update gene symbols
  updated <- HGNChelper::checkGeneSymbols(rownames(data), unmapped.as.na = FALSE, species = "human")
  updated$Suggested.Symbol <- make.unique(updated$Suggested.Symbol)
  rownames(data) <- updated$Suggested.Symbol
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(data, project = sample_name, min.features = 500)
  seurat_obj <- scCustomize::Add_Cell_QC_Metrics(seurat_obj, add_top_pct = FALSE, species = "human")
  gc()
  # Calculate QC thresholds using 3 MAD
  median_percent_MT <- median(seurat_obj$percent_mito)
  mad_percent_MT <- mad(seurat_obj$percent_mito)
  high_threshold_MT <- median_percent_MT + 3 * mad_percent_MT
  
  median_percent_Count <- median(seurat_obj$nCount_RNA)
  mad_percent_Count <- mad(seurat_obj$nCount_RNA)
  high_threshold_Count <- median_percent_Count + 3 * mad_percent_Count
  low_threshold_Count <- median_percent_Count - 3 * mad_percent_Count
  
  median_percent_Feature <- median(seurat_obj$nFeature_RNA)
  mad_percent_Feature <- mad(seurat_obj$nFeature_RNA)
  high_threshold_Feature <- median_percent_Feature + 3 * mad_percent_Feature
  low_threshold_Feature <- median_percent_Feature - 3 * mad_percent_Feature
  gc()
  # Generate and save QC plots
  QC_plot_1 <- QC_Plot_UMIvsGene(seurat_object = seurat_obj, 
                                 meta_gradient_name = "percent_mito", 
                                 low_cutoff_gene = low_threshold_Feature,
                                 high_cutoff_UMI = high_threshold_Count,
                                 meta_gradient_low_cutoff = high_threshold_MT) + 
    ggtitle(sample_name) + 
    labs(caption = "Used 3 MAD for thresholds")
  
  pdf(file = paste0(plot_dir, sample_name, "_QC.pdf"), height = 5, width = 11)
  print(QC_plot_1)
  dev.off()
  
  # Subset the data based on QC thresholds
  seurat_obj <- subset(seurat_obj, subset = percent_mito < high_threshold_MT & 
                         nCount_RNA < high_threshold_Count & 
                         nFeature_RNA > low_threshold_Feature & 
                         nFeature_RNA < high_threshold_Feature)
  gc()
  # Detect multiplets
  sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
  sce <- scDblFinder(sce)
  seurat_obj$Multiplet <- sce$scDblFinder.class
  gc()
  # Add metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = ifelse(grepl("Epi", sample_name), "Epithelial", "Mesenchymal"), col.name = "Cell_State")
  seurat_obj <- AddMetaData(seurat_obj, metadata = sample_name, col.name = "Sample")
  seurat_obj <- AddMetaData(seurat_obj, metadata = "HMLE", col.name = "Cell_line")
  gc()
  # Store QC metrics
  qc_metrics[[sample_name]] <- data.frame(
    Sample = sample_name,
    high_threshold_MT = high_threshold_MT,
    high_threshold_Count = high_threshold_Count,
    low_threshold_Count = low_threshold_Count,
    high_threshold_Feature = high_threshold_Feature,
    low_threshold_Feature = low_threshold_Feature,
    Cells_Before_Filtering = ncol(data),
    Cells_After_Filtering = ncol(seurat_obj)
  )
  gc()
  # Store the filtered object in the environment
  assign(sample_name, seurat_obj, envir = filtered_objects)
  
  # Clean up
  rm(sce, seurat_obj, data)
  gc()
}

# Combine QC metrics into a single dataframe and save
qc_summary <- do.call(rbind, qc_metrics)
write.csv(qc_summary, file = paste0(save_dir, "HMLE_QC_Summary.csv"), row.names = FALSE)

# Save all filtered Seurat objects in a combined RData file
save(list = ls(envir = filtered_objects), envir = filtered_objects, file = paste0(save_dir, "Filtered_HMLE.RData"))


