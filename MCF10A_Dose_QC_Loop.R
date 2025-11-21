# Improved and Annotated Code
library(Seurat)
library(tidyverse)
library(scCustomize)
library(HGNChelper)
library(scDblFinder)
set.seed(123)
getwd()

# Input and Output Paths
data.dir <- "~/EMT Model Cell lines/"
plot_dir <- "~/EMT Model Cell lines/Plots/"
save_dir <- "~/EMT Model Cell lines/"

gc()

# Read Data
dose.data <- Read10X("~/MCF10A_data_3")
seurat <- CreateSeuratObject(dose.data, project = "MCF10A_dose", min.features = 500)

# Annotation Data from Authors
annotation2 <- read.csv("~/KazuConcentration_GBC_annotation_Samples.csv", as.is = TRUE, check.names = FALSE)
names(annotation2) <- c("Cell", "GBC_annot", "Values")
seurat$Cell <- rownames(seurat@meta.data)

# Merge Metadata
merged <- left_join(seurat@meta.data, annotation2, by = "Cell")
rownames(merged) <- merged$Cell
stopifnot(all(rownames(seurat@meta.data) == rownames(merged)))  

seurat@meta.data <- merged


seurat$GBC_annot[is.na(seurat$GBC_annot)] <- "Unlabelled"

# Add Sample Metadata
seurat <- seurat %>%
  AddMetaData(object = ., metadata = case_when(
    .$GBC_annot == "GBC1" ~ "MCF10A_Control",
    .$GBC_annot == "GBC2" ~ "MCF10A_12.5pM",
    .$GBC_annot == "GBC3" ~ "MCF10A_25pM",
    .$GBC_annot == "GBC4" ~ "MCF10A_50pM",
    .$GBC_annot == "GBC8" ~ "MCF10A_100pM",
    .$GBC_annot == "GBC9" ~ "MCF10A_200pM",
    .$GBC_annot == "GBC12" ~ "MCF10A_400pM",
    .$GBC_annot == "GBC13" ~ "MCF10A_800pM",
    .$GBC_annot == "Unlabelled" ~ "MCF10A_Unlabelled",
    TRUE ~ as.character(.$GBC_annot)  # Default Case
  ), col.name = "Sample")

# Initialize QC Metrics and Filtered Objects
qc_metrics <- list()
filtered_objects <- new.env()
sample_names <- unique(seurat$Sample)

# Process Each Sample
for (sample in sample_names) {
  cat("Processing:", sample, "\n")
  
  # Subset Data
  cell_IDs <- seurat@meta.data %>% filter(Sample == sample) %>% pull(Cell)
  sample_data <- seurat[["RNA"]]$counts[, cell_IDs]
  sample_meta <- seurat@meta.data[cell_IDs, ]
  cells_before <- ncol(sample_data)
  
  # Filter Genes
  o <- order(rowSums(sample_data), decreasing = TRUE)
  sample_data <- sample_data[o, ]
  
  
  # Update Gene Symbols
  updated <- HGNChelper::checkGeneSymbols(rownames(sample_data), unmapped.as.na = FALSE, species = "human")
  updated$Suggested.Symbol <- make.unique(updated$Suggested.Symbol)
  rownames(sample_data) <- updated$Suggested.Symbol
  
  # Create Seurat Object
  sample_seurat <- CreateSeuratObject(sample_data, project = sample, min.features = 500, meta.data = sample_meta)
  sample_seurat <- scCustomize::Add_Cell_QC_Metrics(sample_seurat, add_top_pct = FALSE, species = "human")
  
  # Calculate QC Thresholds
  median_percent_MT <- median(sample_seurat$percent_mito)
  mad_percent_MT <- mad(sample_seurat$percent_mito)
  high_threshold_MT <- median_percent_MT + 3 * mad_percent_MT
  
  median_percent_Count <- median(sample_seurat$nCount_RNA)
  mad_percent_Count <- mad(sample_seurat$nCount_RNA)
  high_threshold_Count <- median_percent_Count + 3 * mad_percent_Count
  low_threshold_Count <- median_percent_Count - 3 * mad_percent_Count
  
  median_percent_Feature <- median(sample_seurat$nFeature_RNA)
  mad_percent_Feature <- mad(sample_seurat$nFeature_RNA)
  high_threshold_Feature <- median_percent_Feature + 3 * mad_percent_Feature
  low_threshold_Feature <- median_percent_Feature - 3 * mad_percent_Feature
  
  # Filter Cells
  sample_seurat <- subset(sample_seurat, subset = percent_mito < high_threshold_MT &
                                         nCount_RNA < high_threshold_Count &
                                         nFeature_RNA > low_threshold_Feature &
                                         nFeature_RNA < high_threshold_Feature)
  # Detect multiplets
  sce <- as.SingleCellExperiment(DietSeurat(sample_seurat))
  sce <- scDblFinder(sce)
  sample_seurat$Multiplet <- sce$scDblFinder.class
  cells_after <- ncol(sample_seurat)
  
  # Store QC Metrics
  qc_metrics[[sample]] <- data.frame(
    Sample = sample,
    high_threshold_MT = high_threshold_MT,
    high_threshold_Count = high_threshold_Count,
    low_threshold_Count = low_threshold_Count,
    high_threshold_Feature = high_threshold_Feature,
    low_threshold_Feature = low_threshold_Feature,
    Cells_Before_Filtering = cells_before,
    Cells_After_Filtering = cells_after
  )
  
  # Save Filtered Object
  assign(sample, sample_seurat, envir = filtered_objects)
}

# Combine and Save QC Metrics
qc_summary <- do.call(rbind, qc_metrics)
write.csv2(qc_summary, file = paste0("MCF10A_QC_Summary.csv"), row.names = FALSE)

# Save All Filtered Objects
save(list = ls(envir = filtered_objects), envir = filtered_objects, file = "Filtered_MCF10A_Dose.RData")
