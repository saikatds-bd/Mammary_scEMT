# Working remotely, so adding checkpoints
library(Seurat)
library(dplyr)
library(SCpubr)
library(scCustomize)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(plyr)
library(glmGamPoi)

setwd("/home/erikk/Desktop/Saikat_sc_data/Breast_Cancer_Patients/Tapsi_Data")


if (!file.exists("checkpoint_dietseurat.rds")) {
  seurat <- readRDS("/home/erikk/Desktop/Saikat_sc_data/Breast_Cancer_Patients/for_proteoglycan/A_single-cell_RNA_expression_atlas_of_normal_prene/Normal/Tapsi_Normal_Total.rds")
  gc()
  seurat <- DietSeurat(seurat, assays = "RNA")
  gc()
  seurat[["RNA"]] <- as(seurat[["RNA"]], "Assay5")
  gc()
  
  meta <- seurat@meta.data[, c("sample_id", "reported_diseases", "disease", "procedure_group")]
  unique_combos <- unique(meta) %>% filter(procedure_group %in% c("Reduction Mammoplasty"))
  seurat <- subset(seurat, subset = sample_id %in% unique_combos$sample_id)
  gc()
  
  seurat$Author_Annotation <- recode(as.character(seurat$broad_cell_type),
                                     "basal" = "Basal", "bcells" = "B-Cells", "fibroblasts" = "Fibroblast",
                                     "lumhr" = "LumHR", "lumsec" = "LumSec", "lymphatic" = "Lymphatic",
                                     "myeloid" = "Myeloid", "pericytes" = "Pericytes", "tcells" = "T-Cells",
                                     "vascular" = "Vascular")
  
  seurat$Final_Annotation <- ifelse(seurat$Author_Annotation %in% c("LumHR", "LumSec", "Basal"), "Epithelial", seurat$Author_Annotation)
  seurat$Patient_ID <- seurat$sample_id
  seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$Patient_ID)
  saveRDS(seurat, "checkpoint_dietseurat.rds")
} else {
  seurat <- readRDS("checkpoint_dietseurat.rds")
}

# Checkpoint: SCTransform
if (!file.exists("checkpoint_sctransform.rds")) {
  options(future.globals.maxSize = 128 * 1024^3)
  seurat <- SCTransform(seurat, verbose = TRUE, conserve.memory=TRUE)
  saveRDS(seurat, "checkpoint_sctransform.rds")
  
} else {
  seurat <- readRDS("checkpoint_sctransform.rds")
}

gc()

# Checkpoint: Dimensionality reduction + Harmony
if (!file.exists("checkpoint_harmony.rds")) {
  seurat <- RunPCA(seurat)
  gc()
  seurat <- RunUMAP(seurat, dims = 1:20, reduction = "pca", reduction.name = "PCA_UMAP")
  gc()
  seurat <- JoinLayers(seurat, assay = "RNA")
  gc()
  library(harmony)
  seurat <- RunHarmony(seurat, "Patient_ID")
  gc()
  seurat <- RunUMAP(seurat, dims = 1:20, reduction = "harmony", reduction.name = "Harmony_UMAP")
  gc()
  seurat <- PrepSCTFindMarkers(seurat, assay = "SCT", verbose = TRUE)
  save(seurat, file = "Processed_110pt.RData")
  saveRDS(seurat, "checkpoint_harmony.rds")
} else {
  seurat <- readRDS("checkpoint_harmony.rds")
}


#scATOMIC
unique_patients <- unique(seurat$Patient_ID)

for (patient_id in unique_patients) {
  # Subset the Seurat object for each unique Patient_ID
  subset_seurat <- subset(seurat, subset = Patient_ID == patient_id)
  
  # Set the default assay and convert to SingleCellExperiment
  DefaultAssay(subset_seurat) <- "RNA"

  gc()
  # Run scATOMIC for copy number prediction
  counts <- subset_seurat[["RNA"]]@counts
  cell_predictions <- run_scATOMIC(counts, mc.cores = 16, breast_mode = TRUE)
  gc()
  # Summarize the results
  results <- create_summary_matrix(
    prediction_list = cell_predictions, use_CNVs = FALSE, modify_results = TRUE, 
    mc.cores = 16, raw_counts = counts, min_prop = 0.5, breast_mode = TRUE
  )
  gc()
  # Generate the tree results
  tree_results_non_interactive <- scATOMICTree(
    predictions_list = cell_predictions, summary_matrix = results, 
    interactive_mode = FALSE, save_results = TRUE,
    project_name = paste0("Scatomic_", patient_id, "_Tumor")
  )
  gc()
  # Add the results as metadata
  subset_seurat <- AddMetaData(subset_seurat, results)
  subset_seurat$Cell_ID <- rownames(subset_seurat@meta.data)
  # Save subset metadata as a CSV file
  write.csv(subset_seurat@meta.data, paste0(patient_id, "_scatomic_metadata.csv"))
  rm(subset_seurat)
  gc()
  
}

seurat$Cell_ID <- rownames(seurat@meta.data)
metadata <- seurat@meta.data
files <- list.files(pattern = "scatomic_metadata\\.csv$", ignore.case = TRUE)


data_frames <- list()

# Loop over the files and read each one
for (file in files) {
  # Read the CSV file using read.csv2
  df <- read.csv(file, as.is = TRUE, check.names = FALSE, row.names = 1)
  # Append the dataframe to the list
  data_frames[[length(data_frames) + 1]] <- df
}

# Combine all dataframes into one using do.call with rbind
combined_df <- do.call(rbind, data_frames)
combined_df$Cell_ID <- rowanames(combined_df)

unique <- c("Cell_ID", setdiff(names(combined_df, names(metadata()))))
filtered_combined <- combined_df %>% dplyr::select(unique)

combined_metadata <- merge(metadata, filtered_combined, by="Cell_ID")
rownames(combined_metadata) <- combined_metadata$Cell_ID
combined_metadata <- combined_metadata[Cells(seurat),]
all(rownames(combined_metadata) == seurat$Cell_ID)

seurat@meta.data <- combined_metadata

#Import the BPE data
singler_ref <- celldex::BlueprintEncodeData()

#We will subset the ref data with major cell types
singler_ref.sub <- singler_ref[,singler_ref$label.main %in% c("Monocytes", "Macrophages", "DC",
                                                              "B-cells", "Fibroblasts", "CD8+ T-cells",
                                                              "CD4+ T-cells", "Epithelial cells", 
                                                              "Endothelial cells", "NK cells")]

singler_results_BPE_sub <- SingleR::SingleR(
  test = seurat[["SCT"]]$data,
  ref = singler_ref.sub,
  labels = singler_ref.sub@colData@listData$label.main
)
gc()

seurat@meta.data$Cell_Type_BPE <- singler_results_BPE_sub@listData$labels


seurat$Cell_Type_BPE[seurat$Cell_Type_BPE == "CD4+ T-cells" |
                       seurat$Cell_Type_BPE == "CD8+ T-cells"] <- "T-Cells"

seurat$Cell_Type_BPE <- recode(seurat$Cell_Type_BPE,
                               "Monocytes" = "Monocyte",
                               "T-cells" = "T-Cells",
                               "Fibroblasts" = "Fibroblast",
                               "Macrophages" = "Macrophage",
                               "Endothelial cells" = "Endothelial",
                               "B-cells" = "B-Cells",
                               "NK cells" = "NK",
                               "Epithelial cells" = "Epithelial",
                               "DC" = "DC"
)

rm(seurat.sce, singler_ref, singler_ref.sub, singler_results_BPE_sub)



bpal <- read.csv("Major_Cell_Type.csv", as.is = T, check.names = F, row.names = 1)
bpal <- bpal %>% filter(scATOMIC_pred %in% seurat$scATOMIC_pred)

seurat <- Add_Sample_Meta(seurat, meta_data = bpal, join_by_seurat = "scATOMIC_pred",
                          join_by_meta = "scATOMIC_pred")


library(tidyverse)
replace_patterns <- c("/", "Normal_Tissue", "Undefined", "Non-Stromal", "Non-Blood", "Blood", "Stromal")
seurat$General_Cell_Type <- seurat$Broad_Cell_Types


seurat$Cell_Type_SingleR <- seurat$Cell_Type_BPE

seurat$Cell_Type_BPE <- ifelse(seurat$Cell_Type_BPE %in% c("Macrophage", "Monocyte", "DC"),
                               "Myeloid", seurat$Cell_Type_BPE)

# Update Final_Annotation
seurat@meta.data <- seurat@meta.data %>%
  mutate(
    Final_Annotation = as.character(if_else(grepl(paste(replace_patterns, collapse = "|"), Annotation), Cell_Type_BPE, Annotation)),
    Final_Annotation = str_replace_all(Final_Annotation, c("cells" = "Cells", 
                                                           "Plasmablasts" = "Plasmablast"))
  )

seurat$Final_Annotation <- ifelse(seurat$Final_Annotation %in% c("Macrophage", "Monocyte", "DC"),
                                  "Myeloid", seurat$Final_Annotation)

qs_save(seurat, "110pt_Tapsi_Processed.qs2") # Same script used for the Pal dataset/normal samples


