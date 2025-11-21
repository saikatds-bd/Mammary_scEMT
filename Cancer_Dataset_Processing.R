
library(Seurat)
library(dplyr)
library(SCpubr)
library(scDECAF)
library(scCustomize)
library(decoupleR)
library(SingleR)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(scATOMIC)
library(plyr)
library(scDblFinder)
library(SingleCellExperiment)
# Load the Seurat object
load("~/Bhupinder_2021.processed.seurat.RData")
bpal <- DietSeurat(processed.seurat, assays = "RNA")
rm(processed.seurat)
bpal <- UpdateSeuratObject(bpal)
table(bpal$patient)
bpal$Patient_ID <- bpal$patient
bpal$patient <- NULL
bpal$Source_and_Patient <- paste0("BPal_", bpal$Patient_ID)


load("~/Liu_2022.processed.seurat.RData")
liu <- DietSeurat(processed.seurat, assays = "RNA")
rm(processed.seurat)
liu <- UpdateSeuratObject(liu)
table(liu$patient)
liu$Patient_ID <- liu$patient
liu$patient <- NULL
liu$Source_and_Patient <- paste0("Liu_", liu$Patient_ID)

load("~/Wu_2021.processed.seurat.RData")
wu <- DietSeurat(processed.seurat, assays = "RNA")
rm(processed.seurat)
gc()
wu <- UpdateSeuratObject(wu)
table(wu$patient)
wu$Patient_ID <- wu$patient
wu$patient <- NULL
wu$Source_and_Patient <- paste0("Wu_", wu$Patient_ID)

load("~/Xu_2021.processed.seurat.RData")
xu <- DietSeurat(processed.seurat, assays = "RNA")
rm(processed.seurat)
gc()

xu <- UpdateSeuratObject(xu)
table(xu$patient)
xu$Patient_ID <- xu$patient
xu$patient <- NULL
xu$Source_and_Patient <- paste0("Xu_", xu$Patient_ID)

seurat_objects <- list(bpal, liu, wu, xu)

common_metadata_columns <- Reduce(intersect, lapply(seurat_objects, function(x) colnames(x@meta.data)))

common_genes <- Reduce(intersect, lapply(seurat_objects, rownames))

seurat_objects <- lapply(seurat_objects, function(x) {
  x <- x[common_genes, ]
  return(x)
})
rm(bpal, liu, wu, xu, common_genes, common_metadata_columns)
gc()


seurat <- Merge_Seurat_List(seurat_objects,
                            project = "EMT_Signature")
rm(seurat_objects)
gc()

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


singler_ref <- celldex::BlueprintEncodeData()


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

qs_save(seurat, "55pt_Pang_Processed.qs2") 







































unique_patients <- unique(seurat$Patient_ID)

for (patient_id in unique_patients) {
  
  subset_seurat <- subset(seurat, subset = Patient_ID == patient_id)
  
  
  DefaultAssay(subset_seurat) <- "RNA"
  sce <- as.SingleCellExperiment(subset_seurat)
  
  
  sce <- scDblFinder(sce)
  subset_seurat$Multiplet <- sce$scDblFinder.class
  rm(sce)
  gc()
  
  # Run scATOMIC
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


#From authors (receptor based subtypes)
metadata <- seurat@meta.data %>% select(Patient_ID, orig.ident)
subtype <- distinct(metadata, Patient_ID, .keep_all = T)
write.csv2(subtype, "All_Subtype.csv", row.names = F)



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

qs_save(seurat, "55pt_Pang_Processed.qs2")


