library(Seurat)
library(SCpubr)
library(scCustomize)
library(escape)
library(scRepertoire)
library(qs2)
library(tidyverse)

#Defining the gene list for mypepithelials. The Tapsi Kumar dataset had ensembl ids

myoepithelials <- list("Myoepithelial" = c("ENSG00000186847", "ENSG00000186081", "ENSG00000073282", "ENSG00000107796", 
                                           "ENSG00000149591", "ENSG00000101335", "ENSG00000091409", "ENSG00000132470", 
                                           "ENSG00000065618"))

all <- qs_read("110pt_Tapsi_Processed.qs2")
all_meta <- all@meta.data
epithelial <-subset(all, subset = Final_Annotation == "Epithelial")
rm(all)
gc()

epithelial <- droplevels(epithelial)
epi_metadata <- epi_metadata[Cells(epithelial),]
all(rownames(epi_metadata) == rownames(epithelial@meta.data))
epithelial@meta.data <- epi_metadata

DefaultAssay(epithelial) <- "RNA"
epithelial <- DietSeurat(epithelial, assays = "RNA", data = FALSE)

epithelial <- SCTransform(epithelial, verbose = TRUE, conserve.memory=T, variable.features.n = 5000)
gc()
epithelial <- RunPCA(epithelial)
gc()

epithelial <- RunUMAP(epithelial, dims = 1:20, reduction = "pca", reduction.name = "PCA_UMAP")
gc()

library(harmony)
epithelial <- RunHarmony(epithelial, "Patient_ID")
epithelial <- RunUMAP(epithelial, dims = 1:20, reduction = "harmony", reduction.name = "Harmony_UMAP")

epithelial <- FindNeighbors(epithelial, reduction = "harmony", dims = 1:20)
epithelial <- FindClusters(epithelial, resolution = c(0.2, 0.3), graph.name = "SCT_snn")

myo <- list("Myoepithelial"=EMT_Selected_ENSEMBL$Myoepithelial) #Containing the gene list
DefaultAssay(epithelial) <- "RNA"
epithelial <- NormalizeData(epithelial)
epithelial <- AddModuleScore(epithelial, features = myoepithelials, name = "_AMS")

epithelial <- Add_Sample_Meta(epithelial, meta_data = metadata_df_filtered, join_by_meta = "Cell_ID",
                              join_by_seurat = "Cell_ID")

epithelial$Myoepithelial_AMS <- epithelial$`_AMS1`

library(scCustomize)
p2 <- FeaturePlot_scCustom(epithelial, features = "Myoepithelial_AMS", figure_plot = T,  reduction = "Harmony_UMAP", na_cutoff = 0.4)
p2

# Based on the above plot
epithelial@meta.data <- epithelial@meta.data %>%
  mutate(Epi_Lineage = recode(
    SCT_snn_res.0.3,
    `1`  = "Myoepithelial", #
    `3`  = "Luminal Epithelial", 
    `2`  = "Luminal Epithelial", #
    `4`  = "Myoepithelial",#
    `5`  = "Luminal Epithelial", #
    
    `6`  = "Luminal Epithelial", #
    `7`  = "Luminal Epithelial", #
    `8`  = "Luminal Epithelial", #
    `9`  = "Luminal Epithelial", #
    `10` = "Myoepithelial",#
    `11` = "Luminal Epithelial", #
    `12` = "Myoepithelial",#
    `13` = "Luminal Epithelial", #
    `14` = "Luminal Epithelial", #
    `15` = "Luminal Epithelial",#
    `16` = "Luminal Epithelial" ,
    .default = "Unassigned"
  ))


epithelial$Refined_Broad <- ifelse(epithelial$Epi_Lineage %in% c("Luminal Progenitor", 
                                                                 "Mature Luminal", "Luminal (LP/ML Mixed)"),
                                   "Luminal Epithelial", "Myoepithelial")


p3 <- DimPlot_scCustom(epithelial, group.by = "Refined_Broad", figure_plot = T,  reduction = "Harmony_UMAP", colors_use = Cell_Type_Colors)
p3

epi_meta <- epithelial@meta.data %>% dplyr::select(Cell_ID, SCT_snn_res.0.3, 
                                       Myoepithelial_AMS, 
                                       Epi_Lineage, Refined_Broad)
merged_metadata <- left_join(all_metadata, epi_meta, by = "Cell_ID")

merged_metadata <- merged_metadata  %>%
  mutate(
    Refined_Final = case_when(
      Refined_Broad == "Luminal Epithelial" ~ "Luminal Epithelial",
      Refined_Broad  == "Myoepithelial" ~ "Myoepithelial",
      TRUE ~ Final_Annotation
    )
  )
rownames(merged_metadata) <- merged_metadata$Cell_ID

qs_save(merged_metadata, file = "20250923_Tapsi_all_metadata.qs2") # Same script was used for the Pal et al dataset
