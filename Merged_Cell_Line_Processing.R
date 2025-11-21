library(Seurat)
library(tidyverse)
library(limma)
library(clustree)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(scCustomize)
library(SCpubr)
library(patchwork)
library(scCustomize)
library(reticulate)
library(harmony)
library(msigdbr)

load("~/Cell_lines_EMT_Project_Colors_New.RData")

load("~/Filtered_HMLE.RData")
load("~/Filtered_D492.RData")

load("~/Filtered_MCF10A_Dose.RData")

set.seed(123)

data.dir <- paste0(getwd())
plot_dir <- paste0(data.dir, "/Plots/")
save_dir <- data.dir

gc()

seurat_objects <- list(D492_Epi, D492_Mes,  HMLE_Epi, HMLE_Mes, 
                       MCF10A_100pM, MCF10A_12.5pM, MCF10A_200pM, MCF10A_25pM, 
                       MCF10A_400pM, MCF10A_50pM, MCF10A_800pM, MCF10A_Control)


gene_lists <- lapply(seurat_objects, function(obj) {
  rownames(obj[["RNA"]]$counts)
})

common_genes <- Reduce(intersect, gene_lists)
seurat_objects <- lapply(seurat_objects, function(obj) {
  subset(obj, features = common_genes)
})

names(seurat_objects) <- c("D492_Epi", "D492_Mes",  "HMLE_Epi", "HMLE_Mes", 
                           "MCF10A_100pM", "MCF10A_12.5pM", "MCF10A_200pM", "MCF10A_25pM", 
                           "MCF10A_400pM", "MCF10A_50pM", "MCF10A_800pM", "MCF10A_Control")


seurat <- Merge_Seurat_List(seurat_objects, add.cell.ids = c("D492_Epi", "D492_Mes",  "HMLE_Epi", "HMLE_Mes", 
                                                             "MCF10A_100pM", "MCF10A_12.5pM", "MCF10A_200pM", "MCF10A_25pM", 
                                                             "MCF10A_400pM", "MCF10A_50pM", "MCF10A_800pM", "MCF10A_Control"),
                            project = "EMT_Models")

gc()

gc()
rm(D492_Epi, D492_Mes, HMLE_Epi, HMLE_Mes, 
   MCF10A_100pM, MCF10A_12.5pM, MCF10A_200pM, MCF10A_25pM, 
   MCF10A_400pM, MCF10A_50pM, MCF10A_800pM, MCF10A_Control, seurat_objects)

seurat <- JoinLayers(seurat, assay = "RNA")



p1 <- QC_Plots_Genes(seurat_object = seurat, group.by = "Sample")
p1
p2 <- QC_Plots_UMIs(seurat_object = seurat,group.by = "Sample")
p2
p3 <- QC_Plots_Mito(seurat_object = seurat, group.by = "Sample")
p3

pdf(paste0(plot_dir, "COmbined_Merged_QC_Plot.pdf"), height = 5, width = 11)
print(p1)
print(p2)
print(p3)
dev.off()

DefaultAssay(seurat) <- "RNA"

seurat$Cell_line <- ifelse(grepl("MCF10A", seurat$Sample), "MCF10A",
                           ifelse(grepl("HMLE", seurat$Sample), "HMLE", ifelse(grepl("D492", seurat$Sample), "D492", "Other")))

seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$Cell_line)
options(future.globals.maxSize = 16 * 1024^3)


seurat <- SCTransform(seurat, variable.features.n = 5000,vars.to.regress = c("percent_mito", "nCount_RNA",
                                                                             "nFeature_RNA",
                                                                             "S.Score", "G2M.Score"), verbose = TRUE)
gc()

sum(is.na(seurat$nCount_RNA))
sum(is.na(seurat$percent_mito))
summary(seurat$nCount_RNA)
summary(seurat$percent_mito)

seurat <- RunPCA(seurat)

seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "PCA_UMAP")

gc()

seurat <- IntegrateLayers(object = seurat, method = CCAIntegration, normalization.method = "SCT", verbose = F)
seurat <- RunUMAP(seurat, dims = 1:20, reduction = "integrated.dr")

seurat <- Joinlayers(seurat)
seurat <- RunHarmony(seurat, "Cell_line")
seurat <- RunUMAP(seurat, dims = 1:20, reduction = "harmony",reduction.name = "Harmony_UMAP")
p1 <- DimPlot_scCustom(seurat, group.by = "Sample", colors_use = EMT_colors, reduction = "PCA_UMAP", figure_plot = T)
p2 <- DimPlot_scCustom(seurat, group.by = "Cell_line", colors_use = EMT_colors, reduction = "umap", figure_plot = T)
p3 <- DimPlot_scCustom(seurat, group.by = "Sample", colors_use = EMT_colors, reduction = "umap", figure_plot = T)
p4 <- DimPlot_scCustom(seurat, group.by = "Sample", colors_use = EMT_colors, reduction = "Harmony_UMAP", figure_plot = T)
pdf(paste0(plot_dir, "Combined_Dimplots.pdf"), height = 4, width = 6)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

save(seurat,EMT_colors, file = "EMT_model_with_CCA.RData")