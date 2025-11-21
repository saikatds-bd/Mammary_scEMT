library(Seurat)
library(SCpubr)
library(scCustomize)
library(qs2)
library(edgeR)
#Containing the Necessary Gene lists


emt_signature <- readxl::read_xlsx("Unfiltered_EMT_Signature.xlsx") #Gene list combining hallmark, sajib, and pan EMT
load("Gene_Lists_and_Colors.RData")

seurat <- qs_read("110pt_Tapsi_Processed.qs2")
metadata <- qs_read("20250923_Tapsi_all_metadata.qs2")

seurat <- droplevels(seurat)
DefaultAssay(seurat) <- "RNA"

seurat <- DietSeurat(seurat, assays = "RNA")


Final_anno_y <- Seurat2PB(seurat, sample = "Patient_ID", cluster = "Refined_Final")

get_object_name <- function(x) {
  base_name <- deparse(substitute(x))
  gsub("_y$", "", base_name)  # removing the trailing _y if present
}

process_pseudobulk <- function(pseudobulk_data) {
  data_name <- get_object_name(pseudobulk_data)
  
  pseudobulk_data$samples$Sample_Names <- rownames(pseudobulk_data$samples)
  
  
  keep.genes <- filterByExpr(pseudobulk_data, group = pseudobulk_data$samples$cluster)
  cat("Number of genes kept:", sum(keep.genes), "\n")
  
  
  pseudobulk_data <- pseudobulk_data[keep.genes, , keep.lib.sizes = FALSE]
  pseudobulk_data <- calcNormFactors(pseudobulk_data, method = "TMM")
  
  
  lcpm <- edgeR::cpm(pseudobulk_data, log = TRUE)
  corrected_lcpm <- as.data.frame(removeBatchEffect(lcpm, batch = pseudobulk_data$samples$sample))
  
  
  if (exists("emt_signature")) {
    emt_genes <- intersect(rownames(corrected_lcpm), emt_signature$ENSEMBL)
    corrected_lcpm <- corrected_lcpm[emt_genes, ]
  }
  
  
  group_means <- t(apply(corrected_lcpm, 1, function(x) tapply(x, pseudobulk_data$samples$cluster, mean, na.rm = TRUE)))
  group_means <- as.data.frame(group_means)
  group_means$Highest_Type <- colnames(group_means)[apply(group_means, 1, which.max)]
  names(group_means) <- paste0(data_name, "_", names(group_means))
  group_means$ensembl <- rownames(group_means)
  
  
  as.data.frame(group_means)
}
gc()

Final_anno_results <- process_pseudobulk(Final_anno_y)
names(Final_anno_results) <- gsub("pseudobulk_data", "Final_anno", names(Final_anno_results))

library(dplyr)
library(scCustomize)

final_emt <- merge(Final_anno_results, emt_signature, by.x = "ensembl", by.y = "ENSEMBL", all.x = TRUE)
final_emt$Gene <- final_emt$gene

Tapsi_Filtered_Epithelial <- final_emt %>%
  filter(Group == "Epithelial") %>%
  filter(Final_anno_Highest_Type %in% c("Luminal Epithelial")) %>%
  filter(`Final_anno_Epithelial` > 2) %>%
  pull(Gene)

#Common_Normal_Epithelial <- intersect(EMT_gene_list_Symbol[["Tapsi_Normal_Epithelial"]], Tapsi_Normal_Epithelial)
Tapsi_Extended_Mesenchymal <- final_emt %>%
  filter(Group == "Mesenchymal") %>%
  filter(Final_anno_Highest_Type == "Myoepithelial" | Final_anno_Highest_Type == "Fibroblast") %>%
  filter(Final_anno_Myoepithelial > 2 | Final_anno_Fibroblast > 2) %>%
  pull(Gene)


Tapsi_Filtered_Mesenchymal <- final_emt %>%
  filter(Author_Category == "Mesenchymal") %>%
  filter(Final_anno_Highest_Type == "Fibroblast") %>%
  filter(Final_anno_Fibroblast > 2) %>%
  pull(Gene)



Tapsi_genes_symbol <- list("Tapsi_Epithelial" = Tapsi_Filtered_Epithelial,
                           "Tapsi_Extended_Mesenchymal" = Tapsi_Extended_Mesenchymal,
                           "Tapsi_Mesenchymal" = Tapsi_Filtered_Myo_Fibro_Mesenchymal)


###BPal Data
library(Seurat)
library(SCpubr)
library(scCustomize)
library(qs2)
library(edgeR)
#Containing the Necessary Gene lists


emt_signature <- readxl::read_xlsx("Unfiltered_EMT_Signature.xlsx") #Gene list combining hallmark, sajib, and pan EMT
load("Gene_Lists_and_Colors.RData")

seurat <- qs_read("13pt_BPal_Processed.qs2")
metadata <- qs_read("BPal_all_metadata.qs2")

seurat <- droplevels(seurat)
DefaultAssay(seurat) <- "RNA"

seurat <- DietSeurat(seurat, assays = "RNA")


Final_anno_y <- Seurat2PB(seurat, sample = "Patient_ID", cluster = "Refined_Final")

get_object_name <- function(x) {
  base_name <- deparse(substitute(x))
  gsub("_y$", "", base_name)  # removing the trailing _y if present
}

process_pseudobulk <- function(pseudobulk_data) {
  data_name <- get_object_name(pseudobulk_data)
  
  pseudobulk_data$samples$Sample_Names <- rownames(pseudobulk_data$samples)
  
  
  keep.genes <- filterByExpr(pseudobulk_data, group = pseudobulk_data$samples$cluster)
  cat("Number of genes kept:", sum(keep.genes), "\n")
  
  
  pseudobulk_data <- pseudobulk_data[keep.genes, , keep.lib.sizes = FALSE]
  pseudobulk_data <- calcNormFactors(pseudobulk_data, method = "TMM")
  
  
  lcpm <- edgeR::cpm(pseudobulk_data, log = TRUE)
  corrected_lcpm <- as.data.frame(removeBatchEffect(lcpm, batch = pseudobulk_data$samples$sample))
  
  
  if (exists("emt_signature")) {
    emt_genes <- intersect(rownames(corrected_lcpm), emt_signature$ENSEMBL)
    corrected_lcpm <- corrected_lcpm[emt_genes, ]
  }
  
  
  group_means <- t(apply(corrected_lcpm, 1, function(x) tapply(x, pseudobulk_data$samples$cluster, mean, na.rm = TRUE)))
  group_means <- as.data.frame(group_means)
  group_means$Highest_Type <- colnames(group_means)[apply(group_means, 1, which.max)]
  names(group_means) <- paste0(data_name, "_", names(group_means))
  group_means$ensembl <- rownames(group_means)
  
  
  as.data.frame(group_means)
}
gc()

Final_anno_results <- process_pseudobulk(Final_anno_y)
names(Final_anno_results) <- gsub("pseudobulk_data", "Final_anno", names(Final_anno_results))

library(dplyr)
library(scCustomize)

final_emt <- merge(Final_anno_results, emt_signature, by.x = "ensembl", by.y = "ENSEMBL", all.x = TRUE)
final_emt$Gene <- final_emt$gene

BPal_Filtered_Epithelial <- final_emt %>%
  filter(Group == "Epithelial") %>%
  filter(Final_anno_Highest_Type %in% c("Luminal Epithelial")) %>%
  filter(`Final_anno_Epithelial` > 2) %>%
  pull(Gene)

#Common_Normal_Epithelial <- intersect(EMT_gene_list_Symbol[["BPal_Normal_Epithelial"]], BPal_Normal_Epithelial)
BPal_Extended_Mesenchymal <- final_emt %>%
  filter(Group == "Mesenchymal") %>%
  filter(Final_anno_Highest_Type == "Myoepithelial" | Final_anno_Highest_Type == "Fibroblast") %>%
  filter(Final_anno_Myoepithelial > 2 | Final_anno_Fibroblast > 2) %>%
  pull(Gene)

BPal_Filtered_Mesenchymal <- final_emt %>%
  filter(Author_Category == "Mesenchymal") %>%
  filter(Final_anno_Highest_Type == "Fibroblast") %>%
  filter(Final_anno_Fibroblast > 2) %>%
  pull(Gene)


BPal_genes_symbol <- list("BPal_Epithelial" = BPal_Filtered_Epithelial,
                           "BPal_Extended_Mesenchymal" = BPal_Extended_Mesenchymal,
                           "BPal_Mesenchymal" = BPal_Filtered_Myo_Fibro_Mesenchymal)

Final_EMT_List <- list("Epithelial" = intersect(BPal_genes_symbol$BPal_Epithelial,
                                                Tapsi_genes_symbol$Tapsi_Epithelial),
                       "Mesenchymal" = intersect(BPal_genes_symbol$BPal_Mesenchymal,
                                                 Tapsi_genes_symbol$Tapsi_Mesencymal))

save(Final_EMT_List, file = "Gene_Lists.RData")




