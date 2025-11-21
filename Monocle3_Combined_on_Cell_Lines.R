#Monocle3
library(monocle3)
library(msigdbr)
library(dplyr)
library(tibble)
library(fgsea)
library(Seurat)
library(presto)
library(cowplot)
library(scCustomize)
library(SeuratWrappers)
library(tradeSeq)
library(magrittr)
load("C:/Users/ssa214/UiT Office 365/O365-PhD_Saikat - General/Single Cell RNA-seq/EMT Project/Could_be_Final/20250205_EMT_model_with_CCA.RData")

data.dir <- getwd()
plot_dir <- paste0(data.dir, "/Plots/")
save_dir <- data.dir
library(SeuratWrappers)
seurat <- PrepSCTFindMarkers(seurat)

options(future.globals.maxSize = 16 * 1024^3) 

set.seed(123)
DefaultAssay(seurat) <- "SCT"
seurat <- FindNeighbors(seurat, reduction = "integrated.dr", dims = 1:20)
seurat <- FindClusters(seurat, resolution = 0.2)
seurat$SCT_snn_res.0.2 <- factor(as.numeric(seurat$SCT_snn_res.0.2))
seurat$seurat_clusters <- factor(as.numeric(seurat$seurat_clusters))
Idents(seurat) <- "SCT_snn_res.0.2"

library(presto)
de_genes <- RunPrestoAll(seurat)


m_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  mutate(gs_name = gsub("HALLMARK_", "", gs_name))
m_df$gs_name <- gsub("_", " ", m_df$gs_name)

fgsea_sets <- split(m_df$gene_symbol, m_df$gs_name)


fgsea_res_all <- data.frame(pathway=character(), pval=double(),
                            padj=double(), ES=double(),
                            NES=double(), nMoreExtreme=double(),
                            size=integer(), leadingEdge=list(),
                            time=character()
)

time_groups <- unique(de_genes$cluster)
for(time_point in time_groups){
  cluster1.genes <- de_genes %>%
    dplyr::filter(cluster == time_point) %>%
    mutate(rank_stats = sign(avg_log2FC) * -log10(p_val)) %>%
    filter(is.finite(rank_stats)) %>%
    arrange(desc(rank_stats)) %>%
    dplyr::select(gene, rank_stats) 
  ranks <- setNames(cluster1.genes$rank_stats, cluster1.genes$gene)
  
  fgseaRes <- fgseaMultilevel(pathways = fgsea_sets, stats = ranks, nproc = 1, eps = 0)
  
  fgseaResTidy <- as_tibble(fgseaRes) %>%
    arrange(desc(NES)) %>%
    mutate(time = time_point)
  
  fgsea_res_all <- bind_rows(fgsea_res_all, fgseaResTidy)
}

fgsea_res_all_plot <-  fgsea_res_all %>%
  filter(padj < 0.05) #%>%

library(ggpubr)
plot_title="Hallmark Pathway Enrichment"
x_lab_name='Gene expression clusters'
library(RColorBrewer)
obj_plot<-ggplot(fgsea_res_all_plot, aes(x=time, y=reorder(pathway,NES), fill= NES)) +
  geom_tile(aes(fill = NES),colour = "white")+

  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu")), name = "NES")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "right")+
  coord_fixed(ratio=0.4)+
  xlab("Gene Expression Clusters") + ylab("Pathways")
obj_plot

pdf(paste0(plot_dir, "Combined_FGSEA_Heatmap.pdf"), height = 6, width = 8)
print(obj_plot)
dev.off()

#########################Monocle3##################
DimPlot(seurat, label = T, group.by = "SCT_snn_res.0.2")

cds <- SeuratWrappers::as.cell_data_set(seurat,assay = "SCT")

set.seed(42)
cds@int_colData@listData$reducedDims$INTEGRATED.DR <- Embeddings(seurat, reduction = "integrated.dr")

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds_clusters <- clusters(cds)

list_cluster <- seurat@meta.data[[sprintf("seurat_clusters")]]

names(list_cluster) <- seurat@assays[["SCT"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


seurat$Monocle3_Clusters <- cds_clusters[Cells(seurat)]



plot_cells(cds, 
           color_cells_by = "cluster",
           group_label_size=4, graph_label_size=3, 
           label_groups_by_cluster = TRUE, label_cell_groups = T, label_roots = T,
           label_leaves = T, label_branch_points = T, label_principal_points = T)

cds <- learn_graph(cds, use_partition=TRUE, close_loop=FALSE)


root_plot <- plot_cells(cds, color_cells_by="cluster",
                group_label_size=4, graph_label_size=3.5,
                label_cell_groups=F, label_principal_points=TRUE,
                label_groups_by_cluster=FALSE) 
root_plot

# We find all the cells that are close to the starting point
cell_ids <- colnames(cds)[seurat$SCT_snn_res.0.2 ==  "3"]
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

##Visualization
library(viridis)
trajec_plot <- plot_cells(cds, color_cells_by = "cluster",label_groups_by_cluster = T,
                          label_roots = F,
                          label_leaves = F,
                          label_principal_points = F,
                          label_branch_points = F,
                          group_label_size = 5) 

pseudo_plot <- plot_cells(cds,
                          color_cells_by = "pseudotime",
                          group_cells_by = "cluster",
                          label_cell_groups = T,
                          label_groups_by_cluster=T,
                          label_leaves=FALSE,
                          label_branch_points=FALSE,
                          label_roots = FALSE,
                          trajectory_graph_color = "red",
                          group_label_size = 5)



cds$monocle3_pseudotime <- pseudotime(cds)

seurat <- AddMetaData(
  object = seurat,
  metadata = cds$monocle3_pseudotime,
  col.name = "Monocle3_Pseudotime"
)


pseudo_plot
trajec_plot

library(SCpubr)
library(ggpubr)
p1 <- do_BarPlot(seurat, split.by = "SCT_snn_res.0.2", group.by = "Sample", colors.use = EMT_colors, position = "fill")
p2 <- do_BoxPlot(seurat, feature = "Monocle3_Pseudotime", group.by = "SCT_snn_res.0.2", order = T, legend.position = "right") + theme_pubr()
p3 <- DimPlot_scCustom(seurat, group.by = "Sample", colors_use = EMT_colors, figure_plot = T, reduction = "umap")
p4 <- DimPlot_scCustom(seurat, group.by = "SCT_snn_res.0.2", colors_use = EMT_colors, figure_plot = T, reduction = "umap")

pdf(paste0(plot_dir,"EMT_CCA_Cell_lines_Dim_Plots.pdf"), height = 5, width = 7)
print(p2)
print(pseudo_plot)
print(trajec_plot)
dev.off()
gc()
save(seurat, cds, mst,cell_ids, closest_vertex, file = "Monocle3_CDS_EMT_Cell_lines_CCA.RData")
load(file = "Monocle3_CDS_EMT_Cell_lines_CCA.RData")


# Get the closest vertice for every cell
y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

root <- cds@principal_graph_aux$UMAP$root_pr_nodes

##################################Adding nodes back to seurat################
mono <- y_to_cells[,c(2,3)]
mono$Monocle3_Node <- paste0("Y_", mono$Y)

all(mono$cells==rownames(seurat@meta.data))

seurat@meta.data$Monocle3_Node <- mono$Monocle3_Node[match(rownames(seurat@meta.data), mono$cells)]

seurat_subset <- subset(seurat, subset = Monocle3_Node %in% c("Y_54", "Y_58", "Y_39", "Y_35", "Y_1"))

# Plot bar plot with only the selected Monocle3_Node groups
SCpubr::do_BarPlot(seurat_subset, 
                   group.by = "Sample", 
                   split.by = "Monocle3_Node", 
                   position = "fill")


SCpubr::do_BarPlot(seurat, group.by = "Sample", split.by = "Monocle3_Node", position = "fill", colors.use = cols)



seurat$Monocle3_Clusters <- cds_clusters[Cells(seurat)]
##################Continue with Tradeseq#######################


# Get the other endpoints
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

endpoints <- c("Y_54", "Y_58")
# For endpoint Y_14
cellWeights <- lapply(endpoints, function(endpoint) {
  # We find the path between the endpoint and the root
  path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path <- as.character(path)
  # We find the cells that map along that path
  df <- y_to_cells[y_to_cells$Y %in% path, ]
  df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()

rownames(cellWeights) <- colnames(cds)
non_zero_cells <- rowSums(cellWeights) != 0
cellWeights <- as.matrix(cellWeights[non_zero_cells, ]) #Taking the cells that are only in the Y_14 lineage

pseudo <- pseudotime(cds)
pseudo <- pseudo[rownames(cellWeights)]

pseudotime <- matrix(pseudo, ncol = ncol(cellWeights),
                     nrow = nrow(cellWeights), byrow = FALSE)

features <- as.data.frame(VariableFeatures(seurat))

features$rows <- rownames(features)


counts <- seurat[["SCT"]]$counts[,rownames(cellWeights)]


gc()

library(BiocParallel)
nCores <-16
param <- MulticoreParam(workers = nCores)
register(param)



source_batch <- seurat@meta.data[rownames(cellWeights),]
all(rownames(source_batch) == colnames(counts))
source_batch <- source_batch$Cell_line
batch_matrix <- model.matrix(~ source_batch)

set.seed(42)

all(colnames(counts) == rownames(cellWeights))

sce <- fitGAM(counts = counts, 
              pseudotime = pseudotime, 
              cellWeights = cellWeights, 
              nknots = 6, 
              BPPARAM = param,
              gene = features$`VariableFeatures(seurat)`, 
              U = batch_matrix,   
              verbose = TRUE)

save(sce, file="Fitgam_Y_54_Endpoint.RData")


startRes <- startVsEndTest(sce)

startRes$padj <- p.adjust(startRes$pvalue, method = "BH")
de <- startRes %>% filter(padj < 0.05)


save(sce, startRes, de, file="Fitgam_cell_line_monocle3.RData")


#Rerunning this to clean memory and start over again

load(file="Fitgam_2_Endpoints.RData")

load(file = "Filtered_Samples_DE_Results_treat.RData")

filtered_sig <- sig %>% filter(gene %in% rownames(de))

write.csv2(filtered_sig, "Saikat_EMT.csv")
