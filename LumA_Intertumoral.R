library(Seurat)
library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
load("Annotated_Pang.RData")

seurat <- subset(seurat, subset = Refined_Annotation == "Cancer" & Pam50_Subtype = "LumA")
# We will separate the patients into quartiles based on their epithelial score
add_quartiles_by_patient <- function(seurat,
                                     score_col = "Epithelial_AMS",
                                     group_col = "Patient_ID",
                                     new_col   = "EpiQuartile",
                                     labels    = c("Q1_lowest","Q2","Q3","Q4_highest"),
                                     min_cells = 20){
  
  meta <- seurat@meta.data
  
  stopifnot(score_col %in% names(meta), group_col %in% names(meta))
  
  
  meta <- meta %>%
    group_by(.data[[group_col]]) %>%
    group_modify(~{
      x <- .x[[score_col]]
      # if too few cells or all NA, return NA labels
      if (sum(!is.na(x)) < min_cells) {
        .x[[new_col]] <- factor(NA, levels = labels)
        return(.x)
      }
      rnk  <- rank(x, ties.method = "random", na.last = "keep")
      brks <- quantile(rnk, probs = c(0, 0.25, 0.50, 0.75, 1.0),
                       na.rm = TRUE, names = FALSE)
      qfac <- cut(rnk, breaks = brks, include.lowest = TRUE, labels = labels)
      .x[[new_col]] <- factor(qfac, levels = labels)
      .x
    }) %>%
    ungroup()
  
  
  seurat@meta.data[[new_col]] <- meta[[new_col]]
  seurat
}

seurat <- add_quartiles_by_patient(seurat)


patient_epi <- seurat@meta.data %>%
  group_by(Patient_ID) %>%
  summarize(median_epi = median(Common_Filtered_Epithelial_AMS, na.rm = TRUE))

patient_epi <- patient_epi %>%
  mutate(
    Quartile = cut(
      rank(median_epi, ties.method = "average"),
      breaks = quantile(rank(median_epi), probs = c(0, 0.25, 0.5, 0.75, 1)),
      include.lowest = TRUE,
      labels = c("Q1_lowest","Q2","Q3","Q4_highest")
    )
  )

writexl::write_xlsx(patient_epi, "Luminal A Patient_Quartiles.xlsx")
seurat$Patient_Quartile <- patient_epi$Quartile[match(seurat$Patient_ID, patient_epi$Patient_ID)]

Cell_Type_Colors["Q1_lowest"] <-  "#00468BFF"
Cell_Type_Colors["Q2"] <- "#DEB887"
Cell_Type_Colors["Q3"] <- "sienna1"
Cell_Type_Colors["Q4_highest"] <- "#E31A1C"

library(scCustomize)

p1 <- DimPlot_scCustom(seurat, group.by = "Patient_Quartile",
                       reduction = "Harmony_UMAP", colors_use = Cell_Type_Colors, figure_plot = T)
p1

p2 <- DimPlot_scCustom(seurat, group.by = "Patient_ID",
                       reduction = "Harmony_UMAP", figure_plot = T)

p3 <- DimPlot_scCustom(seurat, group.by = "SCT_snn_res.0.2",
                       reduction = "Harmony_UMAP", figure_plot = T, colors_use = Cell_Type_Colors, label = T)

p3

pdf("Intertumoral_UMAP.pdf", height = 5, width = 7)
print(p1)
print(p2)
print(p3)
dev.off()


Quartile_Colors <- c(
  "Q1_lowest"   = "#00468BFF",
  "Q2"          = "#DEB887",
  "Q3"          = "sienna1",
  "Q4_highest"  = "#E31A1C"
)

bar_cl <- do_BarPlot(seurat, group.by = "Patient_Quartile", split.by = "SCT_snn_res.0.2",
                     position = "fill",
                     colors.use = Quartile_Colors)

bar_pt <- do_BarPlot(seurat, group.by = "Patient_ID", split.by = "SCT_snn_res.0.2",
                     position = "fill")

pdf("Intertumoral_Composition.pdf", height = 5, width = 5)
print(bar_cl)
print(bar_pt)
dev.off()


###Pseudobulk for functional analysis

setwd("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Single Cell RNA-seq/EMT Project/20250814/Validation/Helpers")
source("MDS_plot.R")
source("Generate_barplot.R")
source("perform_enrichment_function.R")
source("Updated_Volcano_function.R")
source("perform_GSEA.R")
source("hallmark_gsea_plot_function.R")
source("EMT_Hallmark_Helper.R")

setwd("~/LumA Patients")
library(edgeR)

df_seurat <- subset(seurat, subset = Patient_Quartile %in% c("Q1_lowest", "Q4_highest"))
DefaultAssay(df_seurat) <- "RNA"
df_seurat <- DietSeurat(df_seurat, assays = "RNA")

gc()

y <- Seurat2PB(df_seurat, sample = "Patient_ID", cluster = "Patient_Quartile")
y$samples$Sample_Names <- rownames(y$samples)

keep.genes <- filterByExpr(y, group=y$samples$group)

y <- y[keep.genes, , keep=FALSE]
y <- calcNormFactors(y, method = "TMM")
lcpm <-edgeR::cpm(y, log = T)

gc()
cluster <- as.factor(y$samples$cluster)


library(ggpubr)
library(ggplot2)

lcpm <- edgeR::cpm(y)

p1 <- plotMDS_ggpubr(
  lcpm = lcpm, 
  sample_data = y$samples, 
  color_var = "cluster",
  label_var = "sample", 
  palette = "Set2",         
  title = "LumA MDS Plot",
  dims=c(1,2)
)
p1

pdf("MDS with all LumA.pdf", height = 7, width = 8)
print(p1)
dev.off()

design <- model.matrix(~ 0 + cluster, y$samples)

colnames(design) <- gsub("cluster", "", colnames(design))
design

vfit <- voomLmFit(y, design, plot=TRUE, sample.weights = T)
design

contr <- makeContrasts(Q1_lowest - Q4_highest, levels = colnames(coef(vfit)))
contr

cfit <- contrasts.fit(vfit, contrasts = contr)
efit <- eBayes(cfit, robust = T)

summary(decideTests(efit))

table_fit <- topTable(efit,coef="Q1_lowest - Q4_highest", n=Inf)

sig <- table_fit %>% 
  filter(adj.P.Val < 0.05) %>% 
  arrange(adj.P.Val)

names(table_fit) <- gsub("gene", "Symbol", names(table_fit))
writexl::write_xlsx(table_fit, "LumA_All_Patient_DE.xlsx")

plot_volcano(data=table_fit, plot_title = "Q1 vs Q4 - LumA Patients", 
             output_file = "Q1 vs Q4 - LumA Patients.svg", 
             height = 4, width = 5)




###GSEA to find difference between patients
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(dplyr)
library(stats)
library(readxl)
library(ggplot2)
library(viridis)

m_t2g_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)

go <- gsea_barplot_multi(
  df= table_fit,                        
  wrap_width = 30,
  use = "GO:BP",      
  per_source = TRUE,          
  n_each = 10,               
  add_source_prefix = FALSE,  
  plot_title = "GSEA - GO:BP",
  seed = 123,
  p_label = 0.05,
  palette = c("#00468BFF", "#ED0000FF"),
  # gseGO controls:
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 1
)


x <- go$gsea$GO_BP$res@result


ggsave("LumA All Patients_GO_BP_GSEA.svg", go$plot, height = 8, width = 10, dpi = 300)


chosen_pos <- c("GO:0002181","GO:0042254","GO:0008380","GO:0042775","GO:0007018", "GO:0006364","GO:0034244", "GO:0042073")
chosen_neg <- c("GO:0002250","GO:0006954","GO:0002758","GO:0042113","GO:0098609",
                "GO:0045087","GO:0050819","GO:0030036","GO:0008544","GO:0002687", "GO:0030855", "GO:0010934")

nes_tbl <- df_all %>%
  dplyr::select(ID, NES, p.adjust,
                Description = dplyr::any_of("Description")) %>%
  mutate(
    Label = if (!all(is.na(Description))) Description else ID,
    Label = gsub("_", " ", Label),
    Label = stringr::str_to_title(Label),
    Label = ifelse(p.adjust < 0.05, paste0(Label, " *"), Label)
  ) %>%
  dplyr::select(-Description)
nes_tbl$ID

combined <- bind_rows(
  nes_tbl %>% arrange(desc(NES)) %>% filter(ID %in% chosen_pos),
  nes_tbl %>% arrange(NES)       %>% filter(ID %in% chosen_neg)
) %>% distinct(ID, .keep_all = TRUE)

# --- helper to prettify labels ---
pretty_go <- function(x){
  x %>%
    str_to_title() %>%
    str_replace_all("\\bDna\\b","DNA") %>%
    str_replace_all("\\bRna\\b","RNA") %>%
    str_replace_all("\\bAtp\\b","ATP") %>%
    str_replace_all("\\bNadh\\b","NADH") %>%
    str_replace_all("\\bMirna\\b","miRNA") %>%
    str_replace_all("\\bRrna\\b","rRNA") %>%
    str_replace_all("\\bB Cell\\b","B cell") %>%    
    str_replace_all("\\bT Cell\\b","T cell") %>%
    
    str_replace_all("(.{1,45})(\\s|$)", "\\1\n") %>%
    str_trim()
}
palette = c("#00468BFF", "#ED0000FF")
combined$Label <- pretty_go(combined$Label)

p <- ggplot(combined, aes(x = reorder(Label, NES), y = NES)) +
  geom_col(aes(fill = p.adjust), width = 0.8) +
  scale_fill_gradient(low = "#00468BFF", high = "#ED0000FF", guide = "colorbar") +
  coord_flip() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
  labs(x = "Pathway", y = "Normalized Enrichment Score", fill = "Adjusted p-value") +
  theme_pubr() +
  theme(panel.spacing = grid::unit(0.1, "lines"),
        legend.position = "right",
        
        axis.text.y     = element_text(size = 9, lineheight = 0.8),  # pathway labels
        axis.text.x     = element_text(size = 9))

ggsave("LumA All Patients_Selected_GO_BP_GSEA.svg", p, height = 8, width = 10, dpi = 300)


net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)
deg <- table_fit %>%
  dplyr::select(Symbol, t) %>% 
  dplyr::filter(!is.na(t))

deg$Symbol <- NULL
deg$t <- as.numeric(deg$t)
deg <-  as.matrix(deg)


contrast_acts <- decoupleR::run_ulm(mat = deg[, 't', drop = FALSE], 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor='mor', 
                                    minsize = 5)



contrast_acts <- contrast_acts %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

sig_contrast_acs <- contrast_acts %>% filter(p_adj < 0.05)

tfs <- sig_contrast_acs$source
emt_tfs <- c("SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2")
f_contrast_acts <- contrast_acts %>%
  filter(source %in% tfs)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = f_contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("TFs")

p


ggsave( "TF_Activity_LumA_Patients.svg",p, height = 5, width = 7)

##Pathway Activity Inference
net <- decoupleR::get_progeny(organism = 'human', 
                              top = 500)


sample_acts <- decoupleR::run_mlm(mat = deg[, 't', drop = FALSE], 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'weight', 
                                  minsize = 5)
sample_acts

contrast_acts <- sample_acts %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

sig_contrast_acs <- contrast_acts %>% filter(p_adj < 0.05)

tfs <- sig_contrast_acs$source
emt_tfs <- c("SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2")
f_contrast_acts <- contrast_acts #%>%
#  filter(source %in% tfs)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = f_contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("TFs")

p


ggsave( "All_Pathway_Activity_LumA_Patients.svg",p, height = 5, width = 7)
