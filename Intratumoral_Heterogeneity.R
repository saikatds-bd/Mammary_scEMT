library(Seurat)
library(qs2)
source("~/Helpers.R")
source("MDS_plot.R")
source("Generate_barplot.R")
source("perform_enrichment_function.R")
source("Updated_Volcano_function.R")
source("perform_GSEA.R")
source("hallmark_gsea_plot_function.R")
#rm(gsea_barplot_multi)
source("EMT_Hallmark_Helper.R")

load("Annotated_Pang.RData")

seurat <- subset(seurat, subset = Refined_Annotation == "Cancer" & Pam50_Subtype = "LumA")
pang_meta <- seurat@meta.data

library(dplyr)
library(tibble)

# choose the epithelial score column
score_col <- "Epithelial_AMS"

# sanity check
stopifnot(score_col %in% colnames(seurat@meta.data))

stopifnot("Patient_ID" %in% colnames(seurat@meta.data))

subtype_info <- pang_meta %>% dplyr::select(Patient_ID, PAM50_Subtype) %>% distinct(Patient_ID, .keep_all = T)

# summarize per patient
epi_summary <- pang_meta %>%
  group_by(Patient_ID) %>%
  summarise(
    n_cells = n(),
    epi_min = min(.data[[score_col]], na.rm = TRUE),
    epi_max = max(.data[[score_col]], na.rm = TRUE),
    epi_median_Q1 = quantile(.data[[score_col]], 0.25, na.rm = TRUE),
    epi_median_Q4 = quantile(.data[[score_col]], 0.75, na.rm = TRUE),
    epi_delta = epi_median_Q4 - epi_median_Q1,  # N_epi
    epi_sd = sd(.data[[score_col]], na.rm = TRUE),
    epi_iqr = IQR(.data[[score_col]], na.rm = TRUE)
  ) %>%
  mutate(
    pass_range = epi_delta >= 0.10,
    pass_sd    = epi_sd >= 0.04 | epi_iqr >= 0.06,
    include_for_erosion = pass_range & pass_sd
  ) %>%
  arrange(desc(include_for_erosion), Patient_ID)

# show summary table
epi_summary
epi_summary <- merge(subtype_info,epi_summary,  by = "Patient_ID")
writexl::write_xlsx(epi_summary, "Criteria_for_Inclusion.xlsx")


patients <- c("ER0025", "ER0032", "ER0042", "ER0114", "ER0125", "ER0151", 
              "HER20308", "CA1", "CA2", "CA4", "CA5", "CID4463", "CID4471", 
              "CID4530N", "CID4535", "CID3941", "CID3948", "CID4067", "CID4290A")


marker_list <- list()

for (pid in patients) {
  message("Processing patient: ", pid)
  
  # 1. Subset to patient
  obj <- subset(seurat, subset = Patient_ID == pid)
  
  # 2. Rank and create quartiles
  vals <- obj$Common_Filtered_Epithelial_AMS
  rnk  <- rank(vals, ties.method = "random", na.last = "keep")
  brks <- quantile(rnk, probs = c(0, 0.25, 0.50, 0.75, 1.0), na.rm = TRUE, names = FALSE)
  
  obj$EpiQuartile <- cut(
    rnk,
    breaks = brks,
    include.lowest = TRUE,
    labels = c("Q1_lowest", "Q2", "Q3", "Q4_highest")
  )
  obj$EpiQuartile <- factor(obj$EpiQuartile, 
                            levels = c("Q1_lowest", "Q2", "Q3", "Q4_highest"))
  
  # 3. Run Presto for Q1 vs Q4
  markers <- RunPresto(obj, ident.1 = "Q1_lowest", ident.2 = "Q4_highest", group.by = "EpiQuartile")
  markers$Symbol <- rownames(markers)
  markers$Patient_ID <- pid
  
  marker_list[[pid]] <- markers
  
  # 4. Plot AUC distribution
  quartiles <- quantile(obj$Common_Filtered_Epithelial_AUC, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  p <- ggdensity(
    data = obj@meta.data,
    x = "Epithelial_AMS",
    bins = 50,
    fill = "#00468B",
    color = "white"
  ) +
    geom_vline(xintercept = quartiles, linetype = "dashed", color = "red") +
    annotate("text", x = quartiles, y = 20, label = c("Q1", "Q2", "Q3"), angle = 90, vjust = -0.5) +
    labs(x = "Epithelial AMS score", y = "Cell count", 
         title = paste0(pid, " epithelial score distribution")) +
    theme_pubr()
  
  ggsave(
    filename = paste0(pid, "_Epi_Distribution.pdf"),
    plot = p,
    height = 3, width = 3
  )
  
  gc()
}

# 5. Combine all into one dataframe
all_markers <- bind_rows(marker_list)
qs_save(all_markers, "AllLumA_Patients_Q1vsQ4_Markers.qs2")

###GSEA
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(org.Hs.eg.db)

conflicted::conflict_prefer("summarise", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(base::setdiff)



m_t2g_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)
m_t2go <- msigdbr(species = "Homo sapiens", collection = "C5",
                  subcollection = "GO:BP") %>% 
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name,gene_symbol)
x <- msigdbr_collections()

out <- run_hallmark_per_cluster(
  luma = all_markers,
  m_t2g_hallmark = m_t2g_hallmark,  # Hallmark TERM2GENE (msigdbr)
  use = c("H"),
  cluster_col = "Patient_ID",
  gene_col = "Symbol",
  lfc_col  = "avg_log2FC",
  pval_col = "p_val",
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

out_go <- run_hallmark_per_cluster(
  luma = all_markers,
  use = "GO:BP",
  m_t2g_hallmark = m_t2g_hallmark,
  cluster_col = "Patient_ID",
  gene_col = "Symbol",
  lfc_col  = "avg_log2FC",
  pval_col = "p_val"
)


##Doing ORA
get_common_regulated <- function(df,
                                 patient_col = "Patient_ID",
                                 gene_col    = "Symbol",
                                 lfc_col     = "avg_log2FC",
                                 padj_col    = "p_val_adj",
                                 padj_thr    = 0.05,
                                 lfc_thr     = 0,
                                 min_patients = NULL  
){
  stopifnot(all(c(patient_col,gene_col,lfc_col,padj_col) %in% names(df)))
  if (is.null(min_patients)) min_patients <- length(unique(df[[patient_col]]))
  
  
  core <- df |>
    dplyr::filter(!is.na(.data[[gene_col]]), !is.na(.data[[lfc_col]]), !is.na(.data[[padj_col]])) |>
    dplyr::group_by(.data[[patient_col]], .data[[gene_col]]) |>
    dplyr::summarise(
      padj = min(.data[[padj_col]], na.rm = TRUE),
      lfc  = mean(.data[[lfc_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
 
  up_sets <- core |>
    dplyr::filter(padj < padj_thr, lfc >  lfc_thr) |>
    dplyr::group_by(.data[[patient_col]]) |>
    dplyr::summarise(genes = list(unique(.data[[gene_col]])), .groups="drop")
  
  down_sets <- core |>
    dplyr::filter(padj < padj_thr, lfc < -lfc_thr) |>
    dplyr::group_by(.data[[patient_col]]) |>
    dplyr::summarise(genes = list(unique(.data[[gene_col]])), .groups="drop")
  
 
  tally_common <- function(sets_tbl){
    if (nrow(sets_tbl) == 0) return(list(counts = dplyr::tibble(Gene=character(), n_patients=integer()),
                                         common = dplyr::tibble(Gene=character())))
    all_genes <- unlist(sets_tbl$genes, use.names = FALSE)
    counts <- dplyr::tibble(Gene = all_genes) |>
      dplyr::count(Gene, name = "n_patients") |>
      dplyr::arrange(dplyr::desc(n_patients), Gene)
    
    list(
      counts = counts,
      common = counts |>
        dplyr::filter(n_patients >= min_patients) |>
        dplyr::select(Gene)
    )
  }
  
  up_stats   <- tally_common(up_sets)
  down_stats <- tally_common(down_sets)
  
  
  list(
    params = list(padj_thr=padj_thr, lfc_thr=lfc_thr, min_patients=min_patients),
    up_by_patient   = up_sets,            
    down_by_patient = down_sets,
    up_counts       = up_stats$counts,    
    down_counts     = down_stats$counts,
    up_common       = up_stats$common,    
    down_common     = down_stats$common
  )
}


all_markers <- qs_read("AllLumA_Patients_Q1vsQ4_Markers.qs2")

luma_included <- epi_summary %>% filter(include_for_erosion == "TRUE" & PAM50_Subtype == "LumA") %>% pull(Patient_ID)

luma_markers <- all_markers %>% filter(Patient_ID %in% luma_included)


# Taking genes only present in above 75% of patients
ceil_70 <- ceiling(0.75 * length(unique(luma_markers$Patient_ID)))
res70 <- get_common_regulated(luma_markers, min_patients = ceil_70,
                              lfc_thr = 1)


common_down <- luma_markers %>% filter(Symbol %in% res70[["down_common"]][["Gene"]])
common_down$Presence_in_EMT_Signature <- ifelse(common_down$Symbol %in%
                                                  New_EMT_Symbol$Common_Filtered_Epithelial,
                                                "Present", "Absent")

dataframe <- list("All DE - LumA" = luma_markers,
                  "Common Downregulated" = common_down)

library(org.Hs.eg.db)
library(msigdbr)
writexl::write_xlsx(dataframe,"ST7_Common_LUMA_Downregulated.xlsx")



down_enrichment <- perform_enrichment_analysis(gene_list = unique(common_down$Symbol), 
                                               universe_genes = unique(all_markers$Symbol),
                                               #species = "Homo sapiens",
                                               subset_name = "Common Downregulated in 75% LumA Patients",
                                               output_prefix = "Common Downregulated in LumA_ORA",
                                               height = 4, width = 6)


