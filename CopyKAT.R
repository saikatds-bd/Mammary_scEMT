library(Seurat)
library(copykat)
library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
load("55pt_Pang_Processed.qs2")  

out_dir <- "~/copyKAT"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

DefaultAssay(seurat) <- "RNA"
seurat <- DietSeurat(seurat, assays = "RNA")
gc()

exclude_as_reference <- c("Cancer","Epithelial","CAF","Fibroblast","Stromal","Smooth muscle cells")
normal_cell_types    <- c("T-Cells","B-Cells","Myeloid","Plasmablast","NK","Mast")

eq_safe <- !is.na(seurat$Author_Annotation) & !is.na(seurat$Final_Annotation) &
  seurat$Author_Annotation == seurat$Final_Annotation

non_malig <- !(seurat$Author_Annotation %in% exclude_as_reference)

seurat$Agreement <- ifelse(eq_safe & non_malig, "Yes", "No")

is_norm_author <- seurat$Author_Annotation %in% normal_cell_types
is_norm_final  <- seurat$Final_Annotation  %in% normal_cell_types
seurat$CopyKAT_Ref <- ifelse(seurat$Agreement == "Yes" & is_norm_author & is_norm_final, "Yes", "No")

set.seed(123)
ref_idx <- which(seurat$CopyKAT_Ref == "Yes")
ref_by_patient <- split(ref_idx, seurat$Patient_ID[ref_idx])
ref_keep <- unlist(lapply(ref_by_patient, function(idx){
  if (length(idx) > 2000) sample(idx, 2000) else idx
}))
seurat$CopyKAT_Ref_final <- ifelse(seq_along(colnames(seurat)) %in% ref_keep, "Yes", "No")


run_copykat_per_patient <- function(seu, patient_id,
                                    max_ref_frac = 0.5, max_ref_n = 2000,
                                    ncores = 16, ks_cut = 0.1, win_size = 25, ngene_chr = 5) {
  cells_pt <- colnames(seu)[seu$Patient_ID == patient_id]
  if (length(cells_pt) == 0) stop("No cells for patient: ", patient_id)
  
  # Reference cells in this patient
  ref_pt <- intersect(colnames(seu)[seu$CopyKAT_Ref_final == "Yes"], cells_pt)
  
  # cap reference
  if (length(ref_pt) > max_ref_frac * length(cells_pt)) {
    set.seed(1); ref_pt <- sample(ref_pt, floor(max_ref_frac * length(cells_pt)))
  }
  if (length(ref_pt) > max_ref_n) {
    set.seed(1); ref_pt <- sample(ref_pt, max_ref_n)
  }
  if (length(ref_pt) < 50) {
    message("[", patient_id, "] Reference < 50 (", length(ref_pt), "); using autodetect normals.")
    ref_pt <- NULL
  }
  
  # Get raw counts as read in the github
  raw_counts <- GetAssayData(seu, slot = "counts")[, cells_pt, drop = FALSE]
  
  # Run CopyKAT
  ck <- copykat(
    rawmat          = as.matrix(raw_counts),  
    id.type         = "S",
    norm.cell.names = ref_pt,                 
    ngene.chr       = ngene_chr,
    win.size        = win_size,
    KS.cut          = ks_cut,
    n.cores         = ncores,
    sam.name        = as.character(patient_id)
  )
  
  
  ck_out <- list(
    patient_id = patient_id,
    params = list(id.type="S", ngene.chr=ngene_chr, win.size=win_size, KS.cut=ks_cut, n.cores=ncores,
                  n_cells = length(cells_pt), n_ref = length(ref_pt)),
    prediction = ck$prediction,     
    CNAmat     = ck$CNAmat          
  )
  ck_out
}

patients <- sort(unique(seurat$Patient_ID))


for (pid in patients) {
  msg_prefix <- paste0("[", pid, "] ")
  out_file <- file.path(out_dir, paste0(pid, "_CopyKAT.RData"))
  if (file.exists(out_file)) {
    message(msg_prefix, "exists; skipping.")
    next
  }
  message(">> ", msg_prefix, "running CopyKAT on ", sum(seurat$Patient_ID == pid), " cells")
  
  
  ck_out <- NULL
  try({
    ck_out <- run_copykat_per_patient(seurat, pid, ncores = 8)
    save(ck_out, file = out_file, compress = TRUE)
    rm(ck_out); invisible(gc())
    notify(paste0("[OK] CopyKAT ", pid, " saved ??? ", basename(out_file)))
  }, silent = TRUE)
  
  if (!file.exists(out_file)) {
    notify(paste0("[FAIL] CopyKAT ", pid, " (see console)"))
  }
  
  rm(list = ls(pattern = "^ck_out$"))
  invisible(gc())
}


pred_files <- list.files(out_dir, pattern = "_CopyKAT\\.RData$", full.names = TRUE)

extract_pred <- function(f) {
  env <- new.env(parent = emptyenv())
  load(f, envir = env)  
  if (!exists("ck_out", envir = env)) return(NULL)
  pred <- env$ck_out$prediction
  if (is.null(pred) || nrow(pred) == 0) return(NULL)
  
  # standrdizing column names
  cnl <- tolower(colnames(pred))
  cell_col <- colnames(pred)[which.max(cnl %in% c("cell.name","cell","barcode","cells"))]
  pred_col <- colnames(pred)[which.max(cnl %in% c("copykat.pred","pred"))]
  clst_col <- colnames(pred)[which.max(cnl %in% c("cluster","copykat.cluster"))]
  
  tibble::tibble(
    Cell       = pred[[cell_col]],
    CopyKAT    = as.character(pred[[pred_col]]),
    CK_Cluster = if (!is.null(clst_col) && nzchar(clst_col)) as.character(pred[[clst_col]]) else NA_character_,
    Patient_ID = env$ck_out$patient_id
  )
}

pred_list <- lapply(pred_files, extract_pred)
pred_df <- bind_rows(pred_list)

# Save merged predictions
save(pred_df, file = file.path(out_dir, "CopyKAT_predictions_merged.RData"), compress = TRUE)
readr::write_csv(pred_df, file.path(out_dir, "CopyKAT_predictions_merged.csv"))

### After CopyKAT

library(dplyr)
library(readr) 
library(stringr) 
library(tidyr)

ck_dir <- "~/Pang/CopyKAT_Individual" 

extract_pred <- function(f) {
  env <- new.env(parent = emptyenv())
  load(f, envir = env)
  ck <- env$ck_out
  if (is.null(ck) || is.null(ck$prediction)) return(NULL)
  
  pred <- ck$prediction
  cnl <- tolower(colnames(pred))
  cell_col <- colnames(pred)[which.max(cnl %in% c("cell.name","cell","barcode","cells"))]
  pred_col <- colnames(pred)[which.max(cnl %in% c("copykat.pred","pred"))]
  
  tibble(
    Cell       = pred[[cell_col]],
    CopyKAT_new= as.character(pred[[pred_col]]),
    Patient_ID = ck$patient_id
  )
}

files <- list.files(pattern = "_CopyKAT\\.RData$")
pred_new <- bind_rows(lapply(files, extract_pred))

# --- Summarise per patient ---

print(summary_new, n = Inf)
metadata <- merge(metadata, pred_new, by = "Cell_ID")

library(ggsankey)
metadata$scATOMIC_Annotation <- metadata$Final_Annotation

tumor_like <- c("Cancer","Epithelial")
cnv_levels <- c("Aneuploid","Diploid","Unassigned")

md <- metadata %>%
  as_tibble(rownames = "Cell") %>%
  mutate(
    Author_Annotation = as.character(Author_Annotation),
    scATOMIC_Annotation  = as.character(scATOMIC_Annotation),
    CNV_Status        = factor(CNV_Status, levels = cnv_levels)
  )


md$Refined_Annotation <- md$Final_Annotation

lock_idx <- with(md, Author_Annotation == scATOMIC_Annotation & scATOMIC_Annotation %in% tumor_like)

adj_idx <- which(!lock_idx & md$scATOMIC_Annotation %in% tumor_like)

i_ep_to_ca <- adj_idx[md$scATOMIC_Annotation[adj_idx] == "Epithelial" & md$CNV_Status[adj_idx] == "Aneuploid"]
md$Refined_Annotation[i_ep_to_ca] <- "Cancer"

i_ca_to_ep <- adj_idx[md$scATOMIC_Annotation[adj_idx] == "Cancer" & md$CNV_Status[adj_idx] == "Diploid"]
md$Refined_Annotation[i_ca_to_ep] <- "Epithelial"


merged <- md

rownames(merged) <- merged$Cell_ID
merged <- merged[Cells(seurat),]
all(rownames(merged) == rownames(seurat@meta.data))

seurat@meta.data <- merged
save(seurat, file = "Annotated_Pang.RData")

cancer <- subset(seurat, subset = Refined_Annotation == "Cancer")

library(genefu)
perform_pam50_classification <- function(seurat_object, symbol_to_gene_id_mapping) {
  
  object_avg <- AggregateExpression(seurat_object, return.seurat = TRUE, group.by = "Patient_ID")
  gc()
  expr_matrix <- object_avg[["RNA"]]$data
  t_expr_matrix <- t(expr_matrix)
  gc()
  
  gene_annotations <- data.frame(
    probe = colnames(t_expr_matrix),
    Gene.Symbol = colnames(t_expr_matrix), 
    EntrezGene.ID = symbol_to_gene_id_mapping[match(colnames(t_expr_matrix), symbol_to_gene_id_mapping$symbol), 1]
  )
  
  
  data(pam50.robust)
  object_predictions <- molecular.subtyping(sbt.model = "pam50", data = t_expr_matrix, annot = gene_annotations, do.mapping = FALSE)
  gc()
  object_subtypes <- as.data.frame(object_predictions$subtype)
  object_subtypes$Patient_ID <- rownames(object_subtypes)
  object_subtypes$Patient_ID <- gsub("-", "_", object_subtypes$Patient_ID)
  names(object_subtypes) <- c("Genefu_PAM50", "Patient_ID")
  
  return(object_subtypes)
}

library(org.Hs.eg.db)

s2g <- toTable(org.Hs.egSYMBOL)
subtypes <- perform_pam50_classification(cancer, s2g)
own_subtype <- read.csv(file="All_Subtype.csv", as.is = T, check.names = F)

merged_subtypes <- merge(subtypes, own_subtype, by= "Patient_ID")
merged_subtypes <- merged_subtypes %>% distinct(Patient_ID, .keep_all = T)

library(scCustomize)
seurat <- Add_Sample_Meta(seurat, meta_data = merged_subtypes, join_by_seurat = "Patient_ID", join_by_meta = "Patient_ID")
save(seurat, file = "Annotated_Pang.RData")









