# Load necessary libraries
library(readxl)
library(tidyverse)
library(HGNChelper)
library(msigdbr)

# Load EMT DE data
load(file= "DE_EMT_Filtered.RData")

filtered_emt <- read.csv2("Saikat_EMT.csv", as.is = TRUE, check.names = FALSE, row.names = 1)
filtered_emt$Source <- "Saikat et al."

# Standardize gene symbols
saikat_updated <- HGNChelper::checkGeneSymbols(filtered_emt$gene, unmapped.as.na = FALSE, species = "human")
filtered_emt$gene <- saikat_updated$Suggested.Symbol

# Load Hallmark EMT genes
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
hallmark_emt <- m_t2g %>% filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% pull(gene_symbol)
hallmark_updated <- HGNChelper::checkGeneSymbols(hallmark_emt, unmapped.as.na = FALSE, species = "human")
hallmark_emt <- data.frame(gene = hallmark_updated$Suggested.Symbol, Source = "Hallmark", stringsAsFactors = FALSE)

# Load Tan EMT markers
library(stringr)
tan_cell <- read.csv2("Tan_EMT_Cell_line.csv", as.is = T, check.names = F)
tan_cell <- tan_cell[!apply(tan_cell == "", 1, all),]
head(tan_cell)
tan_cell$Number <- str_extract(tan_cell$Symbol, "^[0-9]+")
tan_cell$GeneSymbol <- sub("^[0-9]+", "", tan_cell$Symbol)

tan_cell$Category <- str_replace(tan_cell$Category, "^[0-9]+\\.?[0-9]*", "")
tan_cell$Symbol <- tan_cell$GeneSymbol
tan_cell <- tan_cell[,c(3,1,2)]
tan_cell$Source <- ifelse(tan_cell$Symbol %in% tan_tumor$Symbol, "Tan_Both", "Tan_Cell")

tan_tumor <- read.csv2("Tan_EMT.csv", as.is = T, check.names = F)
tan_tumor$Source <- ifelse(tan_tumor$Symbol %in% tan_cell$Symbol, "Tan_Both", "Tan_Tumor")


tan_emt <- rbind(tan_tumor, tan_cell)
tan_emt <- distinct(tan_emt, Symbol, .keep_all = T)
write.csv2(tan_emt, "Tan_Combined_EMT.csv", row.names = F)

tan_updated <- HGNChelper::checkGeneSymbols(tan_emt$Symbol, unmapped.as.na = FALSE, species = "human")
tan_updated <- merge(tan_emt, tan_updated, by.x="Symbol", by.y = "x")
tan_emt <- data.frame(gene = tan_updated$Suggested.Symbol, Source = tan_updated$Source, stringsAsFactors = FALSE)
tan_emt$Source <- "Tan et al."

# Merge all EMT lists
emt_combined <- bind_rows(filtered_emt %>% select(gene, Source),
                          hallmark_emt,
                          tan_emt)

emt_final <- emt_combined %>%
  group_by(gene) %>%
  summarise(Source = paste(unique(Source), collapse = " + "), .groups = "drop")

all_emt <- merge(emt_final, mes_vs_ep, by.x="gene", by.y="gene", all.x=T)

write.csv2(all_emt, "Unfiltered_Combined_EMT_Signature.csv")
save(all_emt, file="EMT_Signature.RData")