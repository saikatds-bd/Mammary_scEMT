library(edgeR)
library(decoupleR)
library(msigdbr)
library(dplyr)
library(tibble)
library(fgsea)
library(Seurat)
library(presto)
library(scCustomize)
library(SCpubr)
library(ggpubr)
library(ggplot2)

set.seed(123)

#Plotting Functions
source("Helper/MDS_plot.R")
source("Helper/Generate_barplot.R")
source("Helper/enrichment_plot_dotplot.R")
source("Helper/Volcano_function.R")


options(future.globals.maxSize = 16 * 1024^3)

load("Cell_lines_EMT_Project_Colors_New.RData")

load(file="EMT_model_with_CCA.RData")

gc()

DefaultAssay(seurat) <- "RNA"

Idents(seurat) <- "SCT_snn_res.0.2"


#Subset for the pure clusters

clusters <- c("3","2")
emt <- subset(seurat, subset=SCT_snn_res.0.2 %in% clusters)

emt$SCT_snn_res.0.2 <- droplevels(emt$SCT_snn_res.0.2)

valid_cells_ids <- emt@meta.data %>%
  group_by(SCT_snn_res.0.2, Sample) %>%
  tally() %>%
  filter(n >= 50) %>%
  inner_join(emt@meta.data, by = c("SCT_snn_res.0.2", "Sample")) %>%
  pull(Cell_ID)

emt <- subset(emt, cells = valid_cells_ids)

y <- Seurat2PB(emt, sample = "Sample", cluster = "SCT_snn_res.0.2")

y$samples$Cell_line <- factor(ifelse(grepl("D492", y$samples$sample), 
                                     "D492", 
                                     ifelse(grepl("HMLE", y$samples$sample), 
                                            "HMLE",ifelse(grepl("MCF10A", y$samples$sample), "MCF10A", "Other"))))

y$samples$State <- ifelse(y$samples$cluster== "3", "Epithelial", "Mesenchymal")
State <- as.factor(y$samples$State)

y$samples$Sample_Names <- rownames(y$samples)


################Remove lowly expressed genes#############
keep.genes <- filterByExpr(y, group=y$samples$State)
table(keep.genes)
library(limma)
y <- y[keep.genes, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y, method = "TMM")
lcpm <-edgeR::cpm(y, log = T)
y$samples$Group <- y$sample$State
y$samples$Sample_Names <- rownames(y$samples)

p1 <- plotMDS_ggpubr(
  lcpm = lcpm, 
  sample_data = y$samples, 
  color_var = "State", # Use 'Condition' column for colors
  label_var = "Sample_Names", # Use 'Replicate' column for labels
  palette = "Set2",         # Change color palette
  title = "Original MDS Plot"
)

p1
ggsave("EMT_Original_MDS.svg", plot = p1,height = 5, width = 6, dpi = 300)
ggsave("EMT_Original_MDS_Cell_lines.svg", plot = p1,height = 5, width = 6, dpi = 300)

corrected_lcpm <- removeBatchEffect(lcpm, batch = y$samples$Cell_line)
p2 <- plotMDS_ggpubr(
  lcpm = corrected_lcpm, 
  sample_data = y$samples, 
  color_var = "State", # Use 'Condition' column for colors
  label_var = "Sample_Names", # Use 'Replicate' column for labels
  palette = "Set2",          # Change color palette
  title = "Corrected MDS Plot - Cell Line"
)
print(p2)
ggsave("EMT_Original_Corrected_MDS_Cell_lines.svg", plot = p2,height = 5, width = 6, dpi = 300)

pdf("Reasons for excluding HMLE_MES_cluster3.pdf", height = 7, width = 8)
print(p1)
print(p2)
dev.off()

#Removing HMLE MES Cluster 3
y <- Seurat2PB(emt, sample = "Sample", cluster = "SCT_snn_res.0.2")

y$samples$Cell_line <- factor(ifelse(grepl("D492", y$samples$sample), 
                                     "D492", 
                                     ifelse(grepl("HMLE", y$samples$sample), 
                                            "HMLE",ifelse(grepl("MCF10A", y$samples$sample), "MCF10A", "Other"))))


y$samples$State <- ifelse(y$samples$cluster== "3", "Epithelial", "Mesenchymal")

State <- as.factor(y$samples$State)

y$samples$Sample_Names <- rownames(y$samples)

y <- y[, y$samples$Sample_Names != "HMLE_Mes_cluster3"]

keep.genes <- filterByExpr(y, group=y$samples$State)

library(limma)
y <- y[keep.genes, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y, method = "TMM")

lcpm <-edgeR::cpm(y, log = T)

p1 <- plotMDS_ggpubr(
  lcpm = lcpm, 
  sample_data = y$samples, 
  color_var = "State",
  label_var = "Sample_Names", 
  palette = "Set2",         
  title = "Original MDS Plot"
)

p1

corrected_lcpm <- removeBatchEffect(lcpm, batch = y$samples$Cell_line)
p2 <- plotMDS_ggpubr(
  lcpm = corrected_lcpm, 
  sample_data = y$samples, 
  color_var = "State", # Use 'Condition' column for colors
  label_var = "Sample_Names", # Use 'Replicate' column for labels
  palette = "Set2",         # Change color palette
  title = "Corrected MDS Plot - Replicate"
)
print(p2)


pdf("After excluding HMLE_MES_cluster3.pdf", height = 7, width = 8)
print(p1)
print(p2)
dev.off()


state <- as.factor(y$samples$State)

cell_line <- as.factor(y$samples$Cell_line)

design <- model.matrix(~0 + state + cell_line, y$samples)

colnames(design) <- gsub("state","", colnames(design))
design

#####Limma#######

vfit <- voomLmFit(y, design, plot=TRUE, sample.weights = T)
sample_weights <- vfit[["targets"]]
ggbarplot(sample_weights, x="Sample_Names", y="sample.weight", fill = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

contr <- makeContrasts(Mesenchymal - Epithelial, levels = colnames(coef(vfit)))
contr

cfit <- contrasts.fit(vfit, contrasts = contr)
tfit <- treat(cfit, robust = T, fc = 1.2)
dt <- decideTests(tfit)
summary(dt)

mes_vs_ep <- topTreat(tfit, coef="Mesenchymal - Epithelial", n=Inf)

mes_vs_ep$Group <- ifelse(mes_vs_ep$logFC > 0, "Mesenchymal", "Epithelial")

sig <- mes_vs_ep %>% 
  filter(adj.P.Val < 0.05) %>% 
  arrange(adj.P.Val)


plot_volcano(data=mes_vs_ep, plot_title = "EMT Model Cell Lines - Mesenchymal vs Epithelial Cluster", 
             output_file = "Volcano EMT Model Cell Lines - Mesenchymal vs Epithelial Cluster.svg", 
             height = 4, width = 5)


save(sig, file = "Filtered_Samples_DE_Results_treat.RData")

save(mes_vs_ep,file= "DE_EMT_Filtered.RData")

















dput(unique(combined_emt$Source))
names(combined_emt)
x <- intersect(sig$gene, all_sig$gene)
y <- setdiff( all_sig$gene,sig$gene)
z <- all_sig %>% filter(gene %in% y)
y
all_sig <- sig

x <- mes_vs_ep %>% filter(gene %in% emt_hallmark)
x$Group <- ifelse(x$logFC > 0, "Mesenchymal", "Epithelial")
names(x)

table(x$Group)
hist(mes_vs_ep$P.Val)
#write.csv2(mes_vs_ep, "Combined_top_treat.csv")
#BiocManager::install("EnhancedVolcano")
##EnhancedVolcano
sig <- mes_vs_ep %>% 
  filter(adj.P.Val < 0.05) %>% 
  arrange(adj.P.Val)
hist(sig$adj.P.Val)
sig$Group <- ifelse(sig$logFC > 0,  "Mesenchymal", "Epithelial")
table(sig$Group)
write.csv2(sig, "Cell_line_Pseudo_EMT.csv")

#load("C:/Users/ssa214/UiT Office 365/O365-PhD_Saikat - General/Single Cell RNA-seq/EMT Project/20241022/EMT_Stingent.RData")
#sig <- voom_trade_05
#de <- x %>% filter(gene %in% )
#mt2g <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)
m_t2tf <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% 
  dplyr::select(gs_name,gene_symbol)


background_genes <- rownames(seurat)
ep <- sig %>% filter(logFC < 0) %>% pull(gene)
mes <- sig %>% filter(logFC > 0) %>% pull(gene)

#down_D492 <- D492 %>% filter(avg_log2FC < -1.5)


#BiocManager::install("clusterProfiler")
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

em_up <- enricher(mes, TERM2GENE = m_t2g, universe = background_genes)
em_down <- enricher(ep, TERM2GENE = m_t2g, universe = background_genes)

# Generate and store the dot plot
up_em_plots <- dotplot(em_up) + ggtitle("EMT Up Regulated Genes")
down_em_plots <- dotplot(em_down) + ggtitle("Down Regulated Genes")

# Generate and store the dot plot
up_em_plots <- dotplot(em_up) + ggtitle("EMT Up Regulated Genes")
down_em_plots <- dotplot(em_down) + ggtitle("Down Regulated Genes")

em_up_tf <- enricher(mes, TERM2GENE = m_t2tf, universe = background_genes)
em_down_tf <- enricher(ep, TERM2GENE = m_t2tf, universe = background_genes)

up_em_plots_tf <- dotplot(em_up_tf) + ggtitle("EMT Up Regulated Genes: C3:TF_GTRD")
down_em_plots_tf <- dotplot(em_down_tf) + ggtitle("Down Regulated Genes: C3:TF_GTRD")
#write.csv2(voom_trade_05, "EMT_Stingent_List.csv")
library(org.Hs.eg.db)
genes_map <- clusterProfiler::bitr(geneID = sig$gene, #Specify the column that has the gene symbol
                                     fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL", "UNIPROT"), #You will get both the associated entrez ids and ensembl ids
                                     OrgDb="org.Hs.eg.db") #it is a widely used datbase. Some people use biomart, but that one is pretty slow
genes_map <- distinct(genes_map, SYMBOL, .keep_all = T)
sig <- merge(sig, genes_map, by.x = "gene", by.y = "SYMBOL")

universe <- clusterProfiler::bitr(geneID = background_genes, #Specify the column that has the gene symbol
                                  fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL", "UNIPROT"), #You will get both the associated entrez ids and ensembl ids
                                  OrgDb="org.Hs.eg.db") #it is a widely used datbase. Some people use biomart, but that one is pretty slow
write.csv2(universe, "Background_Genes.csv")
ep_go <- sig %>% filter(logFC < 0) %>% pull(ENTREZID)
mes_go <- sig %>% filter(logFC > 0) %>% pull(ENTREZID)

ep_go_res <- enrichGO(ep_go, OrgDb = org.Hs.eg.db, ont = "BP",
                      universe = universe$ENTREZID)
mes_go_res <- enrichGO(mes_go, OrgDb = org.Hs.eg.db, ont = "BP",
                       universe = universe$ENTREZID)


go_up_em_plots <- dotplot(mes_go_res) + ggtitle("EMT Up Regulated Genes - GO")
go_down_em_plots <- dotplot(ep_go_res) + ggtitle("EMT Down Regulated Genes - GO")

ep_wiki_res <- enrichWP(ep_go, organism = "Homo sapiens",
                      universe = universe$ENTREZID)
mes_wiki_res <- enrichWP(mes_go, organism = "Homo sapiens",
                       universe = universe$ENTREZID)


wiki_up_em_plots <- dotplot(mes_wiki_res) + ggtitle("EMT Up Regulated Genes - Wikipathway")
wiki_down_em_plots <- dotplot(ep_wiki_res) + ggtitle("EMT Down Regulated Genes - Wikipathway")

ep_kegg_res <- enrichKEGG(ep_go, organism = "hsa",
                        universe = universe$ENTREZID,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1)
mes_kegg_res <- enrichKEGG(mes_go, organism = "hsa",
                         universe = universe$ENTREZID)

kegg_up_em_plots <- dotplot(mes_kegg_res) + ggtitle("EMT Up Regulated Genes - KEGG")
kegg_down_em_plots <- dotplot(ep_kegg_res) + ggtitle("EMT Down Regulated Genes - KEGG")


pdf("ORA on Stingent EMT Up and Down.pdf", height = 4, width = 7)
print(up_em_plots)
print(down_em_plots)
print(go_up_em_plots)
print(go_down_em_plots)
print(wiki_up_em_plots)
#print(wiki_down_em_plots)
print(kegg_up_em_plots)
print(kegg_down_em_plots)
dev.off()



#Visualization of lowly expressed genes
library(reshape2)
cpm <- cpm(y, log = F)
corrected_expression_matrix <- removeBatchEffect(cpm, batch = y$samples$Cell_line)

genes_to_check <- sig %>% filter(AveExpr < 0) %>% pull(gene)
emt_genes <- c("CDH1", "CDH2", "EPCAM", "SNAI1", "SNAI2", "ZEB1", "ZEB2")

expr <- corrected_expression_matrix[genes_to_check,]
expr <- cpm[genes_to_check,]

expr <- corrected_expression_matrix[emt_genes,]
expr <- cpm[emt_genes,]
expr_long <- melt(expr)

# Rename columns for better readability
colnames(expr_long) <- c("Gene", "Sample", "Expression")

# Create a separate bar plot for each gene
unique_genes <- unique(expr_long$Gene)

# Loop through each gene and create a bar plot
for (gene in unique_genes) {
  gene_data <- expr_long %>% filter(Gene == gene)
  
  # Create the bar plot for the current gene
  p <- ggplot(gene_data, aes(x = Sample, y = Expression)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    labs(title = paste("Expression of", gene),
         x = "Sample",
         y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample names for better readability
  
  # Print the plot
  print(p)
}

#Checking overlap with FitGAM
library(tradeSeq)
load("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Single Cell RNA-seq/EMT Project/20241002/Fitgam_cell_line_all_lineages_monocle3.RData")
startRes <- startVsEndTest(sce, lineages=TRUE, l2fc = 0.5)
trade_emt <- startRes[,c(10:12,27)]
names(trade_emt) <- c( "waldStat","df", "pvalue","logFC")
trade_emt$padj <- p.adjust(trade_emt$pvalue, method = "BH")
all_trade_emt <- trade_emt %>% filter(padj < 0.05)
all_trade_emt$gene <- rownames(all_trade_emt)
all_common <- sig %>% filter(gene %in% all_trade_emt$gene)


load("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Single Cell RNA-seq/EMT Project/20241002/Fitgam_cell_line_monocle3.RData")
startRes <- startVsEndTest(sce, lineages=TRUE, l2fc = 0.5)
trade_emt <- startRes
trade_emt$padj <- p.adjust(trade_emt$pvalue, method = "BH")
trade_emt <- trade_emt %>% filter(padj < 0.05)

all_trade_emt$gene <- rownames(all_trade_emt)
sig$UCell <- ifelse(sig$Group == "Epithelial", paste0(sig$gene, "-"), paste0(sig$gene, "+"))

voom_trade_05 <- sig %>% filter(gene %in% all_trade_emt$gene)
length(setdiff(sig$gene, filtered_emt$gene))
voom_05 <- sig
save(voom_05, voom_trade_05, file="EMT_Stingent.RData")
