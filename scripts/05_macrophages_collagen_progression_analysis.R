# MPNST Collagen & Cell Cycle Gene Analysis in Cancer-associated Macrophages
# This script performs two parallel analyses:
#   1. Using COL6 high/low status (based on COL6A1/2/3 expression) along with additional marker analyses
#   2. Using TRAP1 high/low status (plus additional marker analyses: invasion, metastasis,
#      angiogenesis, and other collagen genes)
# Differential expression (DEG) is then run for both High vs Low and Low vs High comparisons for each grouping.

library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
library(reshape2)  # For melt function used in correlation heatmap
library(pheatmap)
library(gridExtra)
library(ggrepel)

# Create directories for outputs if they don't exist
if (!dir.exists("pictures")) dir.create("pictures")
if (!dir.exists("statistics")) dir.create("statistics")

# Load the pre-processed Seurat object
data.merged <- readRDS("data_mpnst.rds")

#-----------------------------------------------
# 1. Cell subsetting workflow: MPNST -> Tumor-associated macrophages
#-----------------------------------------------
cat("Subsetting MPNST cells...\n")
mpnst_only <- subset(data.merged, subset = orig.ident == "MPNST")
cat("Total MPNST cells:", ncol(mpnst_only), "\n")

# Preserve cell type annotations from the merged object
Idents(mpnst_only) <- Idents(data.merged)

cat("Subsetting to Tumor-associated macrophages...\n")
tumor_macrophages <- subset(mpnst_only, idents = "Tumor-associated macrophage")
cat("Tumor-associated macrophages:", ncol(tumor_macrophages), "\n")
if(ncol(tumor_macrophages) < 10) {
  stop("Not enough Tumor-associated macrophages for meaningful analysis")
}

# Normalize data and set assay to "RNA"
DefaultAssay(tumor_macrophages) <- "RNA"
tumor_macrophages <- NormalizeData(tumor_macrophages)

#########################################
# Branch 1: Analysis using COL6 high/low status
#########################################

# --- COL6 grouping ---
col6_genes <- c("COL6A1", "COL6A2", "COL6A3")
missing_genes <- col6_genes[!col6_genes %in% rownames(tumor_macrophages)]
if(length(missing_genes) > 0) {
  warning(paste("The following collagen genes are missing:", paste(missing_genes, collapse=", ")))
  col6_genes <- col6_genes[col6_genes %in% rownames(tumor_macrophages)]
  if(length(col6_genes) == 0) { stop("No COL6A1/2/3 genes found in the dataset!") }
}

categorize_by_expression <- function(gene) {
  if(gene %in% rownames(tumor_macrophages)) {
    exp_values <- GetAssayData(tumor_macrophages, layer = "data")[gene, ]
    cutoff <- median(exp_values)
    status <- ifelse(exp_values > cutoff, paste0(gene, "-high"), paste0(gene, "-low"))
    cat(gene, "expression cutoff (median):", round(cutoff, 4), "\n")
    print(table(status))
    return(list(status = status, cutoff = cutoff))
  } else {
    return(NULL)
  }
}

# Record COL6 gene status in metadata
col6_status <- list()
for(gene in col6_genes) {
  col6_status[[gene]] <- categorize_by_expression(gene)
  tumor_macrophages[[paste0(gene, "_status")]] <- col6_status[[gene]]$status
}

# Create combined COL6 score (average z-score)
if(length(col6_genes) > 1) {
  z_scores <- matrix(0, nrow = length(col6_genes), ncol = ncol(tumor_macrophages))
  rownames(z_scores) <- col6_genes
  for(i in 1:length(col6_genes)) {
    gene <- col6_genes[i]
    exp_values <- GetAssayData(tumor_macrophages, layer = "data")[gene, ]
    z_scores[i, ] <- scale(exp_values)
  }
  col6_mean_zscore <- colMeans(z_scores)
  col6_status_combined <- ifelse(col6_mean_zscore > 0, "COL6-high", "COL6-low")
  tumor_macrophages$COL6_status <- col6_status_combined
  cat("Combined COL6 status:\n")
  print(table(tumor_macrophages$COL6_status))
}

# Set identity to COL6_status
if("COL6_status" %in% colnames(tumor_macrophages@meta.data)) {
  Idents(tumor_macrophages) <- "COL6_status"
  current_status <- "COL6_status"
} else {
  Idents(tumor_macrophages) <- paste0(col6_genes[1], "_status")
  current_status <- paste0(col6_genes[1], "_status")
}

# UMAP plot for COL6 grouping
p_status_col6 <- DimPlot(tumor_macrophages, group.by = "COL6_status", 
                         cols = c("COL6-low" = "blue", "COL6-high" = "red")) +
                 ggtitle("Combined COL6 Status Groups") +
                 theme(plot.title = element_text(size = 14, face = "bold"))
ggsave("pictures/COL6_status_groups_tumor_macrophages.pdf", p_status_col6, width = 8, height = 7)

# Cell Cycle Gene Analysis (using same list as before)
list_cell_cycle <- "CDK2,CDK4,CDK6,CDK7,CDKN1A,CDKN1B,STAG1,CDKN1C,CDKN2A,CDKN2B,CDKN2C,CDKN2D,ANAPC10,MAD2L2,STAG2,PTTG2,GADD45G,DBF4,YWHAQ,CHEK1,CHEK2,CREBBP,GADD45A,E2F1,E2F2,E2F3,E2F4,E2F5,EP300,ORC6,ORC3,CDC26,ABL1,ANAPC13,SMC1B,SFN,GSK3B,ANAPC2,ANAPC4,HDAC1,HDAC2,MAD2L1,SMAD2,SMAD3,SMAD4,MCM2,MCM3,MCM4,MCM5,MCM6,MCM7,MDM2,MYC,GADD45B,ATM,WEE2,ORC1,ORC2,ORC4,ORC5,PCNA,FZR1,ANAPC5,ANAPC7,ANAPC11,PLK1,ATR,PRKDC,RAD21,RB1,RBL1,RBL2,CCND1,ANAPC1,SKP1,SKP2,,,BUB1,BUB1B,TFDP1,TFDP2,TGFB1,TGFB2,TGFB3,TP53,TTK,SKP1P2,,WEE1,YWHAB,YWHAE,YWHAG,YWHAH,YWHAZ,ZBTB17,SMC1A,CDC7,CDC45,MAD1L1,CUL1,CCNB3,CDC14B,CDC14A,CDC23,CDC16,CCNA2,CCNA1,CCNB1,CCND2,CCND3,CCNE1,CCNH,PKMYT1,SMC3,CCNB2,CCNE2,BUB3,PTTG1,ESPL1,CDK1,CDC6,CDC20,CDC25A,CDC25B,CDC25C,CDC27,RBX1"
cell_cycle_genes <- trimws(unlist(strsplit(list_cell_cycle, split=",")))
cell_cycle_genes <- cell_cycle_genes[cell_cycle_genes != ""]
present_cc_genes <- cell_cycle_genes[cell_cycle_genes %in% rownames(tumor_macrophages)]
if(length(present_cc_genes)==0) {
  cat("No cell cycle genes found in the dataset (COL6 grouping).\n")
} else {
  p_cellcycle_col6 <- DotPlot(tumor_macrophages, 
                              features = present_cc_genes, 
                              cols = c("lightgrey", "blue"),
                              dot.scale = 8,
                              scale = FALSE) + 
                      RotatedAxis() +
                      ggtitle("Cell Cycle Genes by COL6 Status") +
                      theme(plot.title = element_text(size = 14, face = "bold"),
                            axis.text.x = element_text(angle = 45, hjust = 1))
  plot_width <- max(12, length(present_cc_genes) * 0.4)
  ggsave("pictures/COL6_cell_cycle_genes_tumor_macrophages.pdf", p_cellcycle_col6, width = plot_width, height = 7)
  
  p_cc_heatmap_col6 <- DoHeatmap(tumor_macrophages, 
                                 features = present_cc_genes, 
                                 group.by = current_status,
                                 disp.min = -2.5, 
                                 disp.max = 2.5) + 
                       ggtitle("Cell Cycle Genes Heatmap (COL6 high vs low)") +
                       theme(plot.title = element_text(size = 14, face = "bold"))
  ggsave("pictures/COL6_cell_cycle_genes_heatmap_tumor_macrophages.pdf", p_cc_heatmap_col6, width = 14, height = 10)
}

# --- Additional Marker Analyses for COL6-based grouping ---
# (These analyses mirror the additional markers in the TRAP1 branch,
#  but here we use the COL6 grouping)
# Invasion markers
invasion_markers <- c(
  "MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP13", "MMP14",
  "TIMP1", "TIMP2", "TIMP3", "TIMP4",
  "PLAU", "PLAUR",
  "CTSB", "CTSL", "CTSD", "CTSK",
  "CD44", "ADAM10", "ADAM17"
)
present_invasion <- invasion_markers[invasion_markers %in% rownames(tumor_macrophages)]
if(length(present_invasion) > 0) {
  p_invasion_col6 <- DotPlot(tumor_macrophages, features = present_invasion, group.by = current_status,
                             cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                      RotatedAxis() +
                      ggtitle("Invasion Markers by COL6 Status") +
                      theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_invasion) * 0.4)
  ggsave("pictures/COL6_invasion_markers_tumor_macrophages.pdf", p_invasion_col6, width = plot_width, height = 7)
}

# Metastasis markers
metastasis_markers <- c(
  "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2",
  "RAC1", "RHOA", "CDC42", "S100A4", "MALAT1", "MTA1", "CXCR4", "CCR7"
)
present_metastasis <- metastasis_markers[metastasis_markers %in% rownames(tumor_macrophages)]
if(length(present_metastasis) > 0) {
  p_metastasis_col6 <- DotPlot(tumor_macrophages, features = present_metastasis, group.by = current_status,
                               cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                       RotatedAxis() +
                       ggtitle("Metastasis Markers by COL6 Status") +
                       theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_metastasis) * 0.4)
  ggsave("pictures/COL6_metastasis_markers_tumor_macrophages.pdf", p_metastasis_col6, width = plot_width, height = 7)
}

# Angiogenesis markers
angiogenesis_markers <- c(
  "VEGFA", "VEGFB", "VEGFC", "FIGF", "PGF", "FGF2", "PDGFA", "PDGFB", "ANGPT1", "ANGPT2",
  "KDR", "FLT1", "FLT4", "TEK", "PDGFRA", "PDGFRB",
  "HIF1A", "EPAS1", "ARNT",
  "PECAM1", "CDH5", "ENG",
  "THBS1", "THBS2", "SERPINF1"
)
present_angiogenesis <- angiogenesis_markers[angiogenesis_markers %in% rownames(tumor_macrophages)]
if(length(present_angiogenesis) > 0) {
  p_angiogenesis_col6 <- DotPlot(tumor_macrophages, features = present_angiogenesis, group.by = current_status,
                                 cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                         RotatedAxis() +
                         ggtitle("Angiogenesis Markers by COL6 Status") +
                         theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_angiogenesis) * 0.4)
  ggsave("pictures/COL6_angiogenesis_markers_tumor_macrophages.pdf", p_angiogenesis_col6, width = plot_width, height = 7)
}

# Other collagen genes (excluding COL6 family)
other_collagen_genes <- c(
  "COL1A1", "COL1A2", "COL2A1", "COL3A1", 
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  "COL5A1", "COL5A2", "COL5A3",
  "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3",
  "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", 
  "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1",
  "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1"
)
present_other_collagen <- other_collagen_genes[other_collagen_genes %in% rownames(tumor_macrophages)]
if(length(present_other_collagen) > 0) {
  p_other_collagen_col6 <- DotPlot(tumor_macrophages, features = present_other_collagen, group.by = current_status,
                                   cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                           RotatedAxis() +
                           ggtitle("Other Collagen Genes by COL6 Status") +
                           theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_other_collagen) * 0.4)
  ggsave("pictures/COL6_other_collagen_genes_tumor_macrophages.pdf", p_other_collagen_col6, width = plot_width, height = 7)
}

# DEG Analysis for COL6-based grouping (run before switching identity)
deg_COL6_high_vs_low <- FindMarkers(tumor_macrophages, ident.1 = "COL6-high", ident.2 = "COL6-low", min.pct = 0.25)
deg_COL6_low_vs_high <- FindMarkers(tumor_macrophages, ident.1 = "COL6-low", ident.2 = "COL6-high", min.pct = 0.25)
write.csv(deg_COL6_high_vs_low, "statistics/COL6_high_vs_low_DEGs_tumor_macrophages.csv", row.names = TRUE)
write.csv(deg_COL6_low_vs_high, "statistics/COL6_low_vs_high_DEGs_tumor_macrophages.csv", row.names = TRUE)

#########################################
# Branch 2: Analysis using TRAP1 high/low status
#########################################

# Reset grouping for TRAP1-based analysis
if("TRAP1" %in% rownames(tumor_macrophages)) {
  TRAP1_exp <- GetAssayData(tumor_macrophages, assay = "RNA", layer = "data")["TRAP1", ]
  TRAP1_cutoff <- median(TRAP1_exp)
  tumor_macrophages$TRAP1_status <- ifelse(TRAP1_exp > TRAP1_cutoff, "TRAP1-high", "TRAP1-low")
  cat("TRAP1 cutoff:", round(TRAP1_cutoff, 4), "\n")
} else {
  stop("TRAP1 not found in the dataset!")
}

Idents(tumor_macrophages) <- "TRAP1_status"
current_status <- "TRAP1_status"

# UMAP plot for TRAP1 grouping
p_status_TRAP1 <- DimPlot(tumor_macrophages, group.by = "TRAP1_status", 
                          cols = c("TRAP1-low" = "blue", "TRAP1-high" = "red")) +
                 ggtitle("TRAP1 Status Groups in Tumor-associated Macrophages") +
                 theme(plot.title = element_text(size = 14, face = "bold"))
ggsave("pictures/TRAP1_status_groups_tumor_macrophages.pdf", p_status_TRAP1, width = 8, height = 7)

# Cell Cycle Gene Analysis for TRAP1 grouping
present_cc_genes_TRAP1 <- cell_cycle_genes[cell_cycle_genes %in% rownames(tumor_macrophages)]
if(length(present_cc_genes_TRAP1)==0) {
  cat("No cell cycle genes found in the dataset (TRAP1 grouping).\n")
} else {
  p_cellcycle_TRAP1 <- DotPlot(tumor_macrophages, 
                               features = present_cc_genes_TRAP1, 
                               group.by = "TRAP1_status",
                               cols = c("lightgrey", "blue"),
                               dot.scale = 8,
                               scale = FALSE) + 
                      RotatedAxis() +
                      ggtitle("Cell Cycle Genes by TRAP1 Status") +
                      theme(plot.title = element_text(size = 14, face = "bold"),
                            axis.text.x = element_text(angle = 45, hjust = 1))
  plot_width <- max(12, length(present_cc_genes_TRAP1) * 0.4)
  ggsave("pictures/TRAP1_cell_cycle_genes_dotplot_tumor_macrophages.pdf", p_cellcycle_TRAP1, width = plot_width, height = 7)
  
  p_cc_heatmap_TRAP1 <- DoHeatmap(tumor_macrophages, 
                                  features = present_cc_genes_TRAP1, 
                                  group.by = "TRAP1_status",
                                  disp.min = -2.5, 
                                  disp.max = 2.5) + 
                        ggtitle("Cell Cycle Genes Heatmap (TRAP1 high vs low)") +
                        theme(plot.title = element_text(size = 14, face = "bold"))
  ggsave("pictures/TRAP1_cell_cycle_genes_heatmap_tumor_macrophages.pdf", p_cc_heatmap_TRAP1, width = 14, height = 10)
}

# --- Additional Marker Analyses for TRAP1-based grouping ---
# Invasion markers
invasion_markers <- c(
  "MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP13", "MMP14",
  "TIMP1", "TIMP2", "TIMP3", "TIMP4",
  "PLAU", "PLAUR",
  "CTSB", "CTSL", "CTSD", "CTSK",
  "CD44", "ADAM10", "ADAM17"
)
present_invasion <- invasion_markers[invasion_markers %in% rownames(tumor_macrophages)]
if(length(present_invasion) > 0) {
  p_invasion_TRAP1 <- DotPlot(tumor_macrophages, features = present_invasion, group.by = current_status,
                              cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                       RotatedAxis() +
                       ggtitle("Invasion Markers by TRAP1 Status") +
                       theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_invasion) * 0.4)
  ggsave("pictures/TRAP1_invasion_markers_tumor_macrophages.pdf", p_invasion_TRAP1, width = plot_width, height = 7)
}

# Metastasis markers
metastasis_markers <- c(
  "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2",
  "RAC1", "RHOA", "CDC42", "S100A4", "MALAT1", "MTA1", "CXCR4", "CCR7"
)
present_metastasis <- metastasis_markers[metastasis_markers %in% rownames(tumor_macrophages)]
if(length(present_metastasis) > 0) {
  p_metastasis_TRAP1 <- DotPlot(tumor_macrophages, features = present_metastasis, group.by = current_status,
                                cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                         RotatedAxis() +
                         ggtitle("Metastasis Markers by TRAP1 Status") +
                         theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_metastasis)*0.4)
  ggsave("pictures/TRAP1_metastasis_markers_tumor_macrophages.pdf", p_metastasis_TRAP1, width = plot_width, height = 7)
}

# Angiogenesis markers
angiogenesis_markers <- c(
  "VEGFA", "VEGFB", "VEGFC", "FIGF", "PGF", "FGF2", "PDGFA", "PDGFB", "ANGPT1", "ANGPT2",
  "KDR", "FLT1", "FLT4", "TEK", "PDGFRA", "PDGFRB",
  "HIF1A", "EPAS1", "ARNT",
  "PECAM1", "CDH5", "ENG",
  "THBS1", "THBS2", "SERPINF1"
)
present_angiogenesis <- angiogenesis_markers[angiogenesis_markers %in% rownames(tumor_macrophages)]
if(length(present_angiogenesis) > 0) {
  p_angiogenesis_TRAP1 <- DotPlot(tumor_macrophages, features = present_angiogenesis, group.by = current_status,
                                  cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                           RotatedAxis() +
                           ggtitle("Angiogenesis Markers by TRAP1 Status") +
                           theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_angiogenesis)*0.4)
  ggsave("pictures/TRAP1_angiogenesis_markers_tumor_macrophages.pdf", p_angiogenesis_TRAP1, width = plot_width, height = 7)
}

# Other collagen genes (excluding COL6 family)
other_collagen_genes <- c(
  "COL1A1", "COL1A2", "COL2A1", "COL3A1", 
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  "COL5A1", "COL5A2", "COL5A3",
  "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3",
  "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", 
  "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1",
  "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1"
)
present_other_collagen <- other_collagen_genes[other_collagen_genes %in% rownames(tumor_macrophages)]
if(length(present_other_collagen) > 0) {
  p_other_collagen_TRAP1 <- DotPlot(tumor_macrophages, features = present_other_collagen, group.by = current_status,
                                    cols = c("lightgrey", "blue"), dot.scale = 8, scale = FALSE) +
                             RotatedAxis() +
                             ggtitle("Other Collagen Genes by TRAP1 Status") +
                             theme(plot.title = element_text(size = 14, face="bold"))
  plot_width <- max(12, length(present_other_collagen)*0.4)
  ggsave("pictures/TRAP1_other_collagen_genes_tumor_macrophages.pdf", p_other_collagen_TRAP1, width = plot_width, height = 7)
}
# DEG Analysis for COL6-based grouping
Idents(tumor_macrophages) <- "COL6_status"
deg_COL6_high_vs_low <- FindMarkers(tumor_macrophages, ident.1 = "COL6-high", ident.2 = "COL6-low", min.pct = 0.25)
deg_COL6_low_vs_high <- FindMarkers(tumor_macrophages, ident.1 = "COL6-low", ident.2 = "COL6-high", min.pct = 0.25)
write.csv(deg_COL6_high_vs_low, "statistics/COL6_high_vs_low_DEGs_tumor_macrophages.csv", row.names = TRUE)
write.csv(deg_COL6_low_vs_high, "statistics/COL6_low_vs_high_DEGs_tumor_macrophages.csv", row.names = TRUE)

# Volcano plot for COL6-based DEG (High vs Low)
deg_COL6 <- deg_COL6_high_vs_low
deg_COL6$gene <- rownames(deg_COL6)
sig_deg_COL6 <- subset(deg_COL6, p_val_adj < 0.05)

top100_COL6 <- head(sig_deg_COL6[order(-abs(sig_deg_COL6$avg_log2FC)), ], 100)

p_volcano_COL6 <- ggplot(deg_COL6, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  geom_point(data = sig_deg_COL6, color = "red", alpha = 0.7) +
  geom_text_repel(data = top100_COL6, aes(label = gene), size = 3, max.overlaps = Inf) +
  theme_minimal() +
  ggtitle("Volcano Plot for COL6 DEG (High vs Low)") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value")

ggsave("pictures/COL6_DEG_volcano_top100_tumor_macrophages.pdf", p_volcano_COL6, width = 8, height = 7)

#########################################
# (Remaining COL6 marker dotplots go here ...)
#########################################

# -----------------------------------------------
# Branch 2: Analysis using TRAP1 high/low status (including additional markers)
# -----------------------------------------------

# DEG Analysis for TRAP1-based grouping
Idents(tumor_macrophages) <- "TRAP1_status"
print(unique(Idents(tumor_macrophages)))  # should show "TRAP1-high" and "TRAP1-low"
deg_TRAP1_high_vs_low <- FindMarkers(tumor_macrophages, ident.1 = "TRAP1-high", ident.2 = "TRAP1-low", min.pct = 0.25)
deg_TRAP1_low_vs_high <- FindMarkers(tumor_macrophages, ident.1 = "TRAP1-low", ident.2 = "TRAP1-high", min.pct = 0.25)
write.csv(deg_TRAP1_high_vs_low, "statistics/TRAP1_high_vs_low_DEGs_tumor_macrophages.csv", row.names = TRUE)
write.csv(deg_TRAP1_low_vs_high, "statistics/TRAP1_low_vs_high_DEGs_tumor_macrophages.csv", row.names = TRUE)

# Volcano plot for TRAP1-based DEG (High vs Low)
deg_TRAP1 <- deg_TRAP1_high_vs_low
deg_TRAP1$gene <- rownames(deg_TRAP1)
sig_deg_TRAP1 <- subset(deg_TRAP1, p_val_adj < 0.05)

top100_TRAP1 <- head(sig_deg_TRAP1[order(-abs(sig_deg_TRAP1$avg_log2FC)), ], 100)

p_volcano_TRAP1 <- ggplot(deg_TRAP1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  geom_point(data = sig_deg_TRAP1, color = "red", alpha = 0.7) +
  geom_text_repel(data = top100_TRAP1, aes(label = gene), size = 3, max.overlaps = Inf) +
  theme_minimal() +
  ggtitle("Volcano Plot for TRAP1 DEG (High vs Low)") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value")

ggsave("pictures/TRAP1_DEG_volcano_top100_tumor_macrophages.pdf", p_volcano_TRAP1, width = 8, height = 7)
cat("Analysis complete!\n")