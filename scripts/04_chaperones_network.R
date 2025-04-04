library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)
library(pheatmap)
library(gridExtra)

# Load the pre-processed Seurat object
data.merged <- readRDS("data_mpnst.rds")

# Subset based on orig.ident (assumes NF and MPNST are stored in orig.ident)
nf <- subset(data.merged, subset = orig.ident == "NF")
mpnst <- subset(data.merged, subset = orig.ident == "MPNST")

# Define the list of chaperone genes (your genes of interest)
chaperone_genes <- c("HSPD1", "HSPE1", "CCT6A", "HSPA9", "DNLZ", 
                     "GRPEL1", "GRPEL2", "DNAJA3", "HSP90AA1", "HSP90AB1", 
                     "TRAP1", "DNAJC19", "HSCB", "SDHAF1", "CLPX", 
                     "HSPB1", "HSPB8", "HSPB7", "CRYAB")

# Keep only the genes that are present in both subsets (for comparison)
common_genes <- chaperone_genes[
  chaperone_genes %in% rownames(nf) & chaperone_genes %in% rownames(mpnst)
]
if (length(common_genes) == 0) {
  stop("None of the specified chaperone genes were found in both NF and MPNST subsets.")
}

# Retrieve normalized expression matrices using the layer argument
nf_expr <- as.matrix(GetAssayData(nf, assay = "RNA", layer = "data")[common_genes, ])
mpnst_expr <- as.matrix(GetAssayData(mpnst, assay = "RNA", layer = "data")[common_genes, ])

# Compute the correlation matrices (using Pearson correlation)
nf_cor <- cor(t(nf_expr), method = "pearson")
mpnst_cor <- cor(t(mpnst_expr), method = "pearson")

# Create heatmap objects with silent = TRUE so we can capture the grobs
nf_heat <- pheatmap(nf_cor, 
                    cluster_rows = TRUE, 
                    cluster_cols = TRUE, 
                    display_numbers = TRUE, 
                    main = "Chaperone Correlation in NF",
                    silent = TRUE)
mpnst_heat <- pheatmap(mpnst_cor, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE, 
                       display_numbers = TRUE, 
                       main = "Chaperone Correlation in MPNST",
                       silent = TRUE)

# Combine the two heatmaps side by side and save to PDF
pdf("pictures/chaperone_correlation_combined.pdf", width = 14, height = 6)
grid.arrange(nf_heat$gtable, mpnst_heat$gtable, ncol = 2)
dev.off()