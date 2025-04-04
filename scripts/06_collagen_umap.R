library(Seurat)
library(ggplot2)
library(patchwork)


# Load the pre-processed Seurat object
data.merged <- readRDS("data_mpnst.rds")


# Create the cluster UMAP plot (left panel)
p_cluster <- DimPlot(data.merged, reduction = "umap",
                     group.by = "ident",
                     label = TRUE,
                     repel = TRUE) +
             ggtitle("Cell Clusters") +
             theme(plot.title = element_text(size = 14, face = "bold"),
                   legend.position = "right")

# Define the folder where the plots will be saved
picture_folder <- "pictures"

# List of genes to visualize
gene_list <- c("COL6A1", "COL6A2", "COL6A3")

# Loop over the genes to create and save side-by-side UMAP plots
for (gene in gene_list) {
  # Create feature plot (gene expression) UMAP for the current gene
  p_feature <- FeaturePlot(data.merged,
                           features = gene,
                           reduction = "umap",
                           min.cutoff = "q5",
                           max.cutoff = "q95",
                           order = TRUE) +
              ggtitle(paste(gene, "Expression")) +
              theme(plot.title = element_text(size = 14, face = "bold"))
  
  # Combine the cluster and feature plots side-by-side
  combined_plot <- p_cluster + p_feature + plot_layout(ncol = 2)
  
  # Save the combined plot with a gene-specific filename in the pictures folder
  file_name <- paste0(picture_folder, "/", gene, "_expression_combined_UMAP.pdf")
  ggsave(file_name, combined_plot, width = 16, height = 8)
  
  message("Saved UMAP visualization for ", gene, " as ", file_name)
}