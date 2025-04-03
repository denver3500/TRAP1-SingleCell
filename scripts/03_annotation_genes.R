# MPNST Cluster Annotation & Top Gene Expression Visualization
# For each cluster, create a combined figure with:
#   Row 1: Integrated UMAP map showing cluster annotation (without extra group labels)
#   Row 2: Three FeaturePlots for the top 3 genes (expression plotted on the same UMAP)

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)

# Create output directory if needed
annotation_dir <- "pictures/annotation"
if (!dir.exists(annotation_dir)) dir.create(annotation_dir, recursive = TRUE)

# Load the pre-processed Seurat object
data.merged <- readRDS("data_mpnst.rds")

cluster_name_map <- c(
    "0" = "Cancer-associated fibroblast",
    "1" = "Cancer-associated fibroblast",
    "2" = "Schwann cell-derived cancer",
    "3" = "Cancer-associated fibroblast",
    "4" = "Schwann cell-derived cancer",
    "5" = "Tumor-associated macrophage",
    "6" = "Tumor endothelium",
    "7" = "T lymphocyte",
    "8" = "Dendritic cell",
    "9" = "Cancer-associated fibroblast",
    "10" = "Tumor endothelium",
    "11" = "Plasma cell",
    "12" = "Non-coding RNA cluster",
    "13" = "Tumor-associated macrophage",
    "14" = "Pericyte",
    "15" = "Cytotoxic T cell",
    "16" = "Inflammatory myeloid cell",
    "17" = "Mesenchymal-like cancer cell",
    "18" = "Tumor-associated macrophage",
    "19" = "Antigen-presenting cell",
    "20" = "Tumor endothelium",
    "21" = "Monocyte",
    "22" = "Tumor-associated macrophage",
    "23" = "Proliferating cancer cell",
    "24" = "Plasma cell",
    "25" = "NK cell",
    "26" = "Muscle cell",
    "27" = "Lymphatic endothelium",
    "28" = "Melanocyte-like cancer cell",
    "29" = "Plasmacytoid dendritic cell"
)

# Read the single cluster genes file
cluster_genes_file <- "statistics/cluster_genes.txt"
if (!file.exists(cluster_genes_file)) {
  stop("Cluster genes file not found: ", cluster_genes_file)
}
all_clusters_data <- readLines(cluster_genes_file)

# Parse file lines to create a map of cluster id -> gene list
cluster_gene_map <- list()
for (line in all_clusters_data) {
  parts <- str_split(line, ":")[[1]]
  if (length(parts) == 2) {
    cluster_id <- as.numeric(str_extract(parts[1], "\\d+"))
    genes <- str_split(parts[2], ",")[[1]]
    genes <- trimws(genes)
    cluster_gene_map[[as.character(cluster_id)]] <- genes
  }
}

if (length(cluster_gene_map) == 0) {
  stop("No cluster data found in file")
}
print(paste("Found data for", length(cluster_gene_map), "clusters"))

# For each cluster (from 0 to 29)
for (cluster_id in 0:29) {
  cluster_id_str <- as.character(cluster_id)
  cat("Processing cluster", cluster_id, "\n")
  
  if (!cluster_id_str %in% names(cluster_gene_map)) {
    warning(paste("No data found for cluster", cluster_id, "- skipping"))
    next
  }
  
  # Get top 3 genes for the current cluster
  cluster_genes <- cluster_gene_map[[cluster_id_str]]
  top_genes <- cluster_genes[1:min(3, length(cluster_genes))]
  
  # Verify that these genes exist in the Seurat object
  missing_genes <- top_genes[!top_genes %in% rownames(data.merged)]
  if (length(missing_genes) > 0) {
    warning(paste("The following genes were not found in the dataset for cluster", cluster_id, ":\n", 
                  paste(missing_genes, collapse = ", ")))
    top_genes <- top_genes[top_genes %in% rownames(data.merged)]
    if (length(top_genes) == 0) {
      warning(paste("No valid genes for cluster", cluster_id, "- skipping"))
      next
    }
  }
  
  # Create the integrated UMAP plot (top panel)
  # To avoid extra labels (like "Unselected" and "Group_1"), we plot the basic UMAP.
  p_cluster <- DimPlot(data.merged, reduction = "umap", label = TRUE) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    ggtitle(paste("Integrated Map - Cluster", cluster_id)) +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # Optional: highlight current cluster if the numeric identity is stored in "seurat_clusters"
  if ("seurat_clusters" %in% colnames(data.merged@meta.data)) {
    cells_in_cluster <- rownames(data.merged@meta.data)[
      data.merged@meta.data$seurat_clusters == cluster_id
    ]
    if (length(cells_in_cluster) > 0) {
      p_cluster <- p_cluster + 
        DimPlot(data.merged, cells.highlight = cells_in_cluster, reduction = "umap",
                cols.highlight = "red", sizes.highlight = 1, group.by = "seurat_clusters") +
        NoLegend()
    }
  }
  
  # Create FeaturePlots for each of the top 3 genes (using the same UMAP reduction)
  gene_plots <- list()
  for (gene in top_genes) {
    temp_field <- paste0("temp_expr_", gene)
    data.merged[[temp_field]] <- GetAssayData(data.merged, assay = "RNA", layer = "data")[gene, ]
    
    p_gene <- FeaturePlot(data.merged, features = temp_field, reduction = "umap") +
      labs(x = "UMAP 1", y = "UMAP 2", color = gene) +
      ggtitle(gene) +
      theme(plot.title = element_text(size = 14, face = "bold"))
    
    data.merged[[temp_field]] <- NULL
    gene_plots[[gene]] <- p_gene
  }
  
  # Combine the three gene plots side-by-side (bottom row)
  if (length(gene_plots) == 1) {
    bottom_row <- gene_plots[[1]]
  } else if (length(gene_plots) == 2) {
    bottom_row <- gene_plots[[1]] | gene_plots[[2]]
  } else if (length(gene_plots) == 3) {
    bottom_row <- gene_plots[[1]] | gene_plots[[2]] | gene_plots[[3]]
  }
  
  # Combine top and bottom rows (2 rows total) using patchwork
  combined_plot <- p_cluster / bottom_row
  
  # Include current cluster name if available
  current_name <- ifelse(cluster_id_str %in% names(cluster_name_map),
                         cluster_name_map[[cluster_id_str]],
                         paste("Cluster", cluster_id_str))
  
  combined_plot <- combined_plot + 
    plot_annotation(title = paste(current_name, ":", paste(top_genes, collapse = ", ")),
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  
  # Save the combined plot as PDF and PNG
  ggsave(filename = paste0(annotation_dir, "/cluster_", cluster_id, "_annotation.pdf"),
         plot = combined_plot, width = 18, height = 10)
  ggsave(filename = paste0(annotation_dir, "/cluster_", cluster_id, "_annotation.png"),
         plot = combined_plot, width = 18, height = 10, dpi = 150)
  
  cat("  - Created visualization for", current_name, "with genes:", paste(top_genes, collapse = ", "), "\n")
}

cat("Cluster annotation visualizations complete!\n")