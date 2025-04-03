library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)

# Get all folders in raw_data that end with "filtered_feature_bc_matrix"
filtered_dirs <- list.dirs(path = "raw_data", full.names = TRUE, recursive = FALSE)
filtered_dirs <- filtered_dirs[grepl("filtered_feature_bc_matrix$", filtered_dirs)]

# Create a list to store all Seurat objects
seurat_objects <- list()

# Loop through each filtered directory
for (i in seq_along(filtered_dirs)) {
    # Get the directory name
    dir_name <- basename(filtered_dirs[i])

    # Extract sample ID (assuming it's the prefix before "-filtered_feature_bc_matrix")
    sample_id <- gsub("-filtered_feature_bc_matrix$", "", dir_name)

    # Read the 10X data
    data <- Read10X(data.dir = filtered_dirs[i])

    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
        counts = data,
        project = if (grepl("MPNST", sample_id)) "MPNST" else if (grepl("NF", sample_id)) "NF" else sample_id,
        min.cells = 3,
        min.features = 200
    )

    # Add sample ID as metadata
    seurat_obj$sample <- sample_id

    # Store in list
    seurat_objects[[sample_id]] <- seurat_obj
}

# Merge all Seurat objects
if (length(seurat_objects) > 1) {
    data.merged <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])
} else {
    data.merged <- seurat_objects[[1]]
}

# Calculate mitochondrial percentage
data.merged[["percent.mt"]] <- PercentageFeatureSet(data.merged, pattern = "^MT-")

# Generate and save violin plot
vln_plot <- VlnPlot(data.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("pictures/qc_violin_plots.pdf", vln_plot, width = 15, height = 5)

# Generate scatter plots
plot1 <- FeatureScatter(data.merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Combine and save scatter plots
combined_plots <- CombinePlots(plots = list(plot1, plot2))
ggsave("pictures/qc_scatter_plots.pdf", combined_plots, width = 12, height = 6)

data.merged <- subset(data.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 9500 & percent.mt < 20)
data.merged <- FindVariableFeatures(data.merged, nfeatures = 5000)

# run standard anlaysis workflow un-integrated
data.merged <- NormalizeData(data.merged)
data.merged <- FindVariableFeatures(data.merged)
data.merged <- ScaleData(data.merged)
data.merged <- RunPCA(data.merged)
data.merged <- FindNeighbors(data.merged, dims = 1:30, reduction = "pca")
data.merged <- FindClusters(data.merged, resolution = 2, cluster.name = "unintegrated_clusters")
data.merged <- RunUMAP(data.merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Generate and display the first UMAP plot (by orig.ident and clusters)
p1 <- DimPlot(data.merged, reduction = "umap.unintegrated", 
              group.by = c("orig.ident", "seurat_clusters"))
ggsave("pictures/umap_by_sample_and_clusters.pdf", p1, width = 12, height = 8)

# Generate and display the second UMAP plot (split by orig.ident)
p2 <- DimPlot(data.merged, reduction = "umap.unintegrated", 
              group.by = "seurat_clusters", split.by = "orig.ident")
ggsave("pictures/umap_split_by_sample.pdf", p2, width = 12, height = 8)

# Also save using the original filename in the pictures folder
pdf("pictures/umap_plot.pdf", width = 12, height = 8)
DimPlot(data.merged, reduction = "umap.unintegrated", 
        group.by = "seurat_clusters", split.by = "orig.ident")
dev.off()

print("UMAP plots saved to pictures folder")

# integration

data.merged <- IntegrateLayers(
    object = data.merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE
)

# re-join layers after integration
data.merged[["RNA"]] <- JoinLayers(data.merged[["RNA"]])
data.merged <- FindNeighbors(data.merged, reduction = "integrated.cca", dims = 1:30)
data.merged <- FindClusters(data.merged, resolution = 1)
data.merged <- RunUMAP(data.merged, dims = 1:30, reduction = "integrated.cca")

p3 <- DimPlot(data.merged, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
ggsave("pictures/integrated_umap_by_sample_and_clusters.pdf", p3, width = 12, height = 8)

p4 <- DimPlot(data.merged, reduction = "umap", split.by = "orig.ident")
ggsave("pictures/integrated_umap_split_by_sample.pdf", p4, width = 12, height = 8)

print("Integrated UMAP plots saved to pictures folder")

# Find clusters markers. Display top-10 markers for each cluster
all.markers <- FindAllMarkers(data.merged, only.pos = TRUE)

all.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(data.merged, features = top10$gene) + NoLegend()

filtered_markers <- top10 %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

cluster_genes <- filtered_markers %>%
    group_by(cluster) %>%
    summarise(genes = paste(gene, collapse = ",")) %>% # Removed space after comma
    mutate(formatted = paste0("cluster_", cluster, ":", genes)) # Removed space after colon

# Write the formatted output to a text file
writeLines(cluster_genes$formatted, "statistics/cluster_genes.txt")
print("Cluster genes saved to statistics folder")

data.merged <- RenameIdents(data.merged,
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

# Original integrated UMAP with all cells
integrated_umap <- DimPlot(
    object = data.merged,
    reduction = "umap",
    label = TRUE,
    label.size = 6,
    repel = TRUE
) +
    labs(title = "Integrated Single-Cell RNA-seq Clusters")

ggsave("pictures/integrated_annotated_clusters.pdf", integrated_umap, width = 12, height = 10)

# Subset to only MPNST cells
mpnst_only <- subset(data.merged, subset = orig.ident == "MPNST")

# Preserve the cell type annotations from data.merged
mpnst_only$cell_type <- Idents(data.merged)[colnames(mpnst_only)]
Idents(mpnst_only) <- mpnst_only$cell_type

# Recompute dimensionality reduction specifically for MPNST cells
mpnst_only <- NormalizeData(mpnst_only)
mpnst_only <- FindVariableFeatures(mpnst_only)
mpnst_only <- ScaleData(mpnst_only)
mpnst_only <- RunPCA(mpnst_only)
mpnst_only <- FindNeighbors(mpnst_only, dims = 1:30)
mpnst_only <- FindClusters(mpnst_only, resolution = 0.8)
mpnst_only <- RunUMAP(mpnst_only, dims = 1:30)

# Generate two UMAP plots - one with the new MPNST-specific clusters
mpnst_new_clusters <- DimPlot(
    object = mpnst_only,
    group.by = "seurat_clusters",
    reduction = "umap",
    label = TRUE,
    label.size = 6,
    repel = TRUE
) +
    labs(title = "MPNST-specific Clusters")

# Another with the original cell type annotations
mpnst_cell_types <- DimPlot(
    object = mpnst_only,
    group.by = "cell_type",
    reduction = "umap",
    label = TRUE,
    label.size = 6,
    repel = TRUE
) +
    labs(title = "MPNST Cells - Original Annotations")

# Save both plots
ggsave("pictures/MPNST_reclustered.pdf", mpnst_new_clusters, width = 12, height = 10)
ggsave("pictures/MPNST_original_annotations.pdf", mpnst_cell_types, width = 12, height = 10)

# Create a side-by-side comparison
combined_plot <- mpnst_new_clusters + mpnst_cell_types
ggsave("pictures/MPNST_cluster_comparison.pdf", combined_plot, width = 20, height = 10)

print("MPNST clustering analysis complete!")

saveRDS(data.merged, "data_mpnst.rds")

# Count cells by cluster (cell type) and condition
print("Generating cell count statistics...")
cell_counts <- table(Idents(data.merged), data.merged$orig.ident)
print("Cell counts by cluster and condition:")
print(cell_counts)

# Convert to a data frame and add percentages
cell_count_df <- as.data.frame.matrix(cell_counts)
cell_count_df$Total <- rowSums(cell_count_df)
cell_count_df$Cluster <- rownames(cell_count_df)

# Calculate percentages
if("NF" %in% colnames(cell_count_df)) {
  cell_count_df$NF_percent <- round((cell_count_df$NF / sum(cell_count_df$NF)) * 100, 2)
}
if("MPNST" %in% colnames(cell_count_df)) {
  cell_count_df$MPNST_percent <- round((cell_count_df$MPNST / sum(cell_count_df$MPNST)) * 100, 2)
}
cell_count_df$Total_percent <- round((cell_count_df$Total / sum(cell_count_df$Total)) * 100, 2)

# Reorder columns for better readability
col_order <- c("Cluster")
if("NF" %in% colnames(cell_count_df)) {
  col_order <- c(col_order, "NF", "NF_percent")
}
if("MPNST" %in% colnames(cell_count_df)) {
  col_order <- c(col_order, "MPNST", "MPNST_percent")
}
col_order <- c(col_order, "Total", "Total_percent")
cell_count_df <- cell_count_df[, col_order]

# Add a summary row with totals
total_row <- data.frame(Cluster = "Total")
if("NF" %in% colnames(cell_count_df)) {
  total_row$NF <- sum(cell_count_df$NF)
  total_row$NF_percent <- 100
}
if("MPNST" %in% colnames(cell_count_df)) {
  total_row$MPNST <- sum(cell_count_df$MPNST)
  total_row$MPNST_percent <- 100
}
total_row$Total <- sum(cell_count_df$Total)
total_row$Total_percent <- 100

# Combine with the main dataframe
cell_count_df <- rbind(cell_count_df, total_row)

# Write to CSV
write.csv(cell_count_df, "statistics/cell_counts_by_cluster_and_condition.csv", row.names = FALSE)
print("Cell count statistics saved to statistics/cell_counts_by_cluster_and_condition.csv")

# Create a more detailed summary for research paper tables
summary_stats <- data.frame(
  Metric = c(
    "Total Cells", 
    "Number of Cell Types", 
    "Largest Cell Type",
    "Cells in Largest Type",
    "% in Largest Type"
  )
)

# Calculate stats for each condition and overall
if("NF" %in% colnames(cell_count_df)) {
  nf_cells <- sum(cell_counts[, "NF"])
  nf_types <- sum(cell_counts[, "NF"] > 0)
  nf_largest <- rownames(cell_counts)[which.max(cell_counts[, "NF"])]
  nf_largest_count <- max(cell_counts[, "NF"])
  nf_largest_pct <- round((nf_largest_count / nf_cells) * 100, 2)
  
  summary_stats$NF <- c(
    nf_cells,
    nf_types,
    nf_largest,
    nf_largest_count,
    nf_largest_pct
  )
}

if("MPNST" %in% colnames(cell_count_df)) {
  mpnst_cells <- sum(cell_counts[, "MPNST"])
  mpnst_types <- sum(cell_counts[, "MPNST"] > 0)
  mpnst_largest <- rownames(cell_counts)[which.max(cell_counts[, "MPNST"])]
  mpnst_largest_count <- max(cell_counts[, "MPNST"])
  mpnst_largest_pct <- round((mpnst_largest_count / mpnst_cells) * 100, 2)
  
  summary_stats$MPNST <- c(
    mpnst_cells,
    mpnst_types,
    mpnst_largest,
    mpnst_largest_count,
    mpnst_largest_pct
  )
}

# Write summary stats to CSV
write.csv(summary_stats, "statistics/cell_type_summary_statistics.csv", row.names = FALSE)
print("Summary statistics saved to statistics/cell_type_summary_statistics.csv")