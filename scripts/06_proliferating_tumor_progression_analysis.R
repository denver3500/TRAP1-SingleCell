# MPNST TRAP1 Expression Analysis in Proliferating Cancer Cells
# Analyzes invasion, metastasis, angiogenesis, and collagen genes in TRAP1-high vs TRAP1-low cells

library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
library(ggrepel)

# Create directories for outputs if they don't exist
if (!dir.exists("pictures")) dir.create("pictures")
if (!dir.exists("statistics")) dir.create("statistics")

# Load the pre-processed Seurat object
data.merged <- readRDS("data_mpnst.rds")

#-----------------------------------------------
# 1. Cell subsetting workflow
#-----------------------------------------------

# First subset to only MPNST cells
print("Subsetting MPNST cells...")
mpnst_only <- subset(data.merged, subset = orig.ident == "MPNST")
print(paste("Total MPNST cells:", ncol(mpnst_only)))

# Preserve cell type annotations from the merged object
Idents(mpnst_only) <- Idents(data.merged)

# Next, subset to only Proliferating cancer cells
print("Subsetting to Proliferating cancer cells...")
proliferating_cancer <- subset(mpnst_only, idents = "Proliferating cancer cell")
print(paste("Proliferating cancer cells:", ncol(proliferating_cancer)))

# Check if we have enough cells for analysis
if(ncol(proliferating_cancer) < 10) {
  stop("Not enough Proliferating cancer cells for meaningful analysis")
}

#-----------------------------------------------
# 2. TRAP1 high/low categorization
#-----------------------------------------------

# Verify TRAP1 gene exists
if(!"TRAP1" %in% rownames(proliferating_cancer)) {
  stop("TRAP1 gene not found in dataset. Check gene spelling/capitalization.")
}

# Get TRAP1 expression and determine high/low groups
DefaultAssay(proliferating_cancer) <- "RNA"
# Make sure data is normalized but not scaled to avoid the warning
proliferating_cancer <- NormalizeData(proliferating_cancer)
# Use layer instead of slot (fixes deprecation warning)
trap1_exp <- GetAssayData(proliferating_cancer, layer = "data")["TRAP1", ]
trap1_cutoff <- median(trap1_exp)
print(paste("TRAP1 expression cutoff (median):", round(trap1_cutoff, 4)))

# Add TRAP1 classification as metadata
proliferating_cancer$TRAP1_status <- ifelse(trap1_exp > trap1_cutoff, "TRAP1-high", "TRAP1-low")
table_result <- table(proliferating_cancer$TRAP1_status)
print("TRAP1 expression groups:")
print(table_result)

# Set identity to TRAP1 status for visualization and analysis
Idents(proliferating_cancer) <- "TRAP1_status"

#-----------------------------------------------
# 3. Define gene categories
#-----------------------------------------------

# Invasion markers
invasion_markers <- c(
  # MMPs
  "MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP13", "MMP14",
  # TIMPs
  "TIMP1", "TIMP2", "TIMP3", "TIMP4",
  # uPA system
  "PLAU", "PLAUR",
  # Cathepsins
  "CTSB", "CTSL", "CTSD", "CTSK",
  # Other invasion markers
  "CD44", "ADAM10", "ADAM17"
)

# Metastasis markers
metastasis_markers <- c(
  # EMT markers
  "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2",
  # Cell motility
  "RAC1", "RHOA", "CDC42",
  # Metastasis-specific genes
  "S100A4", "MALAT1", "MTA1",
  # Chemokine receptors
  "CXCR4", "CCR7"
)

# Angiogenesis markers
angiogenesis_markers <- c(
  # Growth factors
  "VEGFA", "VEGFB", "VEGFC", "FIGF", "PGF", "FGF2", "PDGFA", "PDGFB", "ANGPT1", "ANGPT2",
  # Receptors
  "KDR", "FLT1", "FLT4", "TEK", "PDGFRA", "PDGFRB",
  # Hypoxia-related
  "HIF1A", "EPAS1", "ARNT",
  # Endothelial markers
  "PECAM1", "CDH5", "ENG",
  # Inhibitors
  "THBS1", "THBS2", "SERPINF1"
)

# Collagen genes
collagen_genes <- c(
  # COL1 family
  "COL1A1", "COL1A2",
  # COL2 family
  "COL2A1",
  # COL3 family
  "COL3A1",
  # COL4 family
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  # COL5 family
  "COL5A1", "COL5A2", "COL5A3",
  # COL6 family
  "COL6A1", "COL6A2", "COL6A3", "COL6A4", "COL6A5", "COL6A6",
  # COL7-28 families
  "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COL10A1",
  "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", 
  "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1",
  "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1"
)

#-----------------------------------------------
# 4. TRAP1 expression visualization
#-----------------------------------------------

# Create TRAP1 expression UMAP
p1 <- FeaturePlot(proliferating_cancer, features = "TRAP1") + 
  ggtitle("TRAP1 Expression in MPNST Proliferating Cancer Cells") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Create TRAP1 status group UMAP
p2 <- DimPlot(proliferating_cancer, group.by = "TRAP1_status", 
              cols = c("TRAP1-low" = "blue", "TRAP1-high" = "red")) +
  ggtitle("TRAP1 Status Groups in Proliferating Cancer Cells") +
  theme(plot.title = element_text(size = 14, face = "bold"))

# Save visualizations
ggsave("pictures/TRAP1_expression_proliferating.pdf", p1, width = 8, height = 7)
ggsave("pictures/TRAP1_status_groups_proliferating.pdf", p2, width = 8, height = 7)

#-----------------------------------------------
# 5. Gene category analysis function
#-----------------------------------------------

analyze_gene_category <- function(gene_list, category_name) {
  # Check which genes are present in the dataset
  genes_present <- gene_list[gene_list %in% rownames(proliferating_cancer)]
  print(paste("Found", length(genes_present), "out of", length(gene_list), category_name, "genes in dataset"))
  
  if(length(genes_present) == 0) {
    print(paste("No", category_name, "genes found in dataset"))
    return(NULL)
  }
  
  # If category is Collagen, order numerically
  if(category_name == "Collagen") {
    # Attempt to extract numeric parts and order
    tryCatch({
      # Extract the numeric part from each collagen gene
      col_numbers <- as.numeric(str_extract(genes_present, "(?<=COL)\\d+"))
      col_subnumbers <- as.numeric(str_extract(genes_present, "(?<=A)\\d+"))
      
      # Create a dataframe for sorting
      col_df <- data.frame(
        gene = genes_present,
        main_num = col_numbers,
        sub_num = col_subnumbers,
        stringsAsFactors = FALSE
      )
      
      # Sort by main number, then by sub-number
      col_df <- col_df %>% arrange(main_num, sub_num)
      genes_present <- col_df$gene
    }, error = function(e) {
      print("Warning: Could not order collagen genes numerically. Using alphabetical order.")
      genes_present <- sort(genes_present)
    })
  }
  
  # Create dotplot - FIXED with scale=FALSE to prevent scaling warning
  dot_plot <- DotPlot(proliferating_cancer, 
                      features = genes_present, 
                      group.by = "TRAP1_status",
                      cols = c("lightgrey", "blue"),
                      dot.scale = 8,
                      scale = FALSE) + 
             RotatedAxis() +
             ggtitle(paste(category_name, "Genes by TRAP1 Status in Proliferating Cancer")) +
             theme(plot.title = element_text(size = 14, face = "bold"),
                   axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save dotplot
  plot_width <- max(12, length(genes_present) * 0.4)  # Adjust width based on number of genes
  ggsave(paste0("pictures/TRAP1_proliferating_", tolower(category_name), "_genes.pdf"), 
         dot_plot, width = plot_width, height = 7)
  
  # Perform statistical tests
  gene_stats <- data.frame(
    gene = character(),
    p_value = numeric(),
    avg_log2FC = numeric(),
    pct.1 = numeric(),
    pct.2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(gene in genes_present) {
    # Get expression values for high and low groups - using layer instead of slot
    high_exp <- GetAssayData(proliferating_cancer, layer = "data")[gene, proliferating_cancer$TRAP1_status == "TRAP1-high"]
    low_exp <- GetAssayData(proliferating_cancer, layer = "data")[gene, proliferating_cancer$TRAP1_status == "TRAP1-low"]
    
    # Statistical test
    test_result <- wilcox.test(high_exp, low_exp)
    
    # Calculate log2 fold change
    avg_high <- mean(expm1(high_exp))
    avg_low <- mean(expm1(low_exp))
    log2FC <- log2((avg_high + 0.1) / (avg_low + 0.1))  # Add small constant to avoid division by zero
    
    # Calculate percent expressed
    pct.1 <- sum(high_exp > 0) / length(high_exp)
    pct.2 <- sum(low_exp > 0) / length(low_exp)
    
    # Add to stats dataframe
    gene_stats <- rbind(gene_stats, 
                        data.frame(gene = gene, 
                                   p_value = test_result$p.value, 
                                   avg_log2FC = log2FC,
                                   pct.1 = pct.1,
                                   pct.2 = pct.2))
  }
  
  # Adjust p-values
  gene_stats$adj_p_value <- p.adjust(gene_stats$p_value, method = "BH")
  
  # Sort by statistical significance
  gene_stats <- gene_stats %>% arrange(p_value)
  
  # Save statistics
  write.csv(gene_stats, 
            paste0("statistics/TRAP1_proliferating_", tolower(category_name), "_genes_statistics.csv"), 
            row.names = FALSE)
  
  # Print significant genes
  sig_genes <- gene_stats %>% filter(adj_p_value < 0.05)
  if(nrow(sig_genes) > 0) {
    print(paste("Significant", category_name, "genes (adj p < 0.05):"))
    print(sig_genes)
  } else {
    print(paste("No significant", category_name, "genes found"))
  }
  
  return(list(plot = dot_plot, stats = gene_stats, genes = genes_present))
}

#-----------------------------------------------
# 6. Analyze each gene category
#-----------------------------------------------

print("Analyzing invasion markers...")
invasion_results <- analyze_gene_category(invasion_markers, "Invasion")

print("Analyzing metastasis markers...")
metastasis_results <- analyze_gene_category(metastasis_markers, "Metastasis")

print("Analyzing angiogenesis markers...")
angiogenesis_results <- analyze_gene_category(angiogenesis_markers, "Angiogenesis")

print("Analyzing collagen genes...")
collagen_results <- analyze_gene_category(collagen_genes, "Collagen")

#-----------------------------------------------
# 7. Create combined visualization of top genes
#-----------------------------------------------

# Function to extract top significant genes from each category (up to 5 from each)
get_top_genes <- function(results, max_genes = 5) {
  if(is.null(results) || !("stats" %in% names(results))) return(character(0))
  return(results$stats %>% 
         filter(adj_p_value < 0.05) %>% 
         arrange(p_value) %>% 
         head(max_genes) %>% 
         pull(gene))
}

top_invasion <- get_top_genes(invasion_results)
top_metastasis <- get_top_genes(metastasis_results)
top_angiogenesis <- get_top_genes(angiogenesis_results)
top_collagen <- get_top_genes(collagen_results)

# Combine all significant genes
all_top_genes <- unique(c(top_invasion, top_metastasis, top_angiogenesis, top_collagen))

if(length(all_top_genes) > 0) {
  # Create combined dotplot of top significant genes - FIXED with scale=FALSE
  p_top <- DotPlot(proliferating_cancer, 
                   features = all_top_genes, 
                   group.by = "TRAP1_status",
                   cols = c("lightgrey", "blue"),
                   dot.scale = 8,
                   scale = FALSE) + 
          RotatedAxis() +
          ggtitle("Top Significant Genes by TRAP1 Status in Proliferating Cancer") +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save combined plot
  plot_width <- max(12, length(all_top_genes) * 0.4)
  ggsave("pictures/TRAP1_top_significant_genes_proliferating.pdf", 
         p_top, width = plot_width, height = 7)
  
  # Create heatmap of top genes - avoid scaling warning by using normalized data
  p_heatmap <- DoHeatmap(proliferating_cancer, 
                         features = all_top_genes, 
                         group.by = "TRAP1_status",
                         disp.min = -2.5, 
                         disp.max = 2.5) + 
              ggtitle("Top Significant Genes: TRAP1-high vs TRAP1-low in Proliferating Cancer")
  
  ggsave("pictures/TRAP1_top_genes_heatmap_proliferating.pdf", 
         p_heatmap, width = 14, height = 10)
} else {
  print("No significant genes found across categories")
}

#-----------------------------------------------
# 8. Find all differentially expressed genes
#-----------------------------------------------
# Ensure the active identities are set for TRAP1 high/low in your proliferating_cancer object
Idents(proliferating_cancer) <- "TRAP1_status"
print(unique(Idents(proliferating_cancer)))  # should show "TRAP1-high" and "TRAP1-low"

# DEG Analysis for TRAP1-based grouping (High vs Low)
deg_TRAP1_high_vs_low <- FindMarkers(proliferating_cancer, 
                                     ident.1 = "TRAP1-high", 
                                     ident.2 = "TRAP1-low", 
                                     min.pct = 0.25)
deg_TRAP1_high_vs_low$gene <- rownames(deg_TRAP1_high_vs_low)
write.csv(deg_TRAP1_high_vs_low, "statistics/TRAP1_high_vs_low_DEGs_proliferating_cancer.csv", row.names = FALSE)

# DEG for TRAP1-low vs TRAP1-high
deg_TRAP1_low_vs_high <- FindMarkers(proliferating_cancer, 
                                     ident.1 = "TRAP1-low", 
                                     ident.2 = "TRAP1-high", 
                                     min.pct = 0.25)
deg_TRAP1_low_vs_high$gene <- rownames(deg_TRAP1_low_vs_high)
write.csv(deg_TRAP1_low_vs_high, "statistics/TRAP1_low_vs_high_DEGs_proliferating_cancer.csv", row.names = FALSE)

# --- Volcano plot for top 100 DEG genes (TRAP1-high vs TRAP1-low) ---

# Remove TRAP1 gene from DEG results before plotting
deg_TRAP1_no_TRAP1 <- deg_TRAP1_high_vs_low[rownames(deg_TRAP1_high_vs_low) != "TRAP1", ]
deg_TRAP1_no_TRAP1$gene <- rownames(deg_TRAP1_no_TRAP1)

# Filter for significant DEGs (adjusted p-value < 0.05) without TRAP1
sig_deg_TRAP1_no_TRAP1 <- subset(deg_TRAP1_no_TRAP1, p_val_adj < 0.05)

# Order by absolute log2FC and take the top 100 
top100_TRAP1_no_TRAP1 <- head(sig_deg_TRAP1_no_TRAP1[order(-abs(sig_deg_TRAP1_no_TRAP1$avg_log2FC)), ], 100)

# Create volcano plot using the DEG results without TRAP1, but label only the top 100 genes
p_volcano_TRAP1 <- ggplot(deg_TRAP1_no_TRAP1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  geom_point(data = sig_deg_TRAP1_no_TRAP1, color = "red", alpha = 0.7) +
  geom_text_repel(data = top100_TRAP1_no_TRAP1, aes(label = gene), size = 3, max.overlaps = Inf) +
  theme_minimal() +
  ggtitle("TRAP1 DEG Volcano Plot (High vs Low) without TRAP1 in Proliferating Cancer") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value")

ggsave("pictures/TRAP1_DEG_volcano_top100_noTRAP1_proliferating_cancer.pdf", 
       p_volcano_TRAP1, width = 8, height = 7)

cat("TRAP1 DEG volcano plot (excluding TRAP1 gene) complete for Proliferating cancer cells!\n")