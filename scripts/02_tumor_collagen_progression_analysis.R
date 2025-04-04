# MPNST Collagen VI Expression Analysis in Schwann Cell-derived Cancer
# Analyzes invasion, metastasis, angiogenesis, and collagen genes in COL6A1/2/3-high vs low cells

library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
library(reshape2)
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

# Next, subset to only Schwann cell-derived cancer cells
print("Subsetting to Schwann cell-derived cancer cells...")
schwann_cancer <- subset(mpnst_only, idents = "Schwann cell-derived cancer")
print(paste("Schwann cell-derived cancer cells:", ncol(schwann_cancer)))

# Check if we have enough cells for analysis
if(ncol(schwann_cancer) < 10) {
  stop("Not enough Schwann cell-derived cancer cells for meaningful analysis")
}

#-----------------------------------------------
# 2. Collagen VI high/low categorization
#-----------------------------------------------

# Normalize data
DefaultAssay(schwann_cancer) <- "RNA"
schwann_cancer <- NormalizeData(schwann_cancer)

# Check for the presence of each collagen gene
col6_genes <- c("COL6A1", "COL6A2", "COL6A3")
missing_genes <- col6_genes[!col6_genes %in% rownames(schwann_cancer)]

if(length(missing_genes) > 0) {
  warning(paste("The following collagen genes are missing:", paste(missing_genes, collapse=", ")))
  col6_genes <- col6_genes[col6_genes %in% rownames(schwann_cancer)]
  if(length(col6_genes) == 0) {
    stop("No COL6A1/2/3 genes found in the dataset!")
  }
}

# Function to categorize cells based on high/low expression of a gene
categorize_by_expression <- function(gene) {
  if(gene %in% rownames(schwann_cancer)) {
    exp_values <- GetAssayData(schwann_cancer, layer = "data")[gene, ]
    cutoff <- median(exp_values)
    status <- ifelse(exp_values > cutoff, paste0(gene, "-high"), paste0(gene, "-low"))
    print(paste(gene, "expression cutoff (median):", round(cutoff, 4)))
    print(table(status))
    return(list(status = status, cutoff = cutoff))
  } else {
    return(NULL)
  }
}

# Categorize cells for each collagen gene
col6_status <- list()
for(gene in col6_genes) {
  col6_status[[gene]] <- categorize_by_expression(gene)
  schwann_cancer[[paste0(gene, "_status")]] <- col6_status[[gene]]$status
}

# Create a combined collagen VI score (average of z-scores)
if(length(col6_genes) > 1) {
  z_scores <- matrix(0, nrow = length(col6_genes), ncol = ncol(schwann_cancer))
  rownames(z_scores) <- col6_genes
  
  for(i in 1:length(col6_genes)) {
    gene <- col6_genes[i]
    exp_values <- GetAssayData(schwann_cancer, layer = "data")[gene, ]
    z_scores[i, ] <- scale(exp_values)
  }
  
  col6_mean_zscore <- colMeans(z_scores)
  col6_status_combined <- ifelse(col6_mean_zscore > 0, "COL6-high", "COL6-low")
  schwann_cancer$COL6_status <- col6_status_combined
  
  print("Combined COL6 status:")
  print(table(schwann_cancer$COL6_status))
}

# Set default identity to combined COL6 status if available, otherwise use COL6A1
if("COL6_status" %in% colnames(schwann_cancer@meta.data)) {
  Idents(schwann_cancer) <- "COL6_status"
  current_status <- "COL6_status"
} else {
  Idents(schwann_cancer) <- paste0(col6_genes[1], "_status")
  current_status <- paste0(col6_genes[1], "_status")
}

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

# Other collagen genes (excluding COL6 family)
other_collagen_genes <- c(
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
  # COL7-28 families
  "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COL10A1",
  "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", 
  "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1",
  "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1"
)

#-----------------------------------------------
# 4. COL6 expression visualization
#-----------------------------------------------

# Individual FeaturePlots for each COL6 gene - FIXED approach to avoid naming conflicts
feature_plots <- list()
for(gene in col6_genes) {
  # Save current default assay
  current_assay <- DefaultAssay(schwann_cancer)
  
  # Force RNA assay for gene expression
  DefaultAssay(schwann_cancer) <- "RNA"
  
  # Create a temporary metadata field for plotting with a unique name
  temp_field_name <- paste0("temp_expr_", gene)
  schwann_cancer[[temp_field_name]] <- GetAssayData(schwann_cancer, assay = "RNA", layer = "data")[gene, ]
  
  # Plot using the temporary field but label it as the gene
  p <- FeaturePlot(schwann_cancer, features = temp_field_name) +
    ggtitle(paste(gene, "Expression in MPNST Schwann-derived Cancer")) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    labs(color = gene)  # Use the real gene name for the legend
  
  # Clean up the temporary field
  schwann_cancer[[temp_field_name]] <- NULL
  
  # Restore original assay
  DefaultAssay(schwann_cancer) <- current_assay
  
  feature_plots[[gene]] <- p
  ggsave(paste0("pictures/", gene, "_expression.pdf"), p, width = 8, height = 7)
}

# Create status group UMAP
if("COL6_status" %in% colnames(schwann_cancer@meta.data)) {
  p_status <- DimPlot(schwann_cancer, group.by = "COL6_status", 
                cols = c("COL6-low" = "blue", "COL6-high" = "red")) +
    ggtitle("Combined COL6 Status Groups") +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  ggsave("pictures/COL6_status_groups.pdf", p_status, width = 8, height = 7)
}

# Create individual status plots for each gene
for(gene in col6_genes) {
  status_field <- paste0(gene, "_status")
  if(status_field %in% colnames(schwann_cancer@meta.data)) {
    p <- DimPlot(schwann_cancer, group.by = status_field, 
              cols = setNames(c("blue", "red"), c(paste0(gene, "-low"), paste0(gene, "-high")))) +
        ggtitle(paste(gene, "Status Groups")) +
        theme(plot.title = element_text(size = 14, face = "bold"))
    
    ggsave(paste0("pictures/", gene, "_status_groups.pdf"), p, width = 8, height = 7)
  }
}

# Create a correlation heatmap between different collagen VI genes
if(length(col6_genes) > 1) {
  cor_matrix <- matrix(0, nrow = length(col6_genes), ncol = length(col6_genes))
  rownames(cor_matrix) <- col6_genes
  colnames(cor_matrix) <- col6_genes
  
  exp_data <- GetAssayData(schwann_cancer, layer = "data")[col6_genes, ]
  
  for(i in 1:length(col6_genes)) {
    for(j in 1:length(col6_genes)) {
      cor_matrix[i, j] <- cor(exp_data[i, ], exp_data[j, ], method = "spearman")
    }
  }
  
  cor_df <- reshape2::melt(cor_matrix)
  colnames(cor_df) <- c("Gene1", "Gene2", "Correlation")
  
  p_cor <- ggplot(cor_df, aes(x = Gene1, y = Gene2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    geom_text(aes(label = round(Correlation, 2)), color = "black") +
    theme_minimal() +
    ggtitle("Correlation Between Collagen VI Genes") +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  ggsave("pictures/COL6_gene_correlation.pdf", p_cor, width = 8, height = 6)
}

#-----------------------------------------------
# 5. Gene category analysis function
#-----------------------------------------------

analyze_gene_category <- function(gene_list, category_name) {
  # Check which genes are present in the dataset
  genes_present <- gene_list[gene_list %in% rownames(schwann_cancer)]
  print(paste("Found", length(genes_present), "out of", length(gene_list), category_name, "genes in dataset"))
  
  if(length(genes_present) == 0) {
    print(paste("No", category_name, "genes found in dataset"))
    return(NULL)
  }
  
  # If category is collagen, order numerically
  if(category_name == "Collagen") {
    # Attempt to extract numeric parts and order
    tryCatch({
      # Extract the numeric part from each collagen gene
      col_numbers <- as.numeric(str_extract(genes_present, "(?<=COL)\\d+"))
      col_subnumbers <- as.numeric(str_extract(genes_present, "(?<=A)\\d+"))
      
      # Create a data frame for sorting
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
  
  # Create dotplot with scale=FALSE to prevent scaling warning
  dot_plot <- DotPlot(schwann_cancer, 
                     features = genes_present, 
                     cols = c("lightgrey", "blue"),
                     dot.scale = 8,
                     scale = FALSE) + 
            RotatedAxis() +
            ggtitle(paste(category_name, "Genes by", gsub("_", " ", current_status))) +
            theme(plot.title = element_text(size = 14, face = "bold"),
                  axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save dotplot
  plot_width <- max(12, length(genes_present) * 0.4)  # Adjust width based on number of genes
  ggsave(paste0("pictures/COL6_", tolower(category_name), "_genes.pdf"), 
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
  
  # Determine high and low status based on current identity
  if(current_status == "COL6_status") {
    high_status <- "COL6-high"
    low_status <- "COL6-low"
  } else {
    gene_name <- sub("_status$", "", current_status)
    high_status <- paste0(gene_name, "-high")
    low_status <- paste0(gene_name, "-low")
  }
  
  for(gene in genes_present) {
    # Get expression values for high and low groups - FIXED: use @meta.data to access metadata
    high_exp <- GetAssayData(schwann_cancer, layer = "data")[gene, schwann_cancer@meta.data[[current_status]] == high_status]
    low_exp <- GetAssayData(schwann_cancer, layer = "data")[gene, schwann_cancer@meta.data[[current_status]] == low_status]
    
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
           paste0("statistics/COL6_", tolower(category_name), "_genes_statistics.csv"), 
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

print("Analyzing other collagen genes...")
collagen_results <- analyze_gene_category(other_collagen_genes, "Collagen")

#-----------------------------------------------
# 7. Create combined visualization of top genes
#-----------------------------------------------

# Get top significant genes from each category (up to 5 from each)
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
  # Define high/low status names here (outside the function)
  if(current_status == "COL6_status") {
    high_status <- "COL6-high"
    low_status <- "COL6-low"
  } else {
    gene_name <- sub("_status$", "", current_status)
    high_status <- paste0(gene_name, "-high")
    low_status <- paste0(gene_name, "-low")
  }

  # Create combined dotplot of top significant genes with scale=FALSE
  p_top <- DotPlot(schwann_cancer, 
                  features = all_top_genes, 
                  cols = c("lightgrey", "blue"),
                  dot.scale = 8,
                  scale = FALSE) + 
         RotatedAxis() +
         ggtitle(paste("Top Significant Genes by", gsub("_", " ", current_status))) +
         theme(plot.title = element_text(size = 14, face = "bold"),
               axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save combined plot
  plot_width <- max(12, length(all_top_genes) * 0.4)
  ggsave("pictures/COL6_top_significant_genes.pdf", 
         p_top, width = plot_width, height = 7)
  
  # Create heatmap of top genes - now high_status and low_status are defined
  p_heatmap <- DoHeatmap(schwann_cancer, 
                        features = all_top_genes, 
                        group.by = current_status,
                        disp.min = -2.5, 
                        disp.max = 2.5) + 
              ggtitle(paste("Top Significant Genes:", high_status, "vs", low_status))
  
  ggsave("pictures/COL6_top_genes_heatmap.pdf", 
         p_heatmap, width = 14, height = 10)
} else {
  print("No significant genes found across categories")
}

#-----------------------------------------------
# 8. Find all differentially expressed genes for each COL6 gene
#-----------------------------------------------

# Perform differential expression for each collagen VI gene
for(gene in col6_genes) {
  status_field <- paste0(gene, "_status")
  if(status_field %in% colnames(schwann_cancer@meta.data)) {
    print(paste("Finding all differentially expressed genes for", gene, "high vs low..."))
    
    # Set identity to current gene status
    Idents(schwann_cancer) <- status_field
    
    # Find markers
    high_status <- paste0(gene, "-high")
    low_status <- paste0(gene, "-low")
    
    gene_degs <- FindMarkers(schwann_cancer, 
                           ident.1 = high_status, 
                           ident.2 = low_status,
                           min.pct = 0.25)
    
    # Add gene names as column
    gene_degs$gene <- rownames(gene_degs)
    
    # Save to file
    write.csv(gene_degs, 
             paste0("statistics/", gene, "_high_vs_low_all_DEGs.csv"), 
             row.names = FALSE)
    
    print(paste("Found", nrow(gene_degs), "total differentially expressed genes for", gene))
  }
}

# For combined status, if available
if("COL6_status" %in% colnames(schwann_cancer@meta.data)) {
  print("Finding all differentially expressed genes for combined COL6 high vs low...")
  
  # Set identity to combined status
  Idents(schwann_cancer) <- "COL6_status"
  
  # Find markers
  combined_degs <- FindMarkers(schwann_cancer, 
                             ident.1 = "COL6-high", 
                             ident.2 = "COL6-low",
                             min.pct = 0.25)
  
  # Add gene names as column
  combined_degs$gene <- rownames(combined_degs)
  
  # Save all DEGs
  write.csv(combined_degs, 
           "statistics/COL6_combined_high_vs_low_all_DEGs.csv", 
           row.names = FALSE)
  
  print(paste("Found", nrow(combined_degs), "total differentially expressed genes for combined COL6 status"))
}

print("Collagen VI analysis complete!")

#-----------------------------------------------
# 9. Cell Cycle Gene Analysis
#-----------------------------------------------

# Define a list of common cell cycle genes
list_cell_cycle <- "CDK2,CDK4,CDK6,CDK7,CDKN1A,CDKN1B,STAG1,CDKN1C,CDKN2A,CDKN2B,CDKN2C,CDKN2D,ANAPC10,MAD2L2,STAG2,PTTG2,GADD45G,DBF4,YWHAQ,CHEK1,CHEK2,CREBBP,GADD45A,E2F1,E2F2,E2F3,E2F4,E2F5,EP300,ORC6,ORC3,CDC26,ABL1,ANAPC13,SMC1B,SFN,GSK3B,ANAPC2,ANAPC4,HDAC1,HDAC2,MAD2L1,SMAD2,SMAD3,SMAD4,MCM2,MCM3,MCM4,MCM5,MCM6,MCM7,MDM2,MYC,GADD45B,ATM,WEE2,ORC1,ORC2,ORC4,ORC5,PCNA,FZR1,ANAPC5,ANAPC7,ANAPC11,PLK1,ATR,PRKDC,RAD21,RB1,RBL1,RBL2,CCND1,ANAPC1,SKP1,SKP2,,,BUB1,BUB1B,TFDP1,TFDP2,TGFB1,TGFB2,TGFB3,TP53,TTK,SKP1P2,,WEE1,YWHAB,YWHAE,YWHAG,YWHAH,YWHAZ,ZBTB17,SMC1A,CDC7,CDC45,MAD1L1,CUL1,CCNB3,CDC14B,CDC14A,CDC23,CDC16,CCNA2,CCNA1,CCNB1,CCND2,CCND3,CCNE1,CCNH,PKMYT1,SMC3,CCNB2,CCNE2,BUB3,PTTG1,ESPL1,CDK1,CDC6,CDC20,CDC25A,CDC25B,CDC25C,CDC27,RBX1"
cell_cycle_genes <- unlist(strsplit(list_cell_cycle, split = ","))
cell_cycle_genes <- trimws(cell_cycle_genes)         # Remove any leading/trailing whitespace
cell_cycle_genes <- cell_cycle_genes[cell_cycle_genes != ""]  # Remove empty strings

# Check which cell cycle genes are present in the dataset
present_cc_genes <- cell_cycle_genes[cell_cycle_genes %in% rownames(schwann_cancer)]

if(length(present_cc_genes) == 0) {
  print("No cell cycle genes found in the dataset.")
} else {
  print(paste("Found", length(present_cc_genes), "cell cycle genes in the dataset."))
  
  # Create a dotplot for cell cycle genes (using the current grouping, e.g. COL6_status)
  p_cellcycle <- DotPlot(schwann_cancer, 
                         features = present_cc_genes, 
                         cols = c("lightgrey", "blue"),
                         dot.scale = 8,
                         scale = FALSE) + 
                 RotatedAxis() +
                 ggtitle(paste("Cell Cycle Genes by", gsub("_", " ", current_status))) +
                 theme(plot.title = element_text(size = 14, face = "bold"),
                       axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the cell cycle dotplot
  plot_width <- max(12, length(present_cc_genes) * 0.4)
  ggsave(filename = "pictures/COL6_cell_cycle_genes.pdf", 
         plot = p_cellcycle, width = plot_width, height = 7)
  
  # Create a heatmap for cell cycle genes using the current grouping
  p_cc_heatmap <- DoHeatmap(schwann_cancer, 
                            features = present_cc_genes, 
                            group.by = current_status,
                            disp.min = -2.5, 
                            disp.max = 2.5) + 
                  ggtitle("Cell Cycle Genes Heatmap") +
                  theme(plot.title = element_text(size = 14, face = "bold"))
  
  # Save the cell cycle heatmap
  ggsave(filename = "pictures/COL6_cell_cycle_genes_heatmap.pdf", 
         plot = p_cc_heatmap, width = 14, height = 10)
}

#----- Compute TRAP1 high/low grouping -----
# Make sure TRAP1 is present
if("TRAP1" %in% rownames(schwann_cancer)) {
  TRAP1_exp <- GetAssayData(schwann_cancer, assay = "RNA", layer = "data")["TRAP1", ]
  TRAP1_cutoff <- median(TRAP1_exp)
  schwann_cancer$TRAP1_status <- ifelse(TRAP1_exp > TRAP1_cutoff, "TRAP1-high", "TRAP1-low")
  print(paste("TRAP1 cutoff:", round(TRAP1_cutoff, 4)))
} else {
  stop("TRAP1 not found in the dataset!")
}

#----- Cell Cycle Gene Analysis based on TRAP1 high vs low -----
# Define the list of common cell cycle genes
list_cell_cycle <- "CDK2,CDK4,CDK6,CDK7,CDKN1A,CDKN1B,STAG1,CDKN1C,CDKN2A,CDKN2B,CDKN2C,CDKN2D,ANAPC10,MAD2L2,STAG2,PTTG2,GADD45G,DBF4,YWHAQ,CHEK1,CHEK2,CREBBP,GADD45A,E2F1,E2F2,E2F3,E2F4,E2F5,EP300,ORC6,ORC3,CDC26,ABL1,ANAPC13,SMC1B,SFN,GSK3B,ANAPC2,ANAPC4,HDAC1,HDAC2,MAD2L1,SMAD2,SMAD3,SMAD4,MCM2,MCM3,MCM4,MCM5,MCM6,MCM7,MDM2,MYC,GADD45B,ATM,WEE2,ORC1,ORC2,ORC4,ORC5,PCNA,FZR1,ANAPC5,ANAPC7,ANAPC11,PLK1,ATR,PRKDC,RAD21,RB1,RBL1,RBL2,CCND1,ANAPC1,SKP1,SKP2,,,BUB1,BUB1B,TFDP1,TFDP2,TGFB1,TGFB2,TGFB3,TP53,TTK,SKP1P2,,WEE1,YWHAB,YWHAE,YWHAG,YWHAH,YWHAZ,ZBTB17,SMC1A,CDC7,CDC45,MAD1L1,CUL1,CCNB3,CDC14B,CDC14A,CDC23,CDC16,CCNA2,CCNA1,CCNB1,CCND2,CCND3,CCNE1,CCNH,PKMYT1,SMC3,CCNB2,CCNE2,BUB3,PTTG1,ESPL1,CDK1,CDC6,CDC20,CDC25A,CDC25B,CDC25C,CDC27,RBX1"
cell_cycle_genes <- unlist(strsplit(list_cell_cycle, split = ","))
cell_cycle_genes <- trimws(cell_cycle_genes)
cell_cycle_genes <- cell_cycle_genes[cell_cycle_genes != ""]

# Check which cell cycle genes exist in your dataset
present_cc_genes <- cell_cycle_genes[cell_cycle_genes %in% rownames(schwann_cancer)]
if(length(present_cc_genes) == 0) {
  print("No cell cycle genes found in the dataset.")
} else {
  print(paste("Found", length(present_cc_genes), "cell cycle genes in the dataset for TRAP1-based analysis."))
  
  # Create a dotplot for cell cycle genes grouped by TRAP1_status
  p_cellcycle_TRAP1 <- DotPlot(schwann_cancer, 
                               features = present_cc_genes, 
                               group.by = "TRAP1_status",
                               cols = c("lightgrey", "blue"),
                               dot.scale = 8,
                               scale = FALSE) + 
                      RotatedAxis() +
                      ggtitle("Cell Cycle Genes by TRAP1 Status") +
                      theme(plot.title = element_text(size = 14, face = "bold"),
                            axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the dotplot
  plot_width <- max(12, length(present_cc_genes) * 0.4)
  ggsave(filename = "pictures/TRAP1_cell_cycle_genes_dotplot.pdf", 
         plot = p_cellcycle_TRAP1, width = plot_width, height = 7)
  
  # Create a heatmap for cell cycle genes using TRAP1_status grouping
  p_cc_heatmap_TRAP1 <- DoHeatmap(schwann_cancer, 
                                  features = present_cc_genes, 
                                  group.by = "TRAP1_status",
                                  disp.min = -2.5, 
                                  disp.max = 2.5) + 
                        ggtitle("Cell Cycle Genes Heatmap (TRAP1 high vs low)") +
                        theme(plot.title = element_text(size = 14, face = "bold"))
  
  # Save the heatmap
  ggsave(filename = "pictures/TRAP1_cell_cycle_genes_heatmap.pdf", 
         plot = p_cc_heatmap_TRAP1, width = 14, height = 10)
}

# Ensure the active identities are set to COL6_status
Idents(schwann_cancer) <- "COL6_status"
print(unique(Idents(schwann_cancer)))  # Should print "COL6-high" and "COL6-low"

# DEG Analysis for COL6 high vs low
deg_COL6_high_vs_low <- FindMarkers(schwann_cancer, 
                                    ident.1 = "COL6-high", 
                                    ident.2 = "COL6-low", 
                                    min.pct = 0.25)
deg_COL6_high_vs_low$gene <- rownames(deg_COL6_high_vs_low)

# Remove the COL6 genes from the DEG results
col6_genes_to_remove <- c("COL6A1", "COL6A2", "COL6A3")
deg_COL6_no_COL6 <- deg_COL6_high_vs_low[!rownames(deg_COL6_high_vs_low) %in% col6_genes_to_remove, ]
deg_COL6_no_COL6$gene <- rownames(deg_COL6_no_COL6)

# Save DEG table without COL6 genes
write.csv(deg_COL6_no_COL6, "statistics/COL6_high_vs_low_DEGs_no_COL6_schwann_cancer.csv", row.names = FALSE)

# Filter for significant DEGs (adjusted p-value < 0.05)
sig_deg_COL6_no_COL6 <- subset(deg_COL6_no_COL6, p_val_adj < 0.05)

# Order by absolute log2 Fold Change and take the top 100 genes
top100_COL6_no_COL6 <- head(sig_deg_COL6_no_COL6[order(-abs(sig_deg_COL6_no_COL6$avg_log2FC)), ], 100)

# Create volcano plot (plot all points and label only the top 100 significant genes)

p_volcano_COL6 <- ggplot(deg_COL6_no_COL6, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  geom_point(data = sig_deg_COL6_no_COL6, color = "red", alpha = 0.7) +
  geom_text_repel(data = top100_COL6_no_COL6, 
                  aes(label = gene),
                  size = 2.5,                    # slightly smaller label size
                  box.padding   = unit(0.5, "lines"),
                  max.overlaps = 100,
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',
                  force = 2) +                    # increases repelling force
  theme_minimal() +
  ggtitle("COL6 DEG Volcano Plot (High vs Low) without COL6A1/2/3") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value")
  
ggsave("pictures/COL6_DEG_volcano_top100_no_COL6_schwann_cancer.pdf", 
       p_volcano_COL6, width = 8, height = 7)

# Save the volcano plot as a PDF
ggsave("pictures/COL6_DEG_volcano_top100_no_COL6_schwann_cancer.pdf", p_volcano_COL6, width = 8, height = 7)