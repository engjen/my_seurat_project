library(Seurat)
library(vroom)
library(anndata)
library(ggplot2)
library(patchwork)
library(ggplot2)

# Define your adaptive thresholds again
thresholds <- list(
  "7305_D0"    = list(min_feat = 20, max_feat = 12000, max_mt = 15, max_count = 150000),
  "7305_D5__UNTR" = list(min_feat = 20, max_feat = 8000,  max_mt = 15, max_count = 30000),
  "7305_D9_OLAP" = list(min_feat = 20, max_feat = 8000,  max_mt = 15, max_count = 30000),
  "7305_D9_UNTR" = list(min_feat = 20, max_feat = 8000,  max_mt = 15, max_count = 30000)
)
min_genes = 3
for (s_name in names(thresholds)) {
  message("Processing: ", s_name)
  conf <- thresholds[[s_name]]
  
  # 1. LOAD DATA
  dir_path <- paste0("data/raw/", s_name, "/DGE_unfiltered/")
  genes  <- vroom(paste0(dir_path, "all_genes.csv.gz"), col_names = TRUE, show_col_types = FALSE)
  meta   <- vroom(paste0(dir_path, "cell_metadata.csv.gz"), col_names = TRUE, show_col_types = FALSE)
  counts <- readMM(paste0(dir_path, "count_matrix.mtx.gz"))
  
  # Coerce to dgCMatrix immediately (standardizes format and saves memory)
  counts <- as(t(counts), "CsparseMatrix")
  
  colnames(counts) <- meta$bc_wells
  # 1. Clean names IMMEDIATELY after defining them
  rownames(counts) <- gsub("_", "-", make.unique(genes[[2]]))
  
  obj <- CreateSeuratObject(counts = counts, project = s_name, meta.data = as.data.frame(meta))
  
  # Try both common patterns: ^MT- (Human) and ^mt- (Mouse)
  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  
  # knee plot
  # 1. Extract the UMI counts per cell (barcode)
  counts_per_cell <- obj$nCount_RNA
  
  # 2. Sort them in descending order
  df_knee <- data.frame(
    rank = 1:length(counts_per_cell),
    counts = sort(counts_per_cell, decreasing = TRUE)
  )
  
  # 3. Filter out absolute zeros (to avoid log10 infinity issues)
  df_knee <- df_knee[df_knee$counts > 0, ]
  
  k1 <- ggplot(df_knee, aes(x = rank, y = counts)) +
    geom_line(size = 1, color = "steelblue") +
    scale_x_log10() + 
    scale_y_log10() +
    geom_hline(yintercept = 200, color = "red", linetype = "dashed") + # Your min_feat floor
    labs(
      title = "Knee Plot: Cell Barcode Quality",
      x = "Barcode Rank (Log10)",
      y = "UMI Counts per Barcode (Log10)"
    ) +
    theme_minimal()
  ggsave(paste0("figures/",  "qc_knee_plot",s_name,".png"), k1, width = 6, height = 6)
  # 2. GENERATE AND SAVE QC PLOTS
  # We use layer="counts" for Seurat v5
  p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = TRUE) +
    geom_vline(xintercept = conf$max_count, color="red", linetype="dashed") +
    geom_hline(yintercept = c(conf$min_feat, conf$max_feat), color="blue", linetype="dashed") +
    theme_bw()
  
  # A much cleaner Violin plot for 1M+ cells
  p2 <- VlnPlot(obj, features = "percent_mt", pt.size = 0, layer = "counts") + 
    geom_hline(yintercept = conf$max_mt, color="red", linetype="dashed")
  
  # Creating a density-aware scatter plot
  p3 <- ggplot(obj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_bin2d(bins = 100) + # Creates a heat-map style density
    scale_fill_viridis_c() +
    scale_x_log10() + scale_y_log10() +
    geom_vline(xintercept = conf$max_count, color="red", linetype="dashed") +
    geom_hline(yintercept = c(conf$min_feat, conf$max_feat), color="blue", linetype="dashed") +
    labs(title = paste(s_name, "Density Scatter"), fill = "Cell Density") +
    theme_minimal()

  # We set raster = FALSE to silence the raster warning (since we aren't plotting points anyway)
  # We set pt.size = 0 to ensure no points are drawn
  v1 <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0,  layer = "counts") + 
    geom_hline(yintercept = c(conf$min_feat, conf$max_feat), color="red", linetype="dashed") +
    NoLegend() + theme(axis.title.x=element_blank())
  
  v2 <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0, layer = "counts") + 
    geom_hline(yintercept = conf$max_count, color="red", linetype="dashed") +
    NoLegend() + theme(axis.title.x=element_blank())
  
  v3 <- VlnPlot(obj, features = "percent_mt", pt.size = 0, layer = "counts") + 
    geom_hline(yintercept = conf$max_mt, color="red", linetype="dashed") +
    NoLegend() + theme(axis.title.x=element_blank())
  
  # --- FORCE HORIZONTAL LAYOUT ---
  # use plot_layout(ncol = 3) to force them into a single row
  summary_vln <- v1 + v2 + v3 + 
    plot_layout(ncol = 3) + 
    plot_annotation(title = paste("Aggregated QC Summary:", s_name))
    # 3. Save it
    ggsave(paste0("figures/", "qc_summary_violins", s_name,".png"), summary_vln, width = 12, height = 5)
    ggsave(paste0("figures/", "qc_scatter", s_name,".png"), p1, width = 8, height = 6)
    ggsave(paste0("figures/", "qc_mito_vln", s_name,".png"), p2, width = 6, height = 6)
    ggsave(paste0("figures/",  "qc_complexity_contour",s_name,".png"), p3, width = 6, height = 6)
    
    # Calculate how many cells each gene appears in
    gene_counts <- Matrix::rowSums(LayerData(obj, layer = "counts") > 0)
    
    # Create a dataframe for plotting
    df_genes <- data.frame(cells_per_gene = gene_counts)
    
    # Plot
    p4 <- ggplot(df_genes, aes(x = cells_per_gene + 1)) +
      geom_histogram(bins = 100, fill = "steelblue", color = "white") +
      scale_x_log10() + 
      geom_vline(xintercept = min_genes, color = "red", linetype = "dashed") +
      labs(title = "Gene Detection Frequency",
           x = "Number of Cells (Log10 +1)",
           y = "Number of Genes") +
      theme_minimal()
    ggsave(paste0("figures/",  "qc_genes_cells",s_name,".png"), p4, width = 6, height = 6)
    
    # --- SNAPSHOT: BEFORE ---
    n_cells_before <- ncol(obj)
    n_genes_before <- nrow(obj)
  
  # 3. FILTER by max count, max genes, max percent mito
  obj <- subset(obj, subset = nFeature_RNA > conf$min_feat & 
                  nFeature_RNA < conf$max_feat & 
                  nCount_RNA < conf$max_count &
                  percent_mt < conf$max_mt)
  
  # Filter genes: Keep only genes expressed in at least 10 cells
  # We do this by calculating the rowSums of the counts matrix
  selected_genes <- rownames(obj)[Matrix::rowSums(LayerData(obj, layer = "counts") > 0) >= min_genes]
  obj <- subset(obj, features = selected_genes)
  
  # --- SNAPSHOT: AFTER ---
  n_cells_after <- ncol(obj)
  n_genes_after <- nrow(obj)
  
  # --- PRINT SUMMARY ---
  message(paste0(
    "\n--- QC Summary for ", s_name, " ---\n",
    "Cells: ", n_cells_before, " -> ", n_cells_after, 
    " (Removed: ", n_cells_before - n_cells_after, ")\n",
    " (Fraction removed: ", 1 - n_cells_after/n_cells_before, ")\n",
    "Genes: ", n_genes_before, " -> ", n_genes_after, 
    " (Removed: ", n_genes_before - n_genes_after, ")\n",
    " (Fraction removed: ", 1 - n_genes_after/n_genes_before, ")\n",
    "-----------------------------------\n"
  ))
  
  # List of your essential markers (using dashes to match our cleaned names)
  core_markers <- c("PTPRC", "CD3E", "EPCAM", "VWF", "ACTA2")
  
  # Check which are still in the object
  present <- core_markers %in% rownames(obj)
  names(present) <- core_markers
  
  # Print the status and cell counts for each
  message("--- Marker Survival Check ---")
  for(gene in core_markers) {
    if(gene %in% rownames(obj)) {
      # Count how many cells have > 0 counts for this gene
      n_cells <- sum(LayerData(obj, layer = "counts")[gene, ] > 0)
      cat(gene, ": FOUND in", n_cells, "cells\n")
    } else {
      cat(gene, ": MISSING (Filtered out)\n")
    }
  }
  
  # 4. SAVE AS H5AD (FOR PYTHON) - Seurat v5 Layer Syntax
  # We pull the 'counts' layer specifically
  count_layer <- LayerData(obj, layer = "counts")
  
  adata <- AnnData(
    X = t(count_layer), # AnnData still expects Cells x Genes
    obs = obj@meta.data,
    var = data.frame(gene_names = rownames(obj), row.names = rownames(obj))
  )
  adata$write_h5ad(paste0("data/processed/", s_name, "_filtered.h5ad"))
  
  rm(obj, adata, counts, count_layer, genes, meta, p1, p2, p3); gc()
  break
}
  