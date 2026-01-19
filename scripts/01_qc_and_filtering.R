library(Seurat)
library(vroom)
library(anndata)
library(ggplot2)

# Define your adaptive thresholds again
thresholds <- list(
  "7305_D0"    = list(min_feat = 500, max_feat = 12000, max_mt = 10, max_count = 150000),
  "7305_D5__UNTR" = list(min_feat = 200, max_feat = 6000,  max_mt = 15, max_count = 30000),
  "7305_D9_OLAP" = list(min_feat = 200, max_feat = 6000,  max_mt = 15, max_count = 30000),
  "7305_D9_UNTR" = list(min_feat = 200, max_feat = 6000,  max_mt = 15, max_count = 30000)
)

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
  
  rownames(counts) <- make.unique(genes[[2]])
  colnames(counts) <- meta$bc_wells
  
  obj <- CreateSeuratObject(counts = counts, project = s_name, meta.data = as.data.frame(meta))
  
  # --- CRITICAL MITO CHECK ---
  # Let's look at the first few gene names to see why MT- might be failing
  message("First few gene names: ", paste(head(rownames(obj), 3), collapse=", "))
  
  # Try both common patterns: ^MT- (Human) and ^mt- (Mouse)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  
  # 2. GENERATE AND SAVE QC PLOTS
  # We use layer="counts" for Seurat v5
  p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_vline(xintercept = conf$max_count, color="red", linetype="dashed") +
    geom_hline(yintercept = c(conf$min_feat, conf$max_feat), color="red", linetype="dashed") +
    theme_bw()
  
  # A much cleaner Violin plot for 1M+ cells
  p2 <- VlnPlot(obj, features = "percent.mt", pt.size = 0) + 
    geom_hline(yintercept = conf$max_mt, color="red", linetype="dashed")
  
  ggsave(paste0("figures/", s_name, "_qc_scatter.png"), p1, width = 8, height = 6)
  ggsave(paste0("figures/", s_name, "_mito_vln.png"), p2, width = 6, height = 6)
  
  # 1. Ridge Plot for Mitochondrial % (The most important QC check)
  # we use 'idents' to show the sub-pools automatically
  p3 <- RidgePlot(obj, features = "percent.mt", ncol = 1) +
    geom_vline(xintercept = conf$max_mt, color="red", linetype="dashed") +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle(paste(s_name, "Mito % by Sub-pool"))
  
  # 2. Ridge Plot for nFeature (Complexity)
  p4 <- RidgePlot(obj, features = "nFeature_RNA", ncol = 1) +
    scale_x_log10() + # Log scale helps see the immune 'bump' better
    theme_minimal() +
    ggtitle(paste(s_name, "Gene Complexity by Sub-pool"))
  
  # Save these new views
  ggsave(paste0("figures/", s_name, "_mito_ridge.png"), p3, width = 7, height = 8)
  ggsave(paste0("figures/", s_name, "_feature_ridge.png"), p4, width = 7, height = 8)
  
  # 3. FILTER
  obj <- subset(obj, subset = nFeature_RNA > conf$min_feat & 
                  nFeature_RNA < conf$max_feat & 
                  nCount_RNA < conf$max_count &
                  percent.mt < conf$max_mt)
  
  # 4. SAVE AS H5AD (FOR PYTHON) - Seurat v5 Layer Syntax
  # We pull the 'counts' layer specifically
  count_layer <- LayerData(obj, layer = "counts")
  
  adata <- AnnData(
    X = t(count_layer), # AnnData still expects Cells x Genes
    obs = obj@meta.data,
    var = data.frame(gene_names = rownames(obj), row.names = rownames(obj))
  )
  adata$write_h5ad(paste0("data/processed/", s_name, "_filtered.h5ad"))
  
  rm(obj, adata, counts, count_layer, genes, meta); gc()
}
