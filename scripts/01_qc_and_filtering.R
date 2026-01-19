# Load Libraries
library(Seurat)
library(tidyverse) # Helpful for data manipulation and plotting
library(vroom)
library(Matrix)

# Define Paths (Relative to your project root)
raw_data_path <- "data/raw/"
processed_data_path <- "data/processed/"


# 1. Define specific thresholds for your different biology
# 2D cell lines (D0) get higher ceilings; 3D tissues get lower ones
thresholds <- list(
  "7305_D0" = list(min_feat = 500, max_feat = 12000, max_mt = 10, max_count = 150000),
  "7305_D5__UNTR" = list(min_feat = 200, max_feat = 6000,  max_mt = 15, max_count = 30000),
  "7305_D9_OLAP" = list(min_feat = 200, max_feat = 6000,  max_mt = 15, max_count = 30000),
  "7305_D9_UNTR" = list(min_feat = 200, max_feat = 6000,  max_mt = 15, max_count = 30000)
)
######################## Load data #################################
#Loading a Parse matrix
for (s_name in names(thresholds)) {
  message("Processing: ", s_name)
  conf <- thresholds[[s_name]]
  
  # vroom is significantly faster for large compressed CSVs
  genes  <- vroom(paste0(dir_path, "all_genes.csv.gz"), col_names = TRUE, show_col_types = FALSE)
  meta   <- vroom(paste0(dir_path, "cell_metadata.csv.gz"), col_names = TRUE, show_col_types = FALSE)
  counts <- readMM(paste0(dir_path, "count_matrix.mtx.gz"))
  
  counts <- t(counts)
  rownames(counts) <- make.unique(genes[[1]])
  colnames(counts) <- meta$bc_wells
  
  obj <- CreateSeuratObject(counts = counts, project = s_name,
                            meta.data = as.data.frame(meta))
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Apply adaptive filtering based on our list
  obj <- subset(obj, subset = nFeature_RNA > conf$min_feat & 
                  nFeature_RNA < conf$max_feat & 
                  nCount_RNA < conf$max_count &
                  percent.mt < conf$max_mt)
  
  saveRDS(obj, paste0("data/processed/", s_name, "_filtered.rds"))
  rm(obj, counts); gc()
}