# Vibe Coding

*How I learned to stop worrying and talk to computer*

Objective: Use gemini to help set up a computational project without knowing how to code.

**Learning goals:**

- ask gemini how to do everything
- use the terminal
- use project folder with git
- install R and run analysis in seurat

## My queries to gemini

1. I have single cell RNA seq counts data that I want to filter, qc, and analyze in R studio. how do I download R studio onto my mac so I can use the latest version of seurat for analysis?

2. how do i set up a project folder with a readme, git repo, R code, data and keep figures and analysis organized? I want to use the terminal to set up my folder structure and I need to download the single cell data from the server, also using the terminal. I connected to my server and its mounted as a volumes/

3. (specific error) my r studio says: RStudio can't establish a connection to R. but when i command space i see r and version 4.5 opens fine.

4. how do i move r studio to the ideal location and how do i make a project in my new seurat_project folder

5. can you explain how to run the R code? also I got an error when I did run it. Error in library(Seurat) : there is no package called ‘Seurat’. Finally I would like to load my own data. I have 4 samples in subfolders data/raw/<sample>/DGE_unfiltered/count_matrix.mtx.gz

6. I have a all_genes.csv.gz and a cell_metadata.csv.gz. will that work? also I got a error when loading tidyverse: ── Conflicts ── tidyverse_conflicts() ──

✖ dplyr::filter() masks stats::filter()

7. the code you gave me looks different from the intial code, here: # Create the Seurat Object

obj <- CreateSeuratObject(counts = counts, project = sample)

return(obj) Do i need to return the obj/pbmc??? what is lapply and function doijng? it looks different from a simple loop?

8. can you rewrite loading using the R way and the csv loading?

9. I got this error: Error in fixupDN.if.valid(value, x@Dim) : 

  length of Dimnames[[1]] (119765) is not equal to Dim[1] (147432)
  
10. so does counts data already have rownames? I assumed it was unlabeled and that was why i needed row/gene names. how can I "look" at counts? its shape and row and column names?

11. dim(counts)

[1] 147432 119764 this is a transposing issue

12. what about the other metadata? should I put that in my seurat object for QC or normaliztion purposes? these are the columns in meta:   bc_wells     sample  species gene_count tscp_count mread_count GRCh38_gene_count GRCh38_tscp_count GRCh38_mread_count

13. since the first column of meta is bc_wells, are you sure colnames(counts) <- meta[[1]]   is correct? how do i check if these are unique cell IDs or repeated well IDs?

14. okay, I loaded the 4 samples into combined_pbmc. I have about 500K cells per sample, except day zero I have 150K. we had previously tried filtering the cells by number of genes or reads ( we used an "elbow plot"). but we ended up filtering out all of the immune cells because they are smaller and had less total DNA than the other cell types (epithelial, fibroblast, stromal and endothelial). what is the standard QC pipeline and does it account for cell size? also, just FYI, these are not 10x droplet reads but they are from PARSE bioscience which uses combinatorial barcoding and may have special QC considerations.

15. I would like to QC each sample and save the filtered data. then maybe it will be small enough to combine for UMAP. at that point I can re-evaluate "sketching"... please remind me.

16. yes, swap those in. also should I be saving each version of the code to git? i never set up my remote repo. also, my earlier plot of nCount versus nFeature_RNA (what exactly are those) showed D0 had higher nCount (1e5) and nFeature (10K). D0 was the 2D cell lines whereas the otehr days were 3D bioprinted tissues which might explain the difference. I worry that the 20K ceiling does not apply to D0 for doublets.

17. i want to save in a format interoperable with python, not .rds. I want to save plots showing the data and thresholds applied for each sample.

18. got these errors: Warning: Data is of class dgTMatrix. Coercing to dgCMatrix.
Rasterizing points since number of points exceeds 100,000.
To disable this behavior set `raster=FALSE`
Warning: Default search for "data" layer in "RNA" assay yielded no results; utilizing "counts" layer instead.
Error in (function (cond)  : 
  error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject 5.0.0 and is
now defunct.
ℹ Please use the `layer` argument instead.
In addition: Warning messages:
1: In SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents = idents,  :
  All cells have the same value of percent.mt.
2: Removed 1 row containing missing values or values outside the scale range
(`geom_hline()`). 

19.for git, i want to add figures/ to my .gitignore and then add my updated readme.md and code

20. why do my violin plots have more than one x-value "identity" for some samples? D0 has one but the others days have 4, e.g. 58, 59, 60, 61 for my D9 untreated.

21. yes, the ridge plot would be good because the violin plots are overplotted with the black stripplot on top. I want to keep seeing the sub-pools in QC violin/ridge plots of percent MT

22. lets put those suggestions on pause and summarize the qc and filtering that we have done into a document that explains the project set up, git, r installl and project setup, data download with rsync, code steps (specific to PARSE) for non-computational biologists so they can do this on their own.


## Terminal linux commands

- `ls` list files
- `ls -la` list hidden files like .gitignore
- `\tab` to auto complete
- `cd` change directory
- `cd ..` go up one directory
- `mkdir` make directory
- `rm` remove/delete a file
- `rm -r` remove/delete a folder
- `cp` copy a file
- `rsync -avP /Volumes/BCCResearch/danielc/BCC\ Bioprint\ Parse\ data/output_combined/7305_D5__UNTR data/raw/` example rsync command.  *careful no trailing </> after the server path!*

# Single-cell RNA-seq Analysis
**Project Start Date:** 2026-01-19
**Data Source:** PARPi treated bioprint. /Volumes/BCCResearch/danielc/BCC\ Bioprint\ Parse\ data/output_combined/7305_*
**Environment:** Mention you are using Seurat v5 and R 4.5.
**The Task:** Analyzing 1.6M cells from Parse Evercode across Day 0 2D and 3D bioprinted tissues.
**Goal:** Filter and cluster single cells from 4 samples to identify cell populations. 
Then look at cell population and tumor pathway activity changes between treatment conditions.
**Hypothesis:** PARPi reduces MYC pathway activity in treated tumor cells because cells with 
high MYC have more damage and can't survive PARPi.
**Reproducibility:** The script 01_qc_and_filtering.R generates the figures/ and the .h5ad files.

## Technical QC Observations
- **Sub-libraries:** Identified via Parse Pipeline metadata. 
- **Visualization:** Used RidgePlots instead of standard Violin plots to avoid point-overplotting at the >1M cell scale.
- **Thresholds:** Adaptive filtering was applied per sample (D0 vs 3D Tissues).
