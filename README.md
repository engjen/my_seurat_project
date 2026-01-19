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
**Goal:** Filter and cluster samples to identify cell populations. Then look at cell population and tumor pathway activity changes between treatment conditions.
**Hypothesis:** PARPi reduces MYC pathway activity in treated tumor cells because cells with high MYC have more damage and can't survive PARPi.
