
#Install and load remotes packages.
install.packages('remotes')
library(remotes)

#Install and load required packages.
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("zellkonverter", version = "1.14.1")
library(Seurat)
library(zellkonverter)

#Read in AnnData object from H5AD file, with raw count data in adata.X.
h5ad_file <- "/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/dreamlet_NEMP17/NEMP17_counts.h5ad"
adata <- readH5AD(h5ad_file, X_name='counts')

#Convert to RDS file
rds_file <- "/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/dreamlet_NEMP17/NEMP17_dreamlet.rds"
saveRDS(adata, file=rds_file)

sessionInfo()
