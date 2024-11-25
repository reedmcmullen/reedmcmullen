#!/usr/bin/env Rscript
#Before running script, install the required packages.
#install.packages('remotes')
#library(remotes)
#remotes::install_github("satijalab/seurat", "seurat5")
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("zellkonverter")

#Load the required packages.
library(Seurat)
library(zellkonverter)
library(reticulate)

#Read in the AnnData object from H5AD file, with raw count data in adata.X.
h5ad_file <- "/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/dreamlet_NEMP17/NEMP17_counts.h5ad"
adata <- readH5AD(h5ad_file, X_name='counts', reader='R')

#Convert to RDS file
rds_file <- "/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/dreamlet_NEMP17/NEMP17_dreamlet.rds"
saveRDS(adata, file=rds_file)

sessionInfo()
