#! /usr/bin/env python3
#Import required packages and modules.
#Import required packages and modules.
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np

#Define and set the current working directory.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17'
os.chdir(directory_path)

#Read in full dataset as and AnnData object.
anndata = directory_path + '/NEMP17.h5ad'
adata = sc.read_h5ad(anndata)

#KNN, UMAP, and leiden clustering.
sc.pp.neighbors(adata)
sc.tl.umap(adata)
res = 0.5
sc.tl.leiden(adata, resolution=res)

#Save dataset.
anndata = directory_path + '/NEMP17.h5ad'
adata.write(anndata, compression='gzip')

#Subset the AnnData object to 10% of cells for faster testing and visualization.
np.random.seed(42)
n_cells = int(adata.n_obs * 0.1)
random_indices = np.random.choice(adata.obs.index, size=n_cells, replace=False)
adata_sub = adata[random_indices].copy()

#Save subset dataset.
anndata_subset = directory_path + '/NEMP17_subset.h5ad'
adata_sub.write(anndata_subset, compression='gzip')

