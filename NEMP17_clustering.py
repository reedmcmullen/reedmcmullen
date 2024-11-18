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
import random

#Define and set the current working directory.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17'
os.chdir(directory_path)

#Read in full dataset as and AnnData object.
print('Loading dataset...')
anndata = directory_path + '/NEMP17.h5ad'
adata = sc.read_h5ad(anndata)

#Redo KNN and UMAP for full integrated dataset.
print('Running neighbor finding...')
sc.pp.neighbors(adata)
print('Running UMAP...')
sc.tl.umap(adata)

#Redo clustering for full integrated dataset
print('Running leiden clustering...')
res = 0.75
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=res)

#Save dataset
print('Saving dataset...')
anndata = directory_path + '/NEMP17.h5ad'
adata.write(anndata, compression='gzip')

#Subset dataset to 100k cells.
random.seed(0)
adata_sub=adata[np.random.choice(adata.obs.index, 100000, replace=False),:]

#Save dataset subset.
print('Saving dataset subset...')
anndata_subset = directory_path + '/NEMP17_subset.h5ad'
adata_sub.write(anndata_subset)


