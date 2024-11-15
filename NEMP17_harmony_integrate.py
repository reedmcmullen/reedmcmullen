#! /usr/bin/env python3
#Import required packages and modules.
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np

#Define and set the current working directory.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17'
os.chdir(directory_path)

#Read in concatenated, preprocessed data from a H5AD file to an AnnData object.
anndata = directory_path + '/NEMP17.h5ad'
adata = sc.read_h5ad(anndata)

#Integrate full dataset with harmony.
sce.pp.harmony_integrate(adata, key=['GEMwell', 'species_assignment'], basis='X_pca', adjusted_basis='X_pca_harmony')

#Set harmony PCA values to be the default PCA values.
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

anndata = directory_path + '/NEMP17.h5ad'
adata.write(anndata)
adata
