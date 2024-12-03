#! /usr/bin/env python3
#Import required packages and modules.
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np
import random
import seaborn as sns
import tqdm as notebook_tqdm

#Set variables and settings.
sc.settings.set_figure_params(dpi=50, facecolor="white")
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17'
os.chdir(directory_path)
save_name = "NEMP17"

#Load AnnData object.
anndata = directory_path + '/NEMP17.h5ad'
adata = sc.read_h5ad(anndata)

#Subset dataset in thirds
third = adata.n_obs // 3
adata_sub1 = adata[:third, :]
adata_sub2 = adata[third:2*third, :]
adata_sub3 = adata[2*third:, :]

#Set counts layer to adata.X
adata_sub1.X = adata_sub1.layers['counts'].copy()
adata_sub2.X = adata_sub2.layers['counts'].copy()
adata_sub3.X = adata_sub3.layers['counts'].copy()

#Save the subsets as new AnnData objects
adata_sub1.write("adata_sub1.h5ad", compression='gzip')
adata_sub2.write("adata_sub2.h5ad", compression='gzip')
adata_sub3.write("adata_sub3.h5ad", compression='gzip')
