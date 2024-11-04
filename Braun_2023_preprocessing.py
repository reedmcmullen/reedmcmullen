#Load in required packages and modules.
import os
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch

#Set settings.
sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")

#Define and set the current working directory.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/refmap_NEMP17'
os.chdir(directory_path)

#Load in reference dataset as an AnnData object from a H5AD file.
adata_ref = sc.read_10x_h5('/wynton/scratch/rmcmullen/datasets/HumanFetalBrainPool.h5')
adata_ref

#Subset the AnnData object to 20% of cells.
# Set a seed for reproducibility.
np.random.seed(42)
# Calculate 20% of the total number of cells
n_cells = int(adata_ref.n_obs * 0.2)
# Randomly sample indices without replacement
random_indices = np.random.choice(adata_ref.obs.index, size=n_cells, replace=False)
# Subset the AnnData object
adata_ref = adata_ref[random_indices].copy()
adata_ref

#Save subset reference
ref_initial = directory_path + '/ref_initial.h5ad'
adata_ref.write(ref_initial)
