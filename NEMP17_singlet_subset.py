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

#Define string for use in naming saved figures.
save_name = "NEMP17"

#Read in concatenated, preprocessed data from a H5AD file to an AnnData object.
anndata_initial = directory_path + '/NEMP17_preprocessed.h5ad'
adata = sc.read_h5ad(anndata_initial)

#Subset to singlets.
adata = adata[adata.obs['species_droplet_type']=='S']

#Save the AnnData object as an .h5ad file.
anndata = directory_path + '/NEMP17.h5ad'
adata.write(anndata, compression='gzip')
