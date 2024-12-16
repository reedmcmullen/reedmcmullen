#! /usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import scanpy.external as sce
import anndata as ad

#Define settings and variables.
sc.settings.set_figure_params(dpi=50, facecolor="white", figsize=(6,6))
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17/rpc_search'
os.chdir(directory_path)

#Load and subset dataset.
print('Loading and subsetting dataset...')
adata = sc.read_h5ad('/wynton/group/pollen/jding/Sara/linnarsson/human_dev_GRCh38-3.0.0.h5ad')

for col in ['CellClass', 'Region', 'Subdivision', 'Subregion', 'Tissue', 'donor_id', 'dissection', 'sample_id',]:
    adata.obs[col] = adata.obs[col].str[2:-1]
adata.obs['CellID'] = adata.obs.index
adata.obs['CellID'] = adata.obs['CellID'].str[2:-1]
adata.obs.index = adata.obs['CellID']
del adata.obs['CellID']

adata = adata[~adata.obs['CellClass'].isin(['Erythrocyte', 'Fibroblast', 'Immune', 'Vascular'])]
adata = adata[adata.obs['Region'].isin(['Forebrain', 'Telencephalon'])]
adata = adata[adata.obs['Subregion'].isin(['Forebrain', 'Telencephalon', 'Cortex', 'Subcortex', 'Striatum', 'Hippocampus'])]

adata.var.index = adata.var['Gene']
adata.raw = adata
adata.layers["counts"] = adata.X.copy()

adata.obs['cluster_id'] = adata.obs['cluster_id'].astype(str)
adata.obs['TopLevelCluster'] = adata.obs['TopLevelCluster'].astype(str)

# Load in cluster metadata.
cluster_id_df = pd.read_excel('/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/refmap_NEMP17/science.adf1226_table_s2.xlsx', index_col='PoolOrder')
cluster_id_df = cluster_id_df[cluster_id_df['ClusterID (PoolClean)'] != '--']
cluster_id_df['ClusterID (PoolClean)'] = cluster_id_df['ClusterID (PoolClean)'].astype(str)

# Create mappings for AutoAnnotation, MeanAge, and AgeBracket.
auto_annotation_mapping = cluster_id_df.set_index('ClusterID (PoolClean)')['AutoAnnotation']
mean_age_mapping = cluster_id_df.set_index('ClusterID (PoolClean)')['MeanAge']
age_bracket_mapping = cluster_id_df.set_index('ClusterID (PoolClean)')['AgeBracket']

# Map the values to adata_ref.obs.
adata.obs['AutoAnnotation'] = adata.obs['cluster_id'].map(auto_annotation_mapping)
adata.obs['MeanAge'] = adata.obs['cluster_id'].map(mean_age_mapping)
adata.obs['AgeBracket'] = adata.obs['cluster_id'].map(age_bracket_mapping)

# Ensure the data types are consistent.
adata.obs['AutoAnnotation'] = adata.obs['AutoAnnotation'].astype(str)
adata.obs['cluster_id'] = adata.obs['cluster_id'].astype(str)
adata.obs['MeanAge'] = adata.obs['MeanAge'].astype(int)
adata.obs['AgeBracket'] = adata.obs['AgeBracket'].astype(str)

# Check for duplicates and make var_names unique if necessary
if not adata.var_names.is_unique:
    adata.var.index = pd.Index(adata.var.index).astype(str)  # Ensure var.index is a standard Index
    adata.var_names_make_unique()

#Preprocess reference dataset without subsetting to HVGs.
print('Preprocessing dataset...')
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat', batch_key='sample_id', subset=False)

#PCA, KNN, UMAP reference dataset.
print('Running PCA, KNN, and UMAP...')
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata)

#Save reference
print('Saving dataset...')
anndata_invivo_human_ref = directory_path + '/invivo_human_ref.h5ad'
adata.write(anndata_invivo_human_ref)
