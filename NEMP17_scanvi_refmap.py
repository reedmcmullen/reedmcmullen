#Load required packages and modules.
import os
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
import random

#Define working directory and settings.
sc.set_figure_params(figsize=(8, 8), frameon=True)
torch.set_float32_matmul_precision("high")
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/refmap_NEMP17'
os.chdir(directory_path)
save_name='NEMP17'

#Load and format reference dataset.
adata_ref = sc.read_h5ad('/wynton/group/pollen/jding/Sara/linnarsson/human_dev_GRCh38-3.0.0.h5ad')
for col in ['CellClass', 'Region', 'Subdivision', 'Subregion', 'Tissue', 'donor_id', 'dissection', 'sample_id',]:
    adata_ref.obs[col] = adata_ref.obs[col].str[2:-1]
adata_ref.obs['CellID'] = adata_ref.obs.index
adata_ref.obs['CellID'] = adata_ref.obs['CellID'].str[2:-1]
adata_ref.obs.index = adata_ref.obs['CellID']
del adata_ref.obs['CellID']

#Subset reference dataset to ectodermal-lineage derived cell types and then to 100k cells.
adata_ref = adata_ref[~adata_ref.obs['CellClass'].isin(['Erythrocyte', 'Fibroblast', 'Immune', 'Vascular'])]
random.seed(0)
adata_ref=adata_ref[np.random.choice(adata_ref.obs.index,100000,replace=False),:]

#Format reference dataset.
adata_ref.raw = adata_ref
adata_ref.layers["counts"] = adata_ref.X.copy()
adata_ref.var.index = adata_ref.var['Gene']
del adata_ref.var['Gene']
adata_ref.X = adata_ref.raw.X
adata_ref.obs['cluster_id'] = adata_ref.obs['cluster_id'].astype(str)
adata_ref.obs['TopLevelCluster'] = adata_ref.obs['TopLevelCluster'].astype(str)

#Save initial reference dataset.
ref_initial = directory_path + '/ref_initial.h5ad'
adata_ref.write(ref_initial, compression='gzip')

#Preprocess reference dataset
if not adata_ref.var_names.is_unique:
    adata_ref.var.index = pd.Index(adata_ref.var.index).astype(str)  # Ensure var.index is a standard Index
    adata_ref.var_names_make_unique()
sc.pp.filter_cells(adata_ref, min_genes=100)
sc.pp.filter_genes(adata_ref, min_cells=10)
sc.pp.normalize_total(adata_ref, exclude_highly_expressed=True)
sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(adata_ref, flavor='seurat', batch_key='sample_id', subset=False)
sc.tl.pca(adata_ref, svd_solver='arpack')
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)

#Save the preprocessed reference dataset.
ref_preprocessed = directory_path + '/ref_preprocessed.h5ad'
adata_ref.write(ref_preprocessed, compression='gzip')

#Load in query dataset.
adata_query = sc.read_h5ad('/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17/NEMP17_preprocessed_subset_harmony.h5ad')

#Preprocess query dataset.
if not adata_query.var_names.is_unique:
    adata_query.var.index = pd.Index(adata_query.var.index).astype(str)  # Ensure var.index is a standard Index
    adata_query.var_names_make_unique()
sc.pp.highly_variable_genes(adata_query, flavor='seurat', batch_key='GEMwell', subset=False)

#Save preprocessed query dataset
query_preprocessed = directory_path + '/query_preprocessed.h5ad'
adata_query.write(query_preprocessed, compression='gzip')

#Plot both datasets.
sc.pl.umap(adata_query, color='leiden', save=f'_{save_name}_query_X_pca.png')
sc.pl.umap(adata_ref, color=['CellClass', 'Subregion', 'TopLevelCluster'], wspace=0.5, save=f'_{save_name}_ref_X_pca.png')

#Subset the reference and query datasets to the HVGs in the query dataset.
shared_var_names = adata_ref.var_names.intersection(adata_query.var_names[adata_query.var['highly_variable']])
adata_ref = adata_ref[:, shared_var_names].copy()
adata_query = adata_query[:, shared_var_names].copy()

#Save subset and preprocessed reference and query datasets.
ref_preprocessed = directory_path + '/ref_preprocessed.h5ad'
adata_ref.write(ref_preprocessed, compression='gzip')
query_preprocessed = directory_path + '/query_preprocessed.h5ad'
adata_query.write(query_preprocessed, compression='gzip')

#Train reference scVI model.
scvi.model.SCVI.setup_anndata(adata_ref, batch_key='sample_id', layer='counts')
scvi_ref = scvi.model.SCVI(adata_ref, use_layer_norm="both", use_batch_norm="none", encode_covariates=True, dropout_rate=0.2, n_layers=2)
scvi_ref.train()

#Save the reference scVI model.
scvi_ref.save(directory_path + '/ref_scvi_model/', overwrite=True)

#Run neighbor finding (using latent representation), leiden clustering, and umap.
adata_ref.obsm["X_scVI"] = scvi_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scVI")
sc.tl.leiden(adata_ref, flavor='igraph', n_iterations=2)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref, color=['CellClass', 'Region', 'Subregion', 'TopLevelCluster'], ncols=2, wspace=0.5, save=f'_{save_name}_ref_X_scVI.png')

#Save subset and preprocessed reference dataset with scVI latent representations.
ref_preprocessed = directory_path + '/ref_preprocessed.h5ad'
adata_ref.write(ref_preprocessed, compression='gzip')

#Update the reference scVI model with the query dataset
adata_query.obs['sample_id'] = adata_query.obs['GEMwell'].astype(str)
scvi.model.SCVI.prepare_query_anndata(adata_query, directory_path + '/ref_scvi_model/')
scvi_query = scvi.model.SCVI.load_query_data(adata_query, directory_path + '/ref_scvi_model/')

#Train the query data.
scvi_query.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

#Save the reference model
scvi_query.save(directory_path + '/query_scvi_model/', overwrite=True)

#Run neighbor finding (using latent representation), leiden clustering, and umap.
adata_query.obsm["X_scVI"] = scvi_query.get_latent_representation()
sc.pp.neighbors(adata_query, use_rep="X_scVI")
sc.tl.leiden(adata_query, flavor='igraph', n_iterations=2)
sc.tl.umap(adata_query)
sc.pl.umap(adata_query, color='leiden', save=f'_{save_name}_query_X_scVI.png')

#Save the subset and preprocessed query dataset with the scVI latent representation.
query_preprocessed = directory_path + '/query_preprocessed.h5ad'
adata_query.write(query_preprocessed, compression='gzip')

#Concatenate datasets
adata_query.obs['query_leiden'] = adata_query.obs['leiden']
adata_ref.obs['ref_leiden'] = adata_ref.obs['leiden']
adata_concat = anndata.concat([adata_ref, adata_query], join='outer', keys=['reference', 'query'], label='dataset')

#Run neighbor finding (using latent representation), leiden clustering, and umap.
sc.pp.neighbors(adata_concat, use_rep="X_scVI")
sc.tl.leiden(adata_concat, flavor='igraph', n_iterations=2)
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['leiden', 'query_leiden', 'ref_leiden', 'dataset', 'Region', 'Subregion', 'CellClass', 'TopLevelCluster'], ncols=4, wspace=1, save=f'_{save_name}_ref_query_concat_X_scVI.png')

#Load scvi model from reference dataset.
scvi_ref = scvi.model.SCVI.load(directory_path + '/ref_scvi_model/', adata_ref)

#Load in cluster AutoAnnotations.
cluster_id_df = pd.read_excel(directory_path + '/science.adf1226_table_s2.xlsx', index_col='PoolOrder')
cluster_id_df = cluster_id_df[cluster_id_df['ClusterID (PoolClean)']!= '--']
cluster_id_df['ClusterID (PoolClean)'] = cluster_id_df['ClusterID (PoolClean)'].astype(str)
mapping = cluster_id_df.set_index('ClusterID (PoolClean)')['AutoAnnotation']
adata_ref.obs['AutoAnnotation'] = adata_ref.obs['cluster_id'].map(mapping)
adata_ref.obs['AutoAnnotation'] = adata_ref.obs['AutoAnnotation'].astype(str)
adata_ref.obs['cluster_id'] = adata_ref.obs['cluster_id'].astype(str)

#Train reference scANVI model.
adata_ref.obs["labels_scanvi"] = adata_ref.obs["AutoAnnotation"].values
scanvi_ref = scvi.model.SCANVI.from_scvi_model(scvi_ref, unlabeled_category="Unknown", labels_key="labels_scanvi")
scanvi_ref.train(max_epochs=20, n_samples_per_label=100)

#Save the scanvi reference model
scanvi_ref.save(directory_path + '/ref_scanvi_model/', overwrite=True)

#Run neighbor finding (using latent representation), leiden clustering, and umap.
adata_ref.obsm["X_scANVI"] = scanvi_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
sc.tl.leiden(adata_ref, flavor='igraph', n_iterations=2)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref, color=['leiden', 'CellClass', 'Subregion', 'TopLevelCluster'], ncols=2, wspace=0.5, save=f'_{save_name}_ref_X_scANVI.png')

#Save subset and preprocessed reference dataset with scVI and scANVI latent representations.
ref_preprocessed = directory_path + '/ref_preprocessed.h5ad'
adata_ref.write(ref_preprocessed, compression='gzip')

#Set up and train the query scANVI model.
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, directory_path + '/ref_scanvi_model/')
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10)

#Save the scanvi query model.
scanvi_query.save(directory_path + '/query_scanvi_model/', overwrite=True)

#Add scANVI latent represetnations to anndata object.
adata_query.obsm["predictions_scANVI"] = scanvi_query.get_latent_representation()
adata_query.obs["predictions_scANVI"] = scanvi_query.predict()

#Save the subset and preprocessed query dataset with scVI and scANVI latent representations.
query_preprocessed = directory_path + '/query_preprocessed.h5ad'
adata_query.write(query_preprocessed, compression='gzip')

#Concatenate the refernce and query datasets and save.
adata_concat = anndata.concat([adata_ref, adata_query], join='outer', keys=['reference', 'query'], label='dataset')
adata_concat.uns['dataset_colors'] = ['#D3D3D3', '#000000']

#Run neighbor finding (using latent representation), leiden clustering, and umap.
sc.pp.neighbors(adata_concat, use_rep="predictions_scANVI")
sc.tl.leiden(adata_concat, flavor='igraph', n_iterations=2)
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['leiden', 'query_leiden', 'ref_leiden', 'dataset', 'Region', 'Subregion', 'CellClass', 'TopLevelCluster'], ncols=4, wspace=1, save=f'_{save_name}_ref_query_concat_X_scVI.png')

#Save the concatenated reference and query datasets.
ref_query_concat = directory_path + '/ref_query_concat_scanvi.h5ad'
adata_concat.write(ref_query_concat, compression='gzip')
