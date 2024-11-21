#! /usr/bin/env python3
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

#Define working directory and settings.
sc.set_figure_params(figsize=(8, 8), frameon=True)
torch.set_float32_matmul_precision("high")
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/refmap_NEMP17'
os.chdir(directory_path)
save_name='NEMP17'

#Load in formatted, subset, and preprocessed reference dataset.
print('Loading reference dataset...')
ref_preprocessed = directory_path + '/ref_preprocessed.h5ad'
adata_ref = sc.read_h5ad(ref_preprocessed)

#Load in formatted, subset, and preprocessed query dataset.
print('Loading query dataset...')
query_preprocessed = directory_path + '/query_preprocessed.h5ad'
adata_query = sc.read_h5ad(query_preprocessed)

#Train the reference scANVI model and save.
print('Training reference scANVI model...')
scvi_ref = scvi.model.SCVI.load(directory_path + '/ref_scvi_model/', adata_ref)
adata_ref.obs["labels_scanvi_AA"] = adata_ref.obs["AutoAnnotation"].values
scanvi_ref = scvi.model.SCANVI.from_scvi_model(scvi_ref, unlabeled_category="Unknown", labels_key="labels_scanvi_AA")
scanvi_ref.train(max_epochs=20, n_samples_per_label=100)
scanvi_ref.save(directory_path + '/ref_scanvi_model_AA/', overwrite=True)

#Run neighbor finding (using scANVI latent representation), leiden clustering, and umap.
print('Reference dataset neighbor finding, clustering, and UMAP with scANVI latent representation...')
adata_ref.obsm["X_scANVI_AA"] = scanvi_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scANVI_AA")
sc.tl.leiden(adata_ref, flavor='igraph', n_iterations=2)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref, color=['leiden', 'CellClass', 'Subregion', 'TopLevelCluster'], ncols=2, wspace=0.5, save=f'_{save_name}_ref_X_scANVI_AA.png')
ref_preprocessed = directory_path + '/ref_preprocessed.h5ad'
adata_ref.write(ref_preprocessed, compression='gzip')

#Set up and train the query scANVI model and save.
print('Updating and training reference scANVI model with query dataset...')
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, directory_path + '/ref_scanvi_model_AA/')
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10)
scanvi_query.save(directory_path + '/query_scanvi_model_AA/', overwrite=True)

#Add scANVI latent representation embeddings to anndata object.
print('Predicting annotations for query dataset...')
adata_query.obsm["predictions_scANVI_AA"] = scanvi_query.get_latent_representation()
adata_query.obs["predictions_scANVI_AA"] = scanvi_query.predict()
query_preprocessed = directory_path + '/query_preprocessed.h5ad'
adata_query.write(query_preprocessed, compression='gzip')