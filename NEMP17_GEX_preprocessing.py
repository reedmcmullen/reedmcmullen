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
save_name = "_NEMP17"

#Read in the CSV file defining the sample names and the paths to the filtered feature barcode matrices output from 10X Genomics cellranger count pipeline.
matrices = pd.read_csv("/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17/filtered_feature_bc_matrix_paths.csv", index_col="sample")

#Read in the 10X genomics gene expression data as Anndata objects and store in a dictionary with sample_name, sample_adata key, value pairs
#Add the GEMwell metadata to the 'batch' column and make unique indices by appending '-batch' to them.
adata_dict = {}
for idx, (sample_name, row) in enumerate(matrices.iterrows()):
    sample_path = row['path']
    print(f'Reading in data for sample: {sample_name}')
    adata = sc.read_10x_mtx(sample_path, cache=True)
    # Modify the index of adata.obs to append the GEMwell number.
    adata.obs.index = adata.obs.index.str.replace('-1', '', regex=False)
    adata.obs.index = adata.obs.index + f'-{idx+1}'  # Append '-integer' to each index
    adata.obs['GEMwell'] = idx+1
    # Store the modified adata in the dictionary
    adata_dict[sample_name] = adata

#Concatenate AnnData objects, taking the union of variables (i.e. 'outer' join).
print('Concatenating AnnData objects')
adata = ad.concat(list(adata_dict.values()), join='outer')

#Save initial concatenated AnnData object as a H5AD file.
results_file_initial = directory_path + '/NEMP17_initial.h5ad'
adata.write(results_file_initial)

#Identify highly expressed genes.
sc.pl.highest_expr_genes(adata, n_top=20)

#Filter out cells based on a very conservative minimum number of genes and cells.
sc.pp.filter_cells(adata, min_genes=3)
sc.pp.filter_genes(adata, min_cells=3)

#Annotate mitochondrial genes and ribosomal genes and calculate QC metrics for each.
adata.var['mito'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mito'
adata.var['ribo'] = adata.var_names.str.startswith('RPS' or 'RPL') # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mito', 'ribo'], percent_top=None, log1p=False, inplace=True)

#Define QC cutoffs.
n_genes_cutoff = 8000
total_counts_cutoff = 25000
mito_cutoff = 10
ribo_cutoff = 5

#Plot violin plots of QC metrics with cutoffs.
fig, axes = plt.subplots(1, 4, figsize=(16, 6))
sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], stripplot=False, show=False)
axes[0].axhline(n_genes_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'total_counts', ax=axes[1], stripplot=False, show=False)
axes[1].axhline(total_counts_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'pct_counts_mito', ax=axes[2], stripplot=False, show=False)
axes[2].axhline(mito_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'pct_counts_ribo', ax=axes[3], stripplot=False, show=False)
axes[3].axhline(ribo_cutoff, color='red', linestyle='--')
#Save the figure
plt.tight_layout()
plt.savefig(save_name+'_qc_metrics_cutoffs.png')

#Filter outlier cells from each AnnData object based on QC metric plots.
adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff, :]
adata = adata[adata.obs.total_counts < total_counts_cutoff, :]
adata = adata[adata.obs.pct_counts_mito < mito_cutoff, :]
adata = adata[adata.obs.pct_counts_ribo < ribo_cutoff, :]

#Save the concatenated AnnData object with the raw count data and add the raw counts to the layer 'counts'.
results_file_counts = directory_path + '/NEMP17_counts.h5ad' # the file that will store the analysis results
adata.write(results_file_counts)
adata.layers['counts'] = adata.X.copy()

#Normalize to median total counts and logarithmize.
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

#Identify and plot highly variable genes from each AnnData object as an H5AD file.
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch", flavor='seurat')
sc.pl.highly_variable_genes(adata, save=save_name+'.png')

#Perform PCA dimensional reduction.
sc.tl.pca(adata, svd_solver='arpack')

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata)

#Compute the umap embedding based on knn-computed neighbors.
sc.tl.umap(adata)

# Cluster umap embeddings using leiden
res = 1
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=res)

#Save the AnnData object as an H5AD file.
results_file_preprocessed = directory_path + '/NEMP17_preprocessed.h5ad'
adata.write(results_file_preprocessed)

#Subset the AnnData object to 10% of cells for faster DEG testing and visualization.
# Set a seed for reproducibility.
np.random.seed(42)
# Calculate 10% of the total number of cells
n_cells = int(adata.n_obs * 0.1)
# Randomly sample indices without replacement
random_indices = np.random.choice(adata.obs.index, size=n_cells, replace=False)
# Subset the AnnData object
adata_subset = adata[random_indices].copy()
adata_subset

#Save the subset AnnData object.
results_file_preprocessed_subset = directory_path + '/NEMP17_preprocessed_subset.h5ad'
adata_subset.write(results_file_preprocessed_subset)
















