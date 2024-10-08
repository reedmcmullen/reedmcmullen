#Import required packages and modules.
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os

#Define and set the current working directory.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17'
os.chdir(directory_path)

#Define string for use in naming saved figures.
save_name = "_NEMP17"

#Read in the csv defining the sample names and the paths to the filtered feature barcode matrices.
matrices = pd.read_csv("/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17/filtered_feature_bc_matrix_paths.csv", index_col="sample")

#Read in the 10X genomics gene expression data as Anndata objects and store in a dictionary with sample_name, sample_adata key, value pairs.
adata_dict = {}
for sample_name, sample_path in matrices.iterrows():
    print(f'Reading in data for sample: {str(sample)}')
    adata_dict[sample_name] = sc.read_10x_mtx(sample_path, cache=True)

# Concatenate AnnData objects, taking the union of variables (i.e. 'outer' join).
print('Concatenating AnnData objects')
adata = ad.concat(list(adata_dict.values()), join='outer')

# Save initial concatenated AnnData object as a H5AD file.
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

#Visualize QC metrics with violin plots.
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mito', 'pct_counts_ribo'], multi_panel=True, stripplot = False, save = save_name + "_qc_metrics.png")

#Define QC cutoffs.
n_genes_cutoff = 8000
total_counts_cutoff = 25000
mito_cutoff = 9
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
adata

#Save the concatenated AnnData object with the raw counts.
results_file_counts = directory_path + '/NEMP17_counts.h5ad' # the file that will store the analysis results
adata.write(results_file_counts)

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















