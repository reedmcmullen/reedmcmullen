#! /usr/bin/env python3
#Import required packages and modules.
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from cnmf import cNMF

#Set settings and variable names.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/cnmf_NEMP17'
os.chdir(directory_path)
save_name= 'NEMP17'
output_directory = directory_path + '/cnmf'
run_name = 'NEMP17_all'
n_iter=200 #Number of NMF replicates.
num_highvar_genes=2000 #Number of over-dispersed genes to use.
K = ' '.join([str(i) for i in range(15,30)]) #Specify the Ks to use as a space separated list in this case "5 6 7 8 9 10".
seed = 14 #Specify a seed pseudorandom number generation for reproducibility
count_fn= anndata

#Load the dataset.
print('Loading the dataset...')
anndata = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17/NEMP17.h5ad'
adata = sc.read_h5ad(anndata)

#Run cNMF
print('Initializing and preparing the cNMF object...')
cnmf_obj = cNMF(output_dir=output_directory, name=run_name) #Initialize cNMF object.
cnmf_obj.prepare(counts_fn=count_fn, components=np.arange(15,30), n_iter=n_iter, seed=seed, num_highvar_genes=num_highvar_genes)
print('Running cNMF factorization...')
cnmf_obj.factorize(worker_i=0, total_workers=1)
print('Running cNMF combination')
cnmf_obj.combine()
cnmf_obj.k_selection_plot(close_fig=False)
selected_K = 22
density_threshold = 2.00
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
density_threshold = 0.3
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)










