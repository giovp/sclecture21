# file for analysing spatial transcriptomic data of human breast cancer, whole transcriptome analysis

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# check all versions included in scanpy
# set some figure markers
sc.logging.print_versions()
sc.set_figure_params(facecolor='white', figsize=(8,8))
sc.settings.verbosity = 3

# read in the dataset from 10xGenomics website
adata = sc.datasets.visium_sge(sample_id='Parent_Visium_Human_BreastCancer')
adata.var_names_make_unique()
print(f'Size of data set: n_obs x n_vars = {adata.shape[0]} x {adata.shape[1]}')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
print('Amount of mitochondrial cells: ' + str(len(['True' for i in adata.var['mt'] if i != False])))

# calculate standard Quality Control (QC) metrics and update the data set
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# perform some basic filtering for total cell counts and gene counts per cell
fig, axs = plt.subplots(2, 2, figsize=(15, 4))
sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0,0])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[0,1])
sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1,0])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[1,1])

# determine the (min,max)-thresholds by inspecting the filtering
cell_thresh_min = 4000
cell_thresh_max = 45000
gene_thresh_min = 2000
gene_thresh_max = 8000

# filter out the outlier cells based on the previously determined thresholds
sc.pp.filter_cells(adata, min_counts=cell_thresh_min)
sc.pp.filter_cells(adata, max_counts=cell_thresh_max)
adata = adata[adata.obs["pct_counts_mt"] < 20]
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)

# visualize taken outlier by red vertical line & plot distplot again (?)

# normalization of visium counts data to detect highly variable genes
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000
                            )
a = 4
