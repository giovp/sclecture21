# file for analysing spatial transcriptomic data of human breast cancer, whole transcriptome analysis

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import imageio
import tqdm


'''
################# pre-processing ###############
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
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] > 40000], kde=False, bins=40, ax=axs[2])
#plt.show()

fig, axs = plt.subplots(1, 3, figsize=(15, 4))
sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[0])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 2500], kde=False, bins=40, ax=axs[1])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] > 7500], kde=False, bins=40, ax=axs[2])
#plt.show(fig)

# plot the fractional mitochondrial counts per total count
ax = plt.axes()
ax.scatter(adata.obs['total_counts'], adata.obs['total_counts_mt'])
ax.set_xlabel('Total number of counts')
ax.set_ylabel('Total number of mitochondrial counts')

# determine the (min,max)-thresholds by inspecting the filtering
cell_thresh_min = 100 #4000 #100
cell_thresh_max = 45000
gene_thresh_min = 2000 #1500 #2000
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
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)

# manifold embedding and clustering based on transcriptional similarity on different resolutions
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
for iRes in [1]:#0.25, 0.5, 0.75, 1]:
    sc.tl.leiden(adata, resolution=iRes, key_added=f'cluster_{iRes}')

    # plot some covariates to check for structure
    plt.rcParams['figure.figsize'] = (4,4)
    sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', f'cluster_{iRes}'], wspace=0.4)

# visualization in spatial coordinates
plt.rcParams['figure.figsize'] = (8,8)
sc.pl.spatial(adata, img_key='hires', color=['total_counts', 'n_genes_by_counts'])

# check for more structure by clustering with the above resolution, but now spatially
for iRes in [1]:#0.25, 0.5, 0.75, 1]:
    sc.pl.spatial(adata, img_key='hires', color=f'cluster_{iRes}', size=1.5)


    # change the region of interest and the transparency to get a deeper picture of the tissue
    sc.pl.spatial(adata, img_key='hires', color=f'cluster_{iRes}', size=1.5,
                  groups=['0', '5'], crop_coord=tuple([1200, 1700, 1900, 1000]), alpha=0.5)


# cluster marker genes by a t-test and plot via a heatmap
for iRes in [1]:#0.25, 0.5, 0.75, 1]:
    sc.tl.rank_genes_groups(adata, f'cluster_{iRes}', method='t-test')
    sc.pl.rank_genes_groups_heatmap(adata, groups='0', n_genes=10, groupby=f'cluster_{iRes}')

    # plot the specific gene
    sc.pl.spatial(adata, img_key='hires', color=[f'cluster_{iRes}', 'SCGB2A2'])

# image resolution reviisited --- actual image has been safed after import
spatial_data = adata.uns['spatial']['Parent_Visium_Human_BreastCancer']
spot_size = spatial_data['scalefactors']['spot_diameter_fullres']*0.5
img = tif
crop_coord = np.asarray([5000, 7500, 15000, 20000])
img_coord = (
    *crop_coord[:2],
    *np.ceil(img.shape[0] - crop_coord[2:4]).astype(int),
)
fig, ax = plt.subplots()
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

def circles(x, y, s, ax, marker=None, c="b", vmin=None, vmax=None, **kwargs):
    """"""
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.
    zipped = np.broadcast(x, y, s)
    patches = [Circle((x_, y_), s_) for x_, y_, s_ in zipped]
    collection = PatchCollection(patches, **kwargs)
    if isinstance(c, np.ndarray) and np.issubdtype(c.dtype, np.number):
        collection.set_array(c)
        collection.set_clim(vmin, vmax)
    else:
        collection.set_facecolor(c)
    ax.add_collection(collection)
    return collection

circles(xcoord, ycoord, s=spot_size, ax=ax)
plt.imshow(tif)
ax.set_xlim(img_coord[0], img_coord[1])
ax.set_ylim(img_coord[3], img_coord[2])
'''




################# image feature ###############
# read in the data set from 10x genomics as well as the large tif image
img = sq.im.ImageContainer('./data/Parent_Visium_Human_BreastCancer/image.tif')
adata = sc.datasets.visium_sge(sample_id='Parent_Visium_Human_BreastCancer')
adata.var_names_make_unique()

# define different feature calculation combinations
params = {
    # all features, corresponding only to tissue underneath spot
    'features_orig':
    {'features': 'summary', 'size': 1, 'scale': 1.0, 'mask_circle': True},
    # summary and histogram features with a bit more context, original resolution
    'features_context':
    {'features': 'summary', 'size': 2, 'scale': 1.0},
    # summary and histogram features with more context and at lower resolution
    'features_lowres' :
    {'features': 'summary', 'size': 4, 'scale': 0.25}
}

# extract features with the different parameters in a loop
for feature_name, cur_params in tqdm.tqdm(params.items()):
    # features will be saved in `adata.obsm[feature_name]`
    sq.im.calculate_image_features(adata, img, key_added=feature_name, n_jobs=4, **cur_params)

# add all data together
# fill nans
adata.obsm['features_orig'].fillna(value=0, inplace=True)
# combine features in one dataframe
adata.obsm['features'] = pd.concat([adata.obsm[f] for f in params.keys()], axis='columns')
# make sure that we have no duplicated feature names in the combined table
adata.obsm['features'].columns = ad.utils.make_index_unique(adata.obsm['features'].columns)


# helper function returning a clustering
def cluster_features(features, like=None):
    """Calculate leiden clustering of features.

    Specify filter of features using `like`.
    """
    # filter features
    if like is not None:
        features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    # interesting analysis: what scaling works best? use e.g. clusters in gexp as ground truth?
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    return adata.obs['leiden']


adata.obs['features_cluster'] = cluster_features(adata.obsm['features'])
