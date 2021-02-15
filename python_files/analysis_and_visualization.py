# file for analysing spatial transcriptomic data of human breast cancer, whole transcriptome analysis

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import imageio
import tqdm
import anndata as ad
import sys


################# pre-processing & basic visualizations ###############
def pre_processing():
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
        #sc.pl.spatial(adata, img_key='hires', color=f'cluster_{iRes}', size=1.5,
        #              groups=['0', '5'], crop_coord=tuple([1200, 1700, 1900, 1000]), alpha=0.5)

    return adata




################# gene marker visualization by different statistical ranking tests ###############
def statistical_ranking_tests(adata):
    # cluster marker genes by a t-test and plot via a heatmap
    for iRes in [1]:#0.25, 0.5, 0.75, 1]:
        sc.tl.rank_genes_groups(adata, f'cluster_{iRes}', method='t-test')
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby=f'cluster_{iRes}')
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, groupby=f'cluster_{iRes}')
        sc.pl.rank_genes_groups_matrixplot(adata, n_genes=10, groupby=f'cluster_{iRes}')

        # plot the specific gene
        #sc.pl.spatial(adata, img_key='hires', color=[f'cluster_{iRes}', 'SCGB2A2'])

    # cluster marker genes by a wilcoxin test and plot via a heatmap
    for iRes in [1]:  # 0.25, 0.5, 0.75, 1]:
        sc.tl.rank_genes_groups(adata, f'cluster_{iRes}', method='wilcoxon')
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby=f'cluster_{iRes}')
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, groupby=f'cluster_{iRes}')
        sc.pl.rank_genes_groups_matrixplot(adata, n_genes=10, groupby=f'cluster_{iRes}')

        # plot the specific gene
        #sc.pl.spatial(adata, img_key='hires', color=[f'cluster_{iRes}_wilcoxon', 'SCGB2A2'])

    return adata




################# zoom into tissue by performing a coordinate transformation ###############
def zoom_into_tissue(adata):
    # image resolution revisited --- actual image has been saved after import
    spatial_data = adata.uns['spatial']['Parent_Visium_Human_BreastCancer']
    xcoord = adata.obsm['spatial'][:, 0]
    ycoord = adata.obsm['spatial'][:, 1]
    spot_size = spatial_data['scalefactors']['spot_diameter_fullres']*0.5
    tif = imageio.imread('./python_files/image.tif')
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

    return adata




################# calculate the image feature clusters ###############
def image_feature_spcae():
    # read in the data set from 10x genomics as well as the large tif image
    img = sq.im.ImageContainer('./data/Parent_Visium_Human_BreastCancer/image.tif')
    adata_2 = sc.datasets.visium_sge(sample_id='Parent_Visium_Human_BreastCancer')
    adata_2.var_names_make_unique()

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
        # features will be saved in `adata_2.obsm[feature_name]`
        sq.im.calculate_image_features(adata_2, img, key_added=feature_name, n_jobs=4, **cur_params)

    # add all data together
    # fill nans
    adata_2.obsm['features_orig'].fillna(value=0, inplace=True)
    # combine features in one dataframe
    adata_2.obsm['features'] = pd.concat([adata_2.obsm[f] for f in params.keys()], axis='columns')
    # make sure that we have no duplicated feature names in the combined table
    adata_2.obsm['features'].columns = ad.utils.make_index_unique(adata_2.obsm['features'].columns)

    # get the umap + leiden analysis for the features
    features = adata_2.obsm['features']
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    # interesting analysis: what scaling works best? use e.g. clusters in gexp as ground truth?
    sc.pp.scale(adata)
    # calculate leiden clustering
    # compute principle component analysis coordinates, loadings and variance composition
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    # compute a neighborhood graph of observations
    sc.pp.neighbors(adata)
    # compute umap
    sc.tl.umap(adata)
    # cluster cells into subgroups using the Leiden algorithm and\
    # plot umap image feature space
    for iRes in [1]:#0.25, 0.5, 0.75, 1]:
        sc.tl.leiden(adata, resolution=iRes, key_added=f'cluster_{iRes}')

        # plot some covariates to check for structure
        plt.rcParams['figure.figsize'] = (4, 4)
        sc.pl.umap(adata, color=f'cluster_{iRes}')

    return adata




################# calculate the image feature clusters in gene expression space ###############
################# this will give a nice comparison to the heatmap 'heatmap_percentages_ge_in_if' ###############
def if_clusters_in_ge_space():
    # read in the data set from 10x genomics as well as the large tif image
    img = sq.im.ImageContainer('./data/Parent_Visium_Human_BreastCancer/image.tif')
    adata = sc.datasets.visium_sge(sample_id='Parent_Visium_Human_BreastCancer')
    adata.var_names_make_unique()

    ####### calculate the gene expression clusters
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    print('Amount of mitochondrial cells: ' + str(len(['True' for i in adata.var['mt'] if i != False])))

    # calculate standard Quality Control (QC) metrics and update the data set
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    cell_thresh_min = 100 #4000 #100
    cell_thresh_max = 45000

    # some filtering
    sc.pp.filter_cells(adata, min_counts=cell_thresh_min)
    sc.pp.filter_cells(adata, max_counts=cell_thresh_max)
    adata = adata[adata.obs["pct_counts_mt"] < 20]
    print(f"#cells after MT filter: {adata.n_obs}")
    sc.pp.filter_genes(adata, min_cells=10)

    # normalization of visium counts data to detect highly variable genes
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)

    # manifold embedding and clustering based on transcriptional similarity on different resolutions
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added=f'gene_expr_cluster')

    # plot some covariates to check for structure
    plt.rcParams['figure.figsize'] = (4, 4)
    sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'gene_expr_cluster'], wspace=0.4)

    ####### calculate the image feature clusters
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
        # compute principle component analysis coordinates, loadings and variance composition
        sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
        # compute a neighborhood graph of observations
        sc.pp.neighbors(adata)
        # cluster cells into subgroups  using the Leiden algorithm
        sc.tl.leiden(adata, key_added='leiden')

        return adata.obs['leiden']

    adata.obs['features_cluster'] = cluster_features(adata.obsm['features'])

    # plot umap image feature clusters in gene expression space
    plt.rcParams['figure.figsize'] = (4, 4)
    sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'features_cluster'])

    # plot the spatial cluster map
    plt.rcParams['figure.figsize'] = (8, 8)
    sc.pl.spatial(adata, img_key='hires', color='features_cluster', size=1.5)

    return adata




################# calculate the gene expression clusters in image feature space ###############
################# this will give a nice comparison to the heatmap 'heatmap_percentages_if_in_ge' ###############
def if_clusters_in_ge_space():
    # read in the data set from 10x genomics as well as the large tif image
    img = sq.im.ImageContainer('./data/Parent_Visium_Human_BreastCancer/image.tif')
    adata = sc.datasets.visium_sge(sample_id='Parent_Visium_Human_BreastCancer')
    adata.var_names_make_unique()

    ####### calculate the image feature clusters
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

    # create features object
    features = adata.obsm['features']
    # create temporary adata to calculate the clustering
    adata2 = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    # interesting analysis: what scaling works best? use e.g. clusters in gexp as ground truth?
    sc.pp.scale(adata2)
    # calculate leiden clustering
    # compute principle component analysis coordinates, loadings and variance composition
    sc.pp.pca(adata2, n_comps=min(10, features.shape[1] - 1))
    # compute a neighborhood graph of observations
    sc.pp.neighbors(adata2)
    # compute umap
    sc.tl.umap(adata2)
    # cluster cells into subgroups  using the Leiden algorithm
    sc.tl.leiden(adata2, key_added='features_cluster')

    # transfer the calculations back to adata
    adata.obs['features_cluster'] = adata2.obs['features_cluster']
    adata.obsm['X_pca'] = adata2.obs['X_pca']
    adata.obs['X_umap'] = adata2.obs['X_umap']

    # plot umap image feature clusters in gene expression space
    plt.rcParams['figure.figsize'] = (4, 4)
    sc.pl.umap(adata, color='features_cluster')

    ####### calculate the gene expression clusters
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    print('Amount of mitochondrial cells: ' + str(len(['True' for i in adata.var['mt'] if i != False])))

    # calculate standard Quality Control (QC) metrics and update the data set
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    cell_thresh_min = 100  # 4000 #100
    cell_thresh_max = 45000

    # some filtering
    sc.pp.filter_cells(adata, min_counts=cell_thresh_min)
    sc.pp.filter_cells(adata, max_counts=cell_thresh_max)
    adata = adata[adata.obs["pct_counts_mt"] < 20]
    print(f"#cells after MT filter: {adata.n_obs}")
    sc.pp.filter_genes(adata, min_cells=10)

    # normalization of visium counts data to detect highly variable genes
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)

    # manifold embedding and clustering based on transcriptional similarity on different resolutions
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, key_added=f'gene_expr_cluster')

    # plot some covariates to check for structure
    plt.rcParams['figure.figsize'] = (4, 4)
    sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'gene_expr_cluster'], wspace=0.4)

    return adata




################# plot heatmap with percentages of ge clusters in if clusters ###############
def heatmap_percentages_ge_if(adata, cluster_direction_factor):
    # cluster_direction_factor identifies the direction the heat map is plotted,
    # e.g. if or ge on x-axis
    # cluster.direction_factor = 'ge_if' xor cluster.direction_factor = 'if_ge' [type: str]

    df = adata.obs
    # get all labels for ge cluster & size of all individual clusters
    cluster_ge = []
    for iCluster in range(0, 10):
        cluster_save = []
        for iRow in range(0, len(df.index)):
            if df['gene_expr_cluster'][iRow] == str(iCluster):
                cluster_save.append(df.index[iRow])
        cluster_ge.append(cluster_save)

    for iCluster in range(0, len(cluster_ge)):
        print(len(cluster_ge[iCluster]))

    # get all labels for if cluster & size of all individual clusters
    cluster_if = []
    for iCluster in range(0, 16):
        cluster_save = []
        for iRow in range(0, len(df.index)):
            if df['features_cluster'][iRow] == str(iCluster):
                cluster_save.append(df.index[iRow])
        cluster_if.append(cluster_save)

    for iCluster in range(0, len(cluster_if)):
        print(len(cluster_if[iCluster]))

    df_2 = calculate_heatmap_percentages(cluster_ge, cluster_if, cluster_direction_factor)

    return df_2


def calculate_heatmap_percentages(cluster_ge, cluster_if, cluster_direction_factor):
    if cluster_direction_factor == 'ge_if':
        # compare those labels with other cluster
        percentages = []
        for iCluster_ge in range(0, len(cluster_ge)):
            clusters = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
            for iLabel_ge in range(0, len(cluster_ge[iCluster_ge])):
                for iCluster_if in range(0, len(cluster_if)):
                    if cluster_ge[iCluster_ge][iLabel_ge] in cluster_if[iCluster_if]:
                        clusters[iCluster_if].append(1)
            per = []
            for cluster in clusters:
                per.append(len(cluster) / len(cluster_ge[iCluster_ge]))
            percentages.append(per)

        # create data frame
        index_column = ['if_0', 'if_1', 'if_2', 'if_3', 'if_4', 'if_5', 'if_6', 'if_7', 'if_8', 'if_9', 'if_10', 'if_11',
                        'if_12', 'if_13', 'if_14', 'if_15']
        df_2 = pd.DataFrame(columns=['ge_0', 'ge_1', 'ge_2', 'ge_3', 'ge_4', 'ge_5', 'ge_6', 'ge_7', 'ge_8', 'ge_9'],
                            data=[])
        for iCol in range(0, np.shape(df_2)[1]):
            df_2[f'ge_{iCol}'] = pd.Series(percentages[iCol])
        df_2['new_index'] = index_column
        df_2 = df_2.set_index('new_index')

    elif cluster_direction_factor == 'if_ge':
        # compare those labels with other cluster
        percentages = []
        for iCluster_if in range(0, len(cluster_if)):
            clusters = [[], [], [], [], [], [], [], [], [], []]
            for iLabel_if in range(0, len(cluster_if[iCluster_if])):
                for iCluster_ge in range(0, len(cluster_ge)):
                    if cluster_if[iCluster_if][iLabel_if] in cluster_ge[iCluster_ge]:
                        clusters[iCluster_ge].append(1)
            per = []
            for cluster in clusters:
                per.append(len(cluster) / len(cluster_if[iCluster_if]))
            percentages.append(per)
    else:
        print('cluster_direction_factor does not have one of the correct input arguments!')
        sys.exit()

    # plot heatmap with values
    sns.heatmap(df_2, annot=True)

    return df_2


################# gene marker visualization by different statistical ranking tests ###############
def visualize_specific_genes():
    # look for specific genes which correlate heavily with cancer
    # decision based on papers 'identifying driver genes involving gene dysregulated expression' and
    # 'Identification of cancer driver genes based on nucleotide context'
    list_of_specific_cancer_related_genes = ['PTEN', 'PIK3CA', 'BRAF', 'NRAS', 'KRAS', 'EGFR',
                                             'TP53', 'KMT2D', 'KDM6A', 'ARIDA1', 'CREBBP',
                                             'GTF2I', 'ZFHX3', 'NFE2L2',
                                             'MAP3K1', 'GATA3', 'CDH1', 'ERBB2', 'PTEN', 'BRCA1',
                                             ]#'MET', 'LKB1', 'ALK', 'RET', 'ROS1','RAC1', 'CTNNB1', 'CDKN2A', 'IDH1', 'NOTCH1']
    # get data
    adata = pre_processing()

    # cluster marker genes by a t-test and wilcoxon test
    for Test in ['t-test', 'wilcoxon']:
        for iRes in [1]:  # 0.25, 0.5, 0.75, 1]:
            sc.tl.rank_genes_groups(adata, f'cluster_{iRes}', method=Test)

            # plot the specific gene, if it is even expressed
            for specific_gene in list_of_specific_cancer_related_genes:
                if specific_gene in adata.var_names:
                    sc.pl.spatial(adata, img_key='hires', color=[f'cluster_{iRes}', specific_gene])

    return adata


adata = if_clusters_in_ge_space()
