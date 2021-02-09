# Image Analysis and Segmentation

This branch contains the workflow and background information to the project 'Image Analysis'.

## 1 Workflow

The python file 'analysis_and_visualization.py' in the python_files subfolder displays all steps of the workflow.
For each step and its description in the workflow, the associated function is stated.

### 1.1 The Data Set
 
The used data set 'Human Breast Cancer: Whole Transcriptome Analysis' has been taken from the 10x genomics website (https://support.10xgenomics.com/spatial-gene-expression/datasets). 
It comes with an image alongside gene counts.

### 1.2 Pre-Processing

#### 1.2.1 Basic Filtering

By calculating basic quality contorl metrics, the total number of cell counts, number of expressed genes as well as the fraction of mitochondrial genes can be determined.
Based on manualy threshold decisions, the data set total cell and not often detected genes can be reduced by omitting the respective outlier regions.
[def pre_processing()]

#### 1.2.2 Normalizaion

Normalize the remaining data set after log transformation. Add the highly variable genes to the data set.
[def pre_processing()]

### 1.3 Clustering

#### 1.3.1 Manifold embedding and clustering

Cluster the manifold encoded by transcriptional similarity (== on gene expression scale) to see the underlying structure given by counts and genes.
[def pre_processing()]

#### 1.3.2 Visualization in spatial coordinates

Visualize how counts and genes by counts behave spatially (== in spatial dimension, in the cell, how the tissue is organized).
Again, the manifold can be clustered to detect the underlying structure, but spatially in the framework of the tissue.
[def pre_processing()]

#### 1.3.3 Marker genes of specific clusters

By taking a specific cluster, one can visualize the expression levels, e.g. of the highest expressed genes.
Furthermore, a gene with optimal expression in the specific cluster can then recapitulate the spatial structure of the cluster. 
[def pre_processing()] & [def statistical_ranking_tests()]

### 1.4 Image Resolution

Directly investigate the image by correctly scaling the spatial coordinates to the image pixels.
[def zoom_into_tissue()]

### 1.5 Image feature space vs Gene expression space

Compare spots with features of the tissue in the feature space with the gene expression space with respect to their clusters.
Visualize in an heat map what percentage of labeled spots from gene expression space clusters can be found in the respective image feature space clusters.
[def image_feature_space()] & [def if_clusters_in_ge_space()] & [def heatmap_percentages_ge_in_if()]

### 1.6 Marker Gene Identification

#### 1.6.1 Highly expressed marker genes found in the tissue

Discover the top expressed genes in the clusters of the tissue by using a statistical t-test or wilcoxon test to rank all genes.
The result is heavily depending on the choice of the statistical ranking test.
[def pre_processing()] & [def statistical_ranking_tests()]

#### 1.6.2 Find marker genes that are highly correlated to cancer in the tissue

Read papers on the current research to learn about those genes who are highly correlated to cancer.
Then, look in the tissue for those genes and visualize their expression level.
[def visualize_specific_genes()]


## 2 Code

All results including images which have been created during the project can be recreated by running the functions in the python file 'analysis_and_visualization.py'.
The appropriate function for each step of the workflow is added in the specific steps at the end.

