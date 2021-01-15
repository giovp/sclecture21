# Image Analysis and Segmentation

This branch contains the workflow and background information to the project 'Image Analysis'.

## 1. Workflow

The python file 'analysis_and_visualization.py' in the python_files subfolder displays all steps of the workflow.

### 1.1 The Data Set
 
The used data set 'Human Breast Cancer: Whole Transcriptome Analysis' has been taken from the 10x genomics website (https://support.10xgenomics.com/spatial-gene-expression/datasets). 
It comes with an image alongside gene counts.

### 1.2 Pre-Processing

#### 1.2.1 Basic Filtering

By calculating basic quality contorl metrics, the total number of cell counts, number of expressed genes as well as the fraction of mitochondrial genes can be determined.
Based on manualy threshold decisions, the data set total cell and not often detected genes can be reduced by omitting the respective outlier regions.

### 1.2.2 Normalizaion

Normalize the remaining data set after log transformation. Add the highly variable genes to the data set.

## 1.3 Clustering

### 1.3.1 Manifold embedding and clustering

Cluster the manifold encoded by transcriptional similarity (== on gene expression scale) to see the underlying structure given by counts and genes.

### 1.3.2 Visualization in spatial coordinates

Visualize how counts and genes by counts behave spatially (== in spatial dimension, in the cell, how the tissue is organized).
Again, the manifold can be clustered to detect the underlying structure, but spatially in the framework of the tissue.

### 1.3.3 Marker genes of specific clusters

By taking a specific cluster, one can visualize the expression levels, e.g. of the highest expressed genes.
Furthermore, a gene with optimal expression in the specific cluster can then recapitulate the spatial structure of the cluster. 

