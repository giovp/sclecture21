# sclecture21
Project repo for sc analysis course.

## Available projects (to be discussed)
### Spatial statistics and global tissue context  
A central goal of spatial transcriptomics is to understand global spatial organization and cellular neighborhood structure. Spatial statistics can give insights on such properties. Tasks include:
- Select 2-3 spatial datasets
- Clustering and diff exp analysis
- Ripley's K stats on clusters
- Neighborhood enrichment analysis
- Moran's I and [Spark](https://github.com/xzhoulab/SPARK) spatially variable genes (SVG) selection
- Compare SVG to DiffExp result
- Interpret results with biological insights

### Image analysis and segmentation
10x genomics Visium data comes weith an image alongside gene counts. An important question that can be asked is how morphology relates to gene expression and tissue organization. Tasks include:
- Select 2-3 spatial datasets (both HnE and Fluorescence)
- Clustering and diff exp analysis
- Compute image features and and feature importance wrt cluster assignment
    - Assess feature importance with different ML/DL models
    - Build simple self-supervised model and compare embeddings with standard image features
- Evaluate segmentation strategies and methods in HnE and Fluorescent, methods such
    - [Cellpose](https://github.com/MouseLand/cellpose)
    - [Stardist](https://github.com/mpicbg-csbd/stardist)
- Compare segmentation results to image features and clusters
