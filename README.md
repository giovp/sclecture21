# sclecture21
Project repo for sc analysis course.

## Useful info for LRZ
login and pass were shared. When logged in:
```bash
module load python
conda create -n py38 python=3.8
source activate py38
# install your software, e.g. 
conda install numpy pandas matplotlib seaborn
#then clone and install scanpy
git clone https://github.com/theislab/scanpy.git
cd scanpy
pip install -e .
# and install useful scanpy modules
pip install leidenalg pynndescent
# and maybe jupyter lab/nb ?
pip install jupyterlab
```
To start jupyter
```bash
# make new login setting the port
ssh -Y lxlogin1.lrz.de -l di82cof -L 6006:localhost:6006
# then load usual stuff
module load python
source activate py38
jupyter-lab --no-browser --port=6006
# copy URL in your browser
```
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
