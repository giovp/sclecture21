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
Based on manualy threshold decisions, the data set cell and gene numbers can be reduced by omitting outlier regions.

