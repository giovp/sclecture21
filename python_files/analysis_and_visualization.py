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