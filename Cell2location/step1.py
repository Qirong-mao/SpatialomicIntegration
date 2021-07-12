import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
data_type = 'float32'
import cell2location
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 3 

RNA_data = sc.read_csv("/tudelft.net/staff-bulk/ewi/insy/DBL/qirongmao/raw_data/scRNA/matrix.csv")
sc_data_folder = '/tudelft.net/staff-bulk/ewi/insy/DBL/qirongmao/processed_data/'
results_folder = '/tudelft.net/staff-bulk/ewi/insy/DBL/qirongmao/results/'

meta = pd.read_csv("/tudelft.net/staff-bulk/ewi/insy/DBL/qirongmao/raw_data/scRNA/metadata.csv")
RNA_data.obs=meta
RNA_data = RNA_data[~RNA_data.obs['cluster_label'].isna(), :]
sc.pp.filter_cells(RNA_data, min_genes=400)
sc.pp.filter_genes(RNA_data, min_cells=3)
RNA_data.raw=RNA_data
RNA_data.X=RNA_data.raw.X.copy()
sc.pp.normalize_total(RNA_data, target_sum=1e6)
sc.pp.log1p(RNA_data)
sc.tl.pca(RNA_data, svd_solver='arpack', n_comps=80, use_highly_variable=False)
RNA_data.obsm['X_pca'] = RNA_data.obsm['X_pca'][:, 1:]
RNA_data.varm['PCs'] = RNA_data.varm['PCs'][:, 1:]
sc.pp.neighbors(RNA_data, n_neighbors=6, n_pcs=79)
sc.tl.umap(RNA_data, min_dist = 0.8, spread = 1.5)
sc.tl.leiden(RNA_data, resolution = 1.5, key_added = "leiden")
RNA_data.var['SYMBOL']=RNA_data.var._stat_axis.values.tolist()
sc.pp.highly_variable_genes(RNA_data,n_top_genes=750)
is_var_gene = RNA_data.var['highly_variable']
RNA_data=RNA_data[:,RNA_data.var.highly_variable]

from cell2location import run_regression


r, RNA_data = run_regression(RNA_data, # input data object]

                   verbose=True, return_all=True,

                   train_args={
                    'covariate_col_names': ['leiden'], # column listing cell type annotation

                    # column listing technology, e.g. 3' vs 5',
                    # when integrating multiple single cell technologies corresponding
                    # model is automatically selected
                    'tech_name_col': None,

                    'n_epochs': 100, 'minibatch_size': 512, 'learning_rate': 0.01,

                    'use_cuda': True, # use GPU?

                    'train_proportion': 0.9, # proportion of cells in the training set (for cross-validation)
                    'l2_weight': True,  # uses defaults for the model

                    'readable_var_name_col': 'SYMBOL', 'use_raw': True},

                   model_kwargs={}, # keep defaults
                   posterior_args={}, # keep defaults

                   export_args={'path': results_folder + 'regression_model/', # where to save results
                                'save_model': True, # save pytorch model?
                                'run_name_suffix': ''})

reg_mod = r['mod']
