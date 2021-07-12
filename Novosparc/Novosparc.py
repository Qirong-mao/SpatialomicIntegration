import os
import numpy as np
import pandas as pd
import scanpy as sc
import novosparc
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist, squareform, pdist
from scipy.stats import ks_2samp

import os
import numpy as np
import pandas as pd
import scanpy as sc
import altair as alt

# The dataset of scRNA-seq dataset repository
data_dir = '/your_repository'
data_path = os.path.join(data_dir, 'sc.h5ad')
dataset1 = sc.read(data_path)
gene_names = dataset1.var.index.tolist()

ncells, ngenes = dataset1.shape

print('number of cells: %d' % ncells)
print('number of genes: %d' % ngenes)


# Preprocessing
spatial_data=sc.read_visium(path='/tudelft.net/staff-bulk/ewi/insy/DBL/qirongmao/raw_data/spatial/V1_Human_Brain_Section_2',count_file='filtered_feature_bc_matrix.h5',load_images=True)
sc.pp.normalize_per_cell(spatial_data)
sc.pp.log1p(spatial_data)
# Loading  prior location information
location_file=pd.read_csv('/tudelft.net/staff-bulk/ewi/insy/DBL/qirongmao/raw_data/spatial/V1_Human_Brain_Section_2/spatial/tissue_positions_list.csv',header=None)
location=spatial_data.obs[['array_col','array_row']]
locations_apriori=location.to_numpy()
# Deriving intersection genes
df_dge = spatial_data.to_df()
dataset_dge = dataset1.to_df()
marker_genes=df_dge.columns.values.tolist()
genes=dataset_dge.columns.values.tolist()
retB = list(set(marker_genes).intersection(set(genes)))

# Running Novosparc
tissue = novosparc.cm.Tissue(dataset=dataset1, locations=locations_apriori)
idx = np.arange(ngenes)
alpha_linear = 1.0
df_marker = df_dge # using expression from data as markers

# using a specific gene, e.g. PKM
marker_name = retB
marker_idx = idx[:597]
marker_names = [marker_name]

# compute linear cost only
tissue.setup_linear_cost(marker_idx, df_marker[marker_names[0]].values)

# reconstruct
tissue.reconstruct(alpha_linear=alpha_linear, verbose=False)

df_sdge = pd.DataFrame(tissue.sdge.T, columns=genenames)

sdge = tissue.sdge
dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=gene_names))
dataset_reconst.obsm['spatial'] = locations_apriori

# Generating Figure2, the format of output figure (how many figures in a row) needs parameter modification of _plotting.py in Novosparc package path.
gw = tissue.gw
ngw = (gw.T / gw.sum(1)).T
cell_idx=pd.read_csv("result.csv")
cell_cols=cell_idx.columns.values.tolist()
title=["L1 (Novosparc)","L1 (Cell2location)","L2 (Novosparc)","L2 (Cell2location)","L3 (Novosparc)","L3 (Cell2location)","L4 (Novosparc)","L4 (Cell2location)","L4ab (Novosparc)","L4ab (Cell2location)","L4c (Novosparc)","L4c (Cell2location)","L5 (Novosparc)","L5 (Cell2location)","L5a (Novosparc)","L5a (Cell2location)","L5b (Novosparc)	L5b (Cell2location)","L6 (Novosparc)","L6 (Cell2location)","L6a (Novosparc)","L6a (Cell2location)","L6b (Novosparc)","L6b (Cell2location)","WM (Novosparc)","WM (Cell2location)"]
dataset_reconst.obs = pd.DataFrame(cell_idx)
novosparc.pl.embedding(dataset_reconst, cell_cols)
from matplotlib import pyplot as plt
plt.savefig('Cellprob1.png')
