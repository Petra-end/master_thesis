print('Importing packages')
# imports
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os

import gc
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import scvi
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from scipy import sparse

# this line forces theano to use the GPU and should go before importing cell2location
#os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'
# if using the CPU uncomment this:
#os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32,openmp=True,force_device=True'

import cell2location
from cell2location.models import RegressionModel
from cell2location.utils import select_slide
from cell2location.utils.filtering import filter_genes
from cell2location.plt import plot_spatial

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

#os.environ["CUDA_VISIBLE_DEVICES"] = "0"

path = '/beegfs/home/pmatyskova/singularity/'
path_binomial = path + 'results/binomialf2'

adata = sc.read_h5ad("/beegfs/scratch/bruening_scratch/lsteuernagel/data/whole_brain_atlas/Macosko/spatial/macosko_hypo_pucks.h5ad")
slices = np.unique(adata.obs['PuckID'])

print('Opening Hypomap')
hm_sc = sc.read_h5ad(path + 'hypoMap.h5ad')
hm_sc.X = hm_sc.raw.X
del hm_sc.raw

#FILTERING
print('Filtering Hypomap')
genes_to_exclude = pd.read_csv('/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/volumetric_analysis/genes_to_exclude.csv',
                                header = None)
hm_sc = hm_sc[:, ~hm_sc.var['features'].isin(genes_to_exclude.iloc[:,0])].copy()

selected = filter_genes(hm_sc, cell_count_cutoff=1, cell_percentage_cutoff2=0.001, nonz_mean_cutoff=1.01)
hm_sc = hm_sc[:, selected].copy()

#NB
print('Binomial')
cell2location.models.RegressionModel.setup_anndata(adata=hm_sc,
                                                    batch_key='Batch_ID',
                                                    labels_key='C185_named')
mod = RegressionModel(hm_sc)
mod.train(max_epochs=250, use_gpu=True) #250

print('Exporting binomial')
hm_sc = mod.export_posterior(hm_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True})


#save data
mod.save(f"{path_binomial}/", overwrite=True)
adata_file = f"{path_binomial}/data_nb.h5ad"
hm_sc.write(adata_file)
mod.plot_QC()
plt.savefig(f"{path_binomial}/nb_taining.png",
             bbox_inches = 'tight')
plt.close()

#NB
if 'means_per_cluster_mu_fg' in hm_sc.varm.keys():
    inf_aver = hm_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in hm_sc.uns['mod']['factor_names']]].copy()
else:
    inf_aver = hm_sc.var[[f'means_per_cluster_mu_fg_{i}' for i in hm_sc.uns['mod']['factor_names']]].copy()
inf_aver.columns = hm_sc.uns['mod']['factor_names']

for i in slices:
    
    path_cell2loc = path + 'results/cell2loc2_{0}'.format(i)
    os.mkdir(path_cell2loc)

    print('Opening data')
    #open data
    adata = sc.read_h5ad("/beegfs/scratch/bruening_scratch/lsteuernagel/data/whole_brain_atlas/Macosko/spatial/macosko_hypo_pucks.h5ad")
    adata = adata[adata.obs['PuckID'] == i] #select just one slice
    adata = adata[(adata.obs['nCount_Spatial'] > 300) & (adata.obs['nCount_Spatial'] < 5000)]
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    del adata.raw

    print('Check if break loop')
    if adata.obs.shape[0] == 0:
        continue

    print('Finding shared genes')
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata.var_names, inf_aver.index)
    adata = adata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    print('Cell2location')
    #Cell2loc
    cell2location.models.Cell2location.setup_anndata(adata=adata)
    mod = cell2location.models.Cell2location(adata, cell_state_df=inf_aver,
                                             N_cells_per_location=2, detection_alpha=200)
    mod.train(max_epochs=10000, train_size=1, use_gpu=True) #max_epochs=30000

    mod.plot_history(1000)
    plt.legend(labels=['full data training'])
    plt.savefig(f"{path_cell2loc}/c2l_loss.png",
                bbox_inches = 'tight')
    plt.close()

    print('Exporting cell2loc')
    adata = mod.export_posterior(adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True})


    #save
    mod.save(f"{path_cell2loc}/", overwrite=True)
    adata_file = f"{path_cell2loc}/sp.h5ad"
    adata.write(adata_file)

    mod.plot_QC()
    plt.savefig(f"{path_cell2loc}/c2l_taining.png",
                bbox_inches = 'tight')
    plt.close()
