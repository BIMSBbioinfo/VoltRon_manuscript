# run with 
# LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6" python main.py
import CAST
import scanpy as sc
import torch, os
import numpy as np
from os.path import join as pj
import warnings
warnings.filterwarnings("ignore")

from CAST.utils import detect_highly_variable_genes
from CAST.utils import extract_coords_exp
from CAST.models.model_GCNII import Args
from CAST import CAST_MARK
from CAST.visualize import kmeans_plot_multiple
from CAST.CAST_Stack import reg_params
from CAST import CAST_STACK

output_path = 'results/' #### input the demo path
os.makedirs(output_path,exist_ok=True)

import time
import tracemalloc
import math

# multithread
import os
n = int(os.environ.get("SLURM_CPUS_PER_TASK") or len(os.sched_getaffinity(0)))
for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[var] = str(n)
torch.set_num_threads(n)
print("torch threads:", torch.get_num_threads())

# get slices
sliceA = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151673/151673_preprocessed.h5")
xy = np.asarray(sliceA.obsm["spatial"])  # works whether it's stored as an array or a DataFrame
sliceA.obs["column"] = xy[:, 0]  # x: horizontal pixel position
sliceA.obs["row"] = xy[:, 1]     # y: vertical pixel position
sliceB = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151674/151674_preprocessed.h5")
xy = np.asarray(sliceB.obsm["spatial"])  # works whether it's stored as an array or a DataFrame
sliceB.obs["column"] = xy[:, 0]  # x: horizontal pixel position
sliceB.obs["row"] = xy[:, 1]     # y: vertical pixel position

# align slices
print("registering images")
tracemalloc.start()
start = time.perf_counter()

def make_unique(s):
    if s in counts:
        counts[s] += 1
    else:
        counts[s] = 0
    return s if counts[s] == 0 else f"{s}_{counts[s]}"

counts = {}
sliceA.var.index = sliceA.var.index.map(make_unique)
counts = {}
sliceB.var.index = sliceB.var.index.map(make_unique)

# combine the two datasets
sample_list= ['Visium1','Visium2'] # [Query, Reference]
sdata = sliceA.concatenate(sliceB)

# rename the dataset labels to Visium1 and Visium2
batch_key = 'batch'
batch_rename = {'0' : sample_list[0],'1' : sample_list[1]}
sdata.obs.replace({batch_key:batch_rename},inplace=True)

# Filter for highly variable genes
sdata.var['highly_variable'] = detect_highly_variable_genes(sdata,batch_key=batch_key,n_top_genes=4000,count_layer='.X')
sdata = sdata[:,sdata.var['highly_variable']]

# Output the combined dataset
sdata.write_h5ad(f'{output_path}/Visium_vs_Visium.h5ad')

coords_raw,exps = extract_coords_exp(sdata, batch_key = 'batch', cols = ['column', 'row'], count_layer = '.X', data_format = 'norm1e4')
torch.save(coords_raw, f'{output_path}/coords_raw.pt')

# set the parameters for CAST Mark
args = Args(
    dataname='task1', # name of the dataset, used to save the log file
    gpu = 0, # gpu id, set to zero for single-GPU nodes
    epochs=400, # number of epochs for training
    lr1= 1e-3, # learning rate
    wd1= 0, # weight decay
    lambd= 1e-3, # lambda in the loss function, refer to online methods
    n_layers=2, # number of GCNII layers, more layers mean a deeper model, larger reception field, at a cost of VRAM usage and computation time
    der=0.5, # edge dropout rate in CCA-SSG
    dfr=0.3, # feature dropout rate in CCA-SSG
    use_encoder=True, # perform a single-layer dimension reduction before the GNNs, helps save VRAM and computation time if the gene panel is large
    encoder_dim=512, # encoder dimension, ignore if `use_encoder` set to `False`
)

# run CAST Mark
embed_dict = CAST_MARK(coords_raw,exps,output_path,args = args,graph_strategy='delaunay')

# set the parameters for CAST Stack
query_sample = sample_list[0]
params_dist = reg_params(dataname = query_sample,
                            gpu = 0 if torch.cuda.is_available() else -1,
                            #### Affine parameters
                            iterations=150,
                            dist_penalty1=0,
                            bleeding=500,
                            d_list = [3,2,1,1/2,1/3],
                            attention_params = [None,3,1,0],
                            #### FFD parameters
                            dist_penalty2 = [0],
                            alpha_basis_bs = [0],
                            meshsize = [8],
                            iterations_bs = [1],
                            attention_params_bs = [[None,3,1,0]],
                            mesh_weight = [None])
# set the alpha basis for the affine transformation
params_dist.alpha_basis = torch.Tensor([1/1000,1/1000,1/50,5,5]).reshape(5,1).to(params_dist.device)

# run CAST Stack
coord_final = CAST_STACK(coords_raw,embed_dict,output_path,sample_list,params_dist,sub_node_idxs = None, rescale = True)

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")
