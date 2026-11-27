# run with 
# LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6" python main.py
import scanpy as sc
import torch
import matplotlib.pyplot as plt
import warnings
import sys
import anndata as ad
from pathlib import Path
from GEASO.utils.utils import intersect, feature_extract
from GEASO.alignment.coarse_alignment import coarse_alignment
from GEASO.alignment.elastic_registration import ElasticRegistration

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
slice1 = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151673/151673_preprocessed.h5")
slice2 = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151674/151674_preprocessed.h5")

# align slices
print("registering images")
tracemalloc.start()
start = time.perf_counter()

# preprocess
sc.pp.filter_cells(slice1, min_genes=10)
sc.pp.filter_genes(slice1, min_cells=3)
slice1.layers['counts'] = slice1.X.copy()
sc.pp.normalize_total(slice1)
sc.pp.log1p(slice1)
sc.pp.highly_variable_genes(slice1, n_top_genes=2000)
sc.pp.filter_cells(slice2, min_genes=10)
sc.pp.filter_genes(slice2, min_cells=3)
slice2.layers['counts'] = slice2.X.copy()
sc.pp.normalize_total(slice2)
sc.pp.log1p(slice2)
sc.pp.highly_variable_genes(slice2, n_top_genes=2000)
slice1.var_names_make_unique()
slice2.var_names_make_unique()
section_ids = ['slice2', 'slice1']
slice2.obs['batch'] = section_ids[0]
slice1.obs['batch'] = section_ids[1]

# feature processing
common_genes = intersect(slice2.var_names, slice1.var_names)
slice2 = slice2[:, common_genes]
slice1 = slice1[:, common_genes]
adata_list = [slice2, slice1]
emb_dict, loss_log, model = feature_extract(adata_list)
torch.cuda.empty_cache()
slice2.obsm['latent'] = emb_dict[section_ids[0]].cpu().numpy()
slice1.obsm['latent'] = emb_dict[section_ids[1]].cpu().numpy()

# caorse reg
R, t, scale, P, spatial_scales, spatial_means, guided_source, guided_target, coarse_source, coarse_target = coarse_alignment(source=slice2, target=slice1, top_K=10, dis_metric='cos', use_latent=True, scale_c=False)
slice2.obsm['spatial_coarse'] = coarse_source.copy()
slice1.obsm['spatial_coarse'] = coarse_target.copy()
adata_cat = ad.concat([slice2, slice1], label='batch')

# fine reg
device = "cpu"
register_model = ElasticRegistration(
    init_P=P,
    guided_source=guided_source,
    guided_target=guided_target,
    source_slice=slice2,
    target_slice=slice1,
    normalize_scales=spatial_scales,
    normalize_means=spatial_means,
    spatial_key='spatial_coarse',
    device=device,
)

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")
