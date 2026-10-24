# run with 
# LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6" python main.py
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squidpy as sq
import pandas as pd
import numpy as np
import paste as pst
import ot
import time
import tracemalloc
import math

ot.backend.get_backend_list()

# get slices
sliceA = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151673/151673_preprocessed.h5")
sliceB = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151674/151674_preprocessed.h5")
sliceA.obsm['spatial'] = sliceA.obsm['spatial'].astype("int32")
sliceB.obsm['spatial'] = sliceB.obsm['spatial'].astype("int32")

# align slices
print("registering images")
tracemalloc.start()
start = time.perf_counter()

pi_AB = pst.pairwise_align(sliceA, sliceB)

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")