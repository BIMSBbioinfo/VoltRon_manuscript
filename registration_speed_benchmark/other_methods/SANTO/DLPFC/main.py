import easydict
from SANTO_utils import santo, simulate_stitching, evaluation
from SANTO_data import intersect
import numpy as np
from tqdm import tqdm
import scanpy as sc
import math
import torch
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import time
import tracemalloc
import math

# get slices
sliceA = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151673/151673_preprocessed.h5")
sliceB = sc.read_h5ad("../../../../../data/10X_Visium_DLPFC/151674/151674_preprocessed.h5")

# align slices
print("registering images")
tracemalloc.start()
start = time.perf_counter()

args = easydict.EasyDict({})
args.epochs = 20
args.lr = 0.001
args.k = 10
args.alpha = 0.1 # weight of transcriptional loss
args.diff_omics = False # whether to use different omics data
args.mode = 'align' # Choose the mode among 'align', 'stitch' and None
args.dimension = 2  # choose the dimension of coordinates (2 or 3)
args.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

align_source_cor,trans_dict = santo(sliceB, sliceA, args)

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")