from pathlib import Path
import sys
import matplotlib.pyplot as plt
import cv2
import os
import anndata
import shutil
from PIL import Image
from PIL import ImageOps
from PIL import ImageEnhance
Image.MAX_IMAGE_PIXELS = 665836957
from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad

# time control
import time
import tracemalloc
import math

# multi-threading
import os
import torch
import SimpleITK as sitk
n = os.environ.get("SLURM_CPUS_PER_TASK") or str(len(os.sched_getaffinity(0)))
for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[var] = n
torch.set_num_threads(int(n))
cv2.setNumThreads(int(n))
sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(int(n))
print("Number of threads: " + n)

START_DIR = "../SPOmiAlign"
PROJECT_ROOT = START_DIR
spomialign_path = PROJECT_ROOT + "/SPOmiAlign"
if str(spomialign_path) not in sys.path:
    sys.path.insert(0, str(spomialign_path))

from tutorial_utils import (
    generate_ssi_visualization,
    get_tutorial_paths,
    run_image_to_image_alignment,
    run_reassignment,
    read_bgr,
)

# DATA_DIR, OUTPUT_ROOT = get_tutorial_paths(PROJECT_ROOT)
DATA_DIR = Path(".")
sample = "lungTMA1"
OUTPUT_ROOT = Path("./results/" + sample)
if not os.path.exists(OUTPUT_ROOT):
    os.mkdir(OUTPUT_ROOT)
else:
    shutil.rmtree(OUTPUT_ROOT)
    os.mkdir(OUTPUT_ROOT)
target_image_path= '../../../../../data/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif'
source_image_path= '../../../../../data/Xenium_R1/image_3_rotated.tif'
valid_source_landmarks = "../../../../../data/Xenium_R1/cells.csv.gz"

# rotate landmarks
_EXACT = {
    0: None,
    90: Image.Transpose.ROTATE_270,
    180: Image.Transpose.ROTATE_180,
    270: Image.Transpose.ROTATE_90,
}

def rotate_image(image, deg, fill=0):
    key = int(round(deg)) % 360
    if key in _EXACT and abs(deg - round(deg)) < 1e-9:
        op = _EXACT[key]
        return image.copy() if op is None else image.transpose(op)
    return image.rotate(-deg, expand=True, resample=Image.BICUBIC,
                        fillcolor=fill)

def _rotation_matrix(deg):
    radii = np.deg2rad(360 - deg)
    s, c = np.sin(radii), np.cos(radii)
    return np.array([[c, -s], [s, c]])

def rotate_keypoints(xy, deg, image=None):
    points = xy
 
    limits = np.asarray(image.size, dtype=float)          # (width, height)
    origin = limits / 2
 
    image = rotate_image(image, deg)
 
    rotated_limits = np.asarray(image.size, dtype=float)
    rotated_origin = rotated_limits / 2
 
    rotation_mat = _rotation_matrix(deg)

    p = points - origin
    p = p / limits                      # to fractional coordinates
    p = p @ rotation_mat.T
    p = p * rotated_limits              # back to pixels
    p = p + rotated_origin
 
    out = pd.DataFrame({"x": p[:, 0], "y": p[:, 1]})
    out.attrs["image"] = image
    return out

for p in [valid_source_landmarks]:
    p = Path(p)
    if not p.exists():
        raise FileNotFoundError(f"Landmark file not found: {p}")

def landmark_csv_to_adata(df, spatial_cols=("x", "y")):

    # try common column names for coordinates
    coord_options = [
        ("x", "y"),
        ("X", "Y"),
        ("x_coord", "y_coord"),
        ("xcoord", "ycoord"),
        ("x_px", "y_px"),
    ]
    for cols in coord_options:
        if all(c in df.columns for c in cols):
            spatial_cols = cols
            break

    if not all(c in df.columns for c in spatial_cols):
        raise ValueError(
            f"Could not find spatial coordinate columns. "
            f"Expected one of: {coord_options}"
        )

    # observation names
    if "name" in df.columns:
        obs_names = df["name"].astype(str)
    elif "cell_id" in df.columns:
        obs_names = df["cell_id"].astype(str)
    else:
        obs_names = [f"Cell_{i}" for i in range(len(df))]

    coords = df[list(spatial_cols)].to_numpy(dtype=float)

    adata = ad.AnnData(
        X=np.zeros((len(df), 0)),
        obs=pd.DataFrame(index=pd.Index(obs_names, name="cell_id")),
        obsm={"spatial": coords},
    )
    return adata

valid_source_landmarks_rot = pd.read_csv(valid_source_landmarks)
valid_source_landmarks_rot = valid_source_landmarks_rot.rename(columns={"x_centroid": "x", "y_centroid": "y"})
print(valid_source_landmarks_rot.head())
valid_source_landmarks_rot = np.array(valid_source_landmarks_rot[['x', 'y']])
query_img = Image.open(source_image_path)
valid_source_landmarks_rot = rotate_keypoints(valid_source_landmarks_rot, 270, image=query_img).to_numpy(dtype=float)
valid_source_landmarks_rot = pd.DataFrame(valid_source_landmarks_rot, columns=["x", "y"])
source_adata = landmark_csv_to_adata(valid_source_landmarks_rot)

# optional: save to disk if you want to pass paths into `PARAMS`
source_h5ad_path = Path(OUTPUT_ROOT) / "source_landmarks.h5ad"
source_adata.write_h5ad(source_h5ad_path)

# source_img = Image.open(source_image_path)
# inverted_source_path = OUTPUT_ROOT / "source_image_inverted.tif"
# inverted_img = ImageOps.invert(source_img)
# inverted_img = ImageEnhance.Contrast(inverted_img).enhance(2.0)
# 
# rotated_query = query_img.invert(source_img)
# rotated_query = rotated_query.transpose(method=Image.Transpose.FLIP_TOP_BOTTOM)
# source_image_path = OUTPUT_ROOT / "source_image_inverted.tif"
# rotated_query.save(str(inverted_source_path))

# align slices
print("registering images")
tracemalloc.start()
start = time.perf_counter()

PARAMS = {
    "sample_id": "XeniumtoHE",
    "target_image_path": target_image_path,
    "source_image_path": source_image_path,
    "source_h5ad_path": source_h5ad_path,
    "source_h5ad_spatial_key": "spatial",
    "Alignment_mode": "affine+bpsline",
    "device": "cpu",
}
PARAMS

result = run_image_to_image_alignment(
    data_root=DATA_DIR,
    output_root=OUTPUT_ROOT,
    **PARAMS,
)

adata = anndata.read_h5ad(result["transformed_h5ad"])
transformed_coords = adata.obsm["spatial"]

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")
