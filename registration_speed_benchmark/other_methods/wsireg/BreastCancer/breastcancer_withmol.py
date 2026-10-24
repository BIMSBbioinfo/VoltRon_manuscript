from wsireg.wsireg2d import WsiReg2D
from PIL import Image
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import os
import shutil
import time
import tracemalloc
import math

Image.MAX_IMAGE_PIXELS = 665836957

# files
sample = "HE_vs_XeniumR1"
output_dir = "results/" + sample + "/"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
else:
    shutil.rmtree(output_dir)
    os.mkdir(output_dir)
query_img = Image.open('../../../../../data/Xenium_R1/image_3.tiff')
ref_img = Image.open('../../../../../data/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif')
cells = pd.read_csv("../../../../../data/Xenium_R1/transcripts.csv.gz")
print(cells.head())
scaleparam = 0.2125 * math.pow(2,3-1)
cells[["x_location"]] = cells[["x_location"]] * scaleparam
cells[["y_location"]] = cells[["y_location"]] * scaleparam
print(cells.head())
cells = np.array(cells[['x_location', 'y_location']])

# initialize registration graph
reg_graph = WsiReg2D("my_reg_project", output_dir= output_dir)

# add registration images (modalities)
rotated_query = query_img.transpose(method=Image.Transpose.FLIP_TOP_BOTTOM)
print(rotated_query.size)
query_img_size = rotated_query.size
rotated_query.save(output_dir + "input.tiff")
reg_graph.add_modality(
    "modality_fluo",
    output_dir + "input.tiff",
    image_res=1,
    preprocessing={
        "image_type": "FL",
        "as_uint8": True,
        "invert_intensity": True,
    },
)

original_size = ref_img.size
print(original_size)
rotated_ref = ref_img
rotated_ref = rotated_ref.convert("L")
rotated_ref = rotated_ref.resize(query_img_size, resample=Image.Resampling.LANCZOS)
rotated_ref.save(output_dir + 'output.tiff')
reg_graph.add_modality(
    "modality_brightfield",
    output_dir + "output.tiff",
    image_res=1,
    preprocessing={
        "image_type": "BF", 
        "as_uint8": True,
        "invert_intensity": False,
    }
)

# add points
reg_graph.add_shape_set("modality_fluo", "Xenium_landmarks", {"array": cells, "shape_type": "Point", "shape_name" : "Xenium_landmarks"}, 1)

# register
print("registering images")
tracemalloc.start()
start = time.perf_counter()

reg_graph.add_reg_path(
    "modality_fluo",
    "modality_brightfield",
    thru_modality=None,
    reg_params=["rigid", "nl"],
)
reg_graph.register_images()

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")

print("saving transformations")
reg_graph.save_transformations()

print("transform images")
tracemalloc.start()
start = time.perf_counter()

reg_graph.transform_images(file_writer="ome.tiff")

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")

print("transform points")
tracemalloc.start()
start = time.perf_counter()

reg_graph.transform_shapes()

elapsed = time.perf_counter() - start
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"    Time: {elapsed:.4f} s")
print(f"    Current memory: {current / 1024**2:.2f} MiB")
print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")
