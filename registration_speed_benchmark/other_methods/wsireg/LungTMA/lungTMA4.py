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

# files
sample = "lungTMA4"
output_dir = "results/" + sample + "/"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
else:
    shutil.rmtree(output_dir)
    os.mkdir(output_dir)
query_img = Image.open('../../../../../data/lungTMA/lungTMA_Assay4.tif')
ref_img = Image.open('../../../../../data/lungTMA/lungTMA_4.jpg')
cells = pd.read_csv("../../../../../data/lungTMA/lungTMA_Assay4coords_unrotated.csv")
print(cells.head())
scaleparam = 0.2125 * math.pow(2,3-1)
cells[["x"]] = cells[["x"]] * scaleparam
cells[["y"]] = cells[["y"]] * scaleparam
print(cells.head())
cells = np.array(cells[['x', 'y']])


# initialize registration graph
reg_graph = WsiReg2D("my_reg_project", output_dir= output_dir)

# add registration images (modalities)
print("writing images")
rotated = query_img.rotate(270, expand=True)
query_img_size = rotated.size
rotated.save(output_dir + "input.tiff")
reg_graph.add_modality(
    "modality_fluo",
    output_dir + "input.tiff",
    image_res=1,
    preprocessing={
        "image_type": "BF",
        "as_uint8": True,
        "invert_intensity": True,
    },
)
original_size = ref_img.size
rotated_ref_img = ref_img.resize(query_img_size) # why did I do this? I think it is not necessary.
rotated_ref_img.save(output_dir + 'output.tiff')
reg_graph.add_modality(
    "modality_brightfield",
    output_dir + "output.tiff",
    image_res=1,
    preprocessing={
        "image_type": "FL", 
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