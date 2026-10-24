''' STAlign
python stalign_template.py config.yml
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
import math
import argparse
import yaml
from STalign import STalign
from PIL import Image
Image.MAX_IMAGE_PIXELS = 665836957
import tifffile
import time
import tracemalloc

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Read arguments from a YAML file.")
    parser.add_argument("yaml_file", type=str, help="Path to the YAML file containing arguments.")
    args = parser.parse_args()

    with open(args.yaml_file, 'r') as file:
        arguments = yaml.safe_load(file)

    # get points
    df = pd.read_csv(arguments["cells"])
    print("Number of Cells:", df.shape[0])
    print(df.head())
    xM = np.array(df['x_centroid'])
    yM = np.array(df['y_centroid'])
    scaleparam = 0.2125 * math.pow(2,3-1)
    xM = xM * scaleparam
    yM = yM * scaleparam
    XJ,YJ,M,fig = STalign.rasterize(xM, yM, dx=30)
    J = np.vstack((M, M, M)) # make into 3xNxM
    J = STalign.normalize(J)

    # get target image
    with tifffile.TiffFile(arguments["target_image"]) as tif:
        base = tif.pages[0]
        # every top-level page the same size as the first is one channel (this skips thumbnails)
        channel_pages = [p for p in tif.pages if p.shape == base.shape]
        print(len(channel_pages), "channel page(s); page 0 shape:", base.shape)
        channels = [[p, *p.pages][2].asarray() for p in channel_pages]   # third resolution of each channel
    V = channels[0] if len(channels) == 1 else np.stack(channels, axis=-1)   # (H, W, 3), like plt.imread
    print("Image size:", V.shape)

    # process image
    Inorm = STalign.normalize(V)
    I = Inorm.transpose(2,0,1)
    # I = Inorm
    YI = np.array(range(I.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    XI = np.array(range(I.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    extentI = STalign.extent_from_x((YI,XI))
    extentI

    # get affine landmarks
    pointsI = pd.read_csv(arguments["affine_target_landmarks"])
    pointsI = pointsI[["y", "x"]].to_numpy()
    pointsJ = pd.read_csv(arguments["affine_source_landmarks"])
    pointsJ = pointsJ[["y", "x"]].to_numpy()
    print("Number of alignment landmarks:", pointsI.shape[0])

    # register
    print("registering points")
    tracemalloc.start()
    start = time.perf_counter()
    
    # affine transformation
    L,T = STalign.L_T_from_points(pointsI,pointsJ)
    affine = np.dot(np.linalg.inv(L), [yM - T[0], xM - T[1]])
    xMaffine = affine[0,:]
    yMaffine = affine[1,:]

    # run stalign
    if torch.cuda.is_available():
        device = 'cuda:0'
    else:
        device = 'cpu'
    params = {'L':L,'T':T,
            'niter': arguments["niter"],
            'pointsI':pointsI,
            'pointsJ':pointsJ,
            'device':device,
            'sigmaM':arguments["sigmaM"],
            'sigmaB':arguments["sigmaB"],
            'sigmaA':arguments["sigmaA"],
            'epV':arguments["epV"],
            'muB': torch.tensor([0,0,0]), # black is background in target
            'muA': torch.tensor([1,1,1]) # use white as artifact
            }
    out = STalign.LDDMM([YI,XI],I,[YJ,XJ],J,**params, )

    elapsed = time.perf_counter() - start
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print(f"    Time: {elapsed:.4f} s")
    print(f"    Current memory: {current / 1024**2:.2f} MiB")
    print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")

    print("transform points")
    tracemalloc.start()
    start = time.perf_counter()

    # get necessary output variables
    A = out['A']
    v = out['v']
    xv = out['xv']

    # transform points
    tpointsJ = STalign.transform_points_target_to_source(xv,v,A,np.stack([yM, xM], -1))

    elapsed = time.perf_counter() - start
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print(f"    Time: {elapsed:.4f} s")
    print(f"    Current memory: {current / 1024**2:.2f} MiB")
    print(f"    Peak memory:    {peak / 1024**2:.2f} MiB")