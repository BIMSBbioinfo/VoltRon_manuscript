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
    xM = np.array(df['x_centroid'])
    yM = np.array(df['y_centroid'])
    XJ,YJ,M,fig = STalign.rasterize(xM, yM, dx=30)
    J = np.vstack((M, M, M)) # make into 3xNxM
    J = STalign.normalize(J)

    # get target image
    V = plt.imread(arguments["target_image"])
    print("Image size:", V.shape)

    # process image
    Inorm = STalign.normalize(V)
    I = Inorm.transpose(2,0,1)
    YI = np.array(range(I.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    XI = np.array(range(I.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    extentI = STalign.extent_from_x((YI,XI))
    extentI

    # get affine landmarks
    pointsI = pd.read_csv(arguments["affine_target_landmarks"])
    pointsI = pointsI.to_numpy()
    pointsJ = pd.read_csv(arguments["affine_source_landmarks"])
    scaleparam = 0.2125 * math.pow(2,3-1)
    pointsJ["x"] = np.array(pointsJ['x']) * scaleparam
    pointsJ["y"] = np.array(pointsJ['y']) * scaleparam
    pointsJ = pointsJ.to_numpy()
    print("Number of alignment landmarks:", pointsI.shape[0])

    # register
    print("registering points")
    tracemalloc.start()
    start = time.perf_counter()

    # affine transformation
    L,T = STalign.L_T_from_points(np.flip(pointsI, axis = 1),np.flip(pointsJ, axis = 1))
    affine = np.dot(np.linalg.inv(L), [yM - T[0], xM - T[1]])
    xMaffine = affine[0,:]
    yMaffine = affine[1,:]

    # run stalign
    if torch.cuda.is_available():
        device = 'cuda:0'
    else:
        device = 'cpu'
    pointsI = np.asarray(np.flip(pointsI, axis = 1)).copy()
    pointsJ = np.asarray(np.flip(pointsJ, axis = 1)).copy()
    print(pointsI, pointsJ)
    params = {'L':L,'T':T,
            'niter': arguments["niter"],
            'pointsI': pointsI,
            'pointsJ': pointsJ,
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
