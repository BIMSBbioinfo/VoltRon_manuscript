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
import os
from STalign import STalign

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Read arguments from a YAML file.")
    parser.add_argument("yaml_file", type=str, help="Path to the YAML file containing arguments.")
    args = parser.parse_args()

    with open(args.yaml_file, 'r') as file:
        arguments = yaml.safe_load(file)

    # get points
    df = pd.read_csv(arguments["cells"])
    print("Number of Cells:", df.shape[0])
    xM = np.array(df['x'])
    yM = np.array(df['y'])
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
    pointsI = pointsI[["y", "x"]].to_numpy()
    pointsJ = pd.read_csv(arguments["affine_source_landmarks"])
    pointsJ = pointsJ[["y", "x"]].to_numpy()
    print("Number of alignment landmarks:", pointsI.shape[0])

    # directory
    if not os.path.exists("results/" + arguments["name"]):
        os.mkdir("results/" + arguments["name"])

    # save affine point comparison
    extentJ = STalign.extent_from_x((YJ,XJ))
    fig,ax = plt.subplots(1,2)
    ax[0].imshow((I.transpose(1,2,0).squeeze()), extent=extentI)
    ax[1].imshow((J.transpose(1,2,0).squeeze()), extent=extentJ)
    ax[0].scatter(pointsI[:,1],pointsI[:,0], c='red')
    ax[1].scatter(pointsJ[:,1],pointsJ[:,0], c='red')
    for i in range(pointsI.shape[0]):
        ax[0].text(pointsI[i,1],pointsI[i,0],f'{i}', c='red')
        ax[1].text(pointsJ[i,1],pointsJ[i,0],f'{i}', c='red')
    fig.savefig("results/" + arguments["name"] + "/compare_affinelandmarks.png")

    # get validation source landmarks
    Xenium_image = plt.imread(arguments["source_image"])
    Xenium_landmarks = pd.read_csv(arguments['valid_source_landmarks'])
    Xenium_landmarks.head()
    xlandxenium = np.array(Xenium_landmarks['x']) # * scaleparam
    ylandxenium = np.array(Xenium_landmarks['y']) # * scaleparam
    print("Number of source validation landmarks:", Xenium_landmarks.shape[0])

    # get validation target landmarks
    HE_landmarks = pd.read_csv(arguments['valid_target_landmarks'])
    HE_landmarks.head()
    xland = np.array(HE_landmarks['x'])
    yland = np.array(HE_landmarks['y'])
    yland =  V.shape[0] - yland

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

    # get necessary output variables
    A = out['A']
    v = out['v']
    xv = out['xv']

    # transform points
    tpointsJ = STalign.transform_points_target_to_source(xv,v,A,np.stack([yM, xM], -1))
    if tpointsJ.is_cuda:
        tpointsJ = tpointsJ.cpu()

    # transform (affine) validation landmark
    affine = np.dot(np.linalg.inv(L), [ylandxenium - T[0], xlandxenium - T[1]])
    xlandxeniumaffine = affine[0,:]
    ylandxeniumaffine = affine[1,:]
    tlandpointsJ_affine = np.stack([xlandxeniumaffine, ylandxeniumaffine], -1)

    # transform validation landmarks
    tlandpointsJ = STalign.transform_points_target_to_source(xv,v,A,np.stack([ylandxenium, xlandxenium], -1))
    if tlandpointsJ.is_cuda:
        tlandpointsJ = tlandpointsJ.cpu()

    # get distance between landmark and transformation
    diff_affine = np.sqrt(np.power(np.array(tlandpointsJ_affine[:,0])-yland,2) + np.power(np.array(tlandpointsJ_affine[:,1])-xland,2))
    diff_nonrigid = np.sqrt(np.power(np.array(tlandpointsJ[:,0])-yland,2) + np.power(np.array(tlandpointsJ[:,1])-xland,2))

    # make result data
    # pd.DataFrame(np.stack([xlandxeniumaffine, ylandxeniumaffine], -1)).to_csv("results/" + arguments["name"] + "/affine_transformed_landmarks.csv")
    # pd.DataFrame(tlandpointsJ).to_csv("results/" + arguments["name"] + "/transformed_landmarks.csv")
    # pd.DataFrame(np.stack([diff_affine, diff_nonrigid], -1)).to_csv("results/" + arguments["name"] + "/diff_landmark.csv")
    df = pd.DataFrame(np.stack([xland, yland, xlandxeniumaffine, ylandxeniumaffine, tlandpointsJ[:,0], tlandpointsJ[:,1], diff_affine, diff_nonrigid], -1))
    df.columns = ['x_landmark', 'y_landmark', 'x_stalign_affine', 'y_staling_affine', 'x_stalign', 'y_stalign', 'diff_affine', 'diff_nonrigid']
    df.to_csv("results/" + arguments["name"] + "/all_landmarks.csv")

    # visualize and save
    fig,ax = plt.subplots()
    ax.imshow((I).transpose(1,2,0),extent=extentI)
    ax.scatter(tlandpointsJ[:,1].detach(),tlandpointsJ[:,0].detach(),alpha=1, c = "red")
    ax.scatter(tlandpointsJ_affine[:,1],tlandpointsJ_affine[:,0],alpha=1, c = "yellow")
    ax.scatter(xland,yland, c='purple')
    fig.savefig("results/" + arguments["name"] + "/compare_validlandmarks.png")