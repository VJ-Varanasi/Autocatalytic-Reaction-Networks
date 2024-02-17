import numpy as np
import pandas as pd
import itertools
import random
import string
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy
import multiprocess
import sys
import os
import ast
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from matplotlib import cm

def len_R (n, k = 2):
    len_r = {2: 8, 3:36, 4: 128, 5: 376, 6: 1004, 7:2528, 8 : 6096, 9:14260, 10:32668}
    # Set of Molecules
    if n not in len_r:
        X = []
        alphabet = string.ascii_uppercase[0:k]
        for i in range(1, n+1):    
            vals = [''.join(m) for m in itertools.product(alphabet, repeat=i)]
            X = X + vals
            
        #print(X)
        #print(F)

        # Reaction (pair of molecules)
        R = {}
        react_count = 1
        for i in range(len(X)):
            cand = X[i] + X[i]
            if len(cand) <= n:
                R[react_count] = [[X[i], X[i]], [cand]]
                react_count +=1
                #Lysis Reaction
                R[react_count] = [[cand],[X[i], X[i]]]
                react_count +=1

            for j in range(i+1, len(X)):
                cand1 = X[i] + X[j]
                cand2 = X[j] + X[i]
                if len(cand1) <= n:
                    if cand2 != cand1:
                        #print(list(R.values()))
                        if [[X[j], X[i]], [cand1]] not in list(R.values()):
                            R[react_count] = [[X[i], X[j]], [cand1]]
                            react_count +=1
                        
                        if [[cand1],[X[j], X[i]]] not in list(R.values()):
                            R[react_count] = [[cand1],[X[i], X[j]]]
                            react_count +=1
                        
                        if [[X[i], X[j]], [cand2]] not in list(R.values()):
                            R[react_count] = [[X[j], X[i]], [cand2]]
                            react_count +=1
                        
                        if [[cand2],[X[i], X[j]]] not in list(R.values()):
                            R[react_count] = [[cand2],[X[j], X[i]]]
                            react_count +=1
                    else:
                        if [[X[j], X[i]], [cand1]] not in list(R.values()):
                            R[react_count] = [[X[i], X[j]], [cand1]]
                            react_count +=1
                        

                        if [[cand1],[X[j], X[i]]] not in list(R.values()):
                            R[react_count] = [[cand1],[X[i], X[j]]]
                            react_count +=1
        return react_count -1
    else:
        return len_r[n]

len_r = {2: 8, 3:36, 4: 128, 5: 376, 6: 1004, 7:2528, 8 : 6096, 9:14260, 10:32668}

non_raf = pd.read_csv("Data/Non-RAF_Perturbations_MLE.csv")
raf = pd.read_csv("Data/RAF_Perturbations_MLE.csv")

non_raf['Size'] = non_raf.apply(lambda x: len_R(x["n"]), axis=1)
raf['Size'] = raf.apply(lambda x: len_R(x["n"]), axis=1)

non_raf["Scaled Exp MLE"] = non_raf["Exp MLE"] / non_raf.Size
raf["Scaled Exp MLE"] = raf["Exp MLE"] / raf.Size
# fig = plt.figure()
# ax = plt.axes(projection='3d')

# angles = [20, 45, 70]

# fig = plt.figure(figsize=(15, 5))
# for i in range(len(angles)):

#     ax = fig.add_subplot(1, 3, i+1, projection='3d')

#     t= ax.plot_trisurf(non_raf.n, non_raf.f, non_raf["Scaled Exp MLE"],cmap = cm.get_cmap('viridis', 12))
#     ax.set_xlabel("Network Size (n)")
#     ax.set_ylabel("Expected Number of Catalyzations (f)")
#     ax.set_zlabel("Estimated MLE")
#     ax.view_init(elev=20, azim=angles[i])

# fig.suptitle('Scaled MLE Surface of Non-RAF Perturbation')
# fig.colorbar(t, ax = ax, shrink = 0.5, aspect = 5)
# fig.tight_layout()

# plt.savefig("Images/Scaled_MLE_Surface_Non-RAF_Perturbation.png")

fig = plt.figure(figsize=(15, 5))
# for i in range(len(angles)):

#     ax = fig.add_subplot(1, 3, i+1, projection='3d')
#     t= ax.plot_trisurf(non_raf.n, non_raf.f, non_raf["Exp MLE"],cmap = cm.get_cmap('viridis', 12))
#     ax.set_xlabel("Network Size (n)")
#     ax.set_ylabel("Expected Number of Catalyzations (f)")
#     ax.set_zlabel("Estimated MLE")
#     ax.view_init(elev=20, azim=angles[i])

# fig.suptitle('MLE Surface of Non-RAF Perturbation')
# fig.colorbar(t, ax = ax, shrink = 0.5, aspect = 5)
# fig.tight_layout()

# plt.savefig("Images/MLE_Surface_Non-RAF_Perturbation.png")



fig = plt.figure(figsize=(15, 5))
plt.title('Non-RAF Perturbation MLE Contour Plot')
plt.xlabel("Network Size (n)")
plt.ylabel("Expected Number of Catalyzations (f)")
contour= plt.tricontourf(non_raf.n,non_raf.f, non_raf["Exp MLE"], cmap='viridis')
colorbar = plt.colorbar(contour)  # Adding a color bar
colorbar.set_label('Perturbation MLE values')
plt.savefig("Images/MLE_Contour_Non-RAF_Perturbation.png")



fig = plt.figure(figsize=(15, 5))
plt.title('RAF Perturbation MLE Contour Plot')
plt.xlabel("Network Size (n)")
plt.ylabel("Expected Number of Catalyzations (f)")
contour= plt.tricontourf(raf.n,raf.f, raf["Exp MLE"], cmap='viridis')
colorbar = plt.colorbar(contour)  # Adding a color bar
colorbar.set_label('Perturbation MLE values')
plt.savefig("Images/MLE_Contour_RAF_Perturbation.png")


# fig = plt.figure()
# ax = plt.axes(projection='3d')

# angles = [30, 60, 90]

# fig = plt.figure(figsize=(15, 5))
# for i in range(len(angles)):
#     ax = fig.add_subplot(1, 3, i+1, projection='3d')
#     t=ax.plot_trisurf(raf.n, raf.f, raf["Exp MLE"], cmap = cm.get_cmap('viridis', 12))
#     #ax.set_title('MLE Surface of RAF Perturbation')
#     ax.set_xlabel("Network Size (n)")
#     ax.set_ylabel("Expected Number of Catalyzations (f)")
#     ax.set_zlabel("Estimated MLE")
#     ax.view_init(elev=10, azim=angles[i])

# fig.suptitle('MLE Surface of RAF Perturbation')
# fig.colorbar(t, ax = ax, shrink = 0.5, aspect = 5)
# fig.tight_layout()

# plt.savefig("Images/MLE_Surface_RAF_Perturbation.png")

# fig = plt.figure(figsize=(15, 5))
# for i in range(len(angles)):
#     ax = fig.add_subplot(1, 3, i+1, projection='3d')
#     t=ax.plot_trisurf(raf.n, raf.f, raf["Scaled Exp MLE"], cmap = cm.get_cmap('viridis', 12))
#     #ax.set_title('MLE Surface of RAF Perturbation')
#     ax.set_xlabel("Network Size (n)")
#     ax.set_ylabel("Expected Number of Catalyzations (f)")
#     ax.set_zlabel("Estimated MLE")
#     ax.view_init(elev=10, azim=angles[i])

# fig.suptitle('Scaled MLE Surface of RAF Perturbation')
# fig.colorbar(t, ax = ax, shrink = 0.5, aspect = 5)
# fig.tight_layout()

# plt.savefig("Images/Scaled_MLE_Surface_RAF_Perturbation.png")

