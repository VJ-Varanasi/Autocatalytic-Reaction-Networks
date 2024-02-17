
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
from numpy import linalg

def norm_dist(data):
    (N,C,R) = mat.shape

    nuc_norms= []
    fro_norms = []
    for i in range(N):
        for j in range(N):
            if i != j:
                nuc_norms.append(linalg.norm(mat[i] - mat[j] , "nuc"))
                fro_norms.append(linalg.norm(mat[i] - mat[j] , "fro"))
    
    nuc= [N, np.mean(nuc_norms), np.std(nuc_norms), np.min(nuc_norms), np.max(nuc_norms)]
    fro = [N, np.mean(fro_norms), np.std(fro_norms), np.min(fro_norms), np.max(fro_norms)]

    return nuc, fro



failure = glob.glob("Data/Matrix-*-Non-RAF.npy")
#print(failure)


non_raf_nuc = []
non_raf_fro = []
for file in failure:
    #run = file.split("-")[0][5:]
    #print(run)
    n = float(file.split("-")[1])
    f = float(file.split("-")[2])
    #print(file)
    #data = pd.to_numeric(pd.read_csv("{}".format(file)).iloc[:,0])

    mat = np.load(file)
    #print(data)
    nuc, fro = norm_dist(mat)
    non_raf_nuc.append([n,f] + nuc)
    non_raf_fro.append([n,f] + fro)

df_nuc= pd.DataFrame(non_raf_nuc, columns = ["n", "f", "N", "Mean", "Std", "Min", "Max"])
df_nuc.to_csv('Data/Non-RAF_nuc_norms.csv')

df_fro= pd.DataFrame(non_raf_fro, columns = ["n", "f", "N", "Mean", "Std", "Min", "Max"])
df_fro.to_csv('Data/Non-RAF_fro_norms.csv')


success = glob.glob("Data/Matrix-*[0-9]-RAF.npy")
#print(success)

raf_nuc = []
raf_fro = []
for file in success:

    #print(run)
    n = float(file.split("-")[1])
    f = float(file.split("-")[2])
    #print(file)
    #data = pd.to_numeric(pd.read_csv("{}".format(file)).iloc[:,0])

    mat = np.load(file)
    #print(data)
    nuc, fro = norm_dist(mat)
    raf_nuc.append([n,f] + nuc)
    raf_fro.append([n,f] + fro)

df_nuc= pd.DataFrame(raf_nuc, columns = ["n", "f", "N", "Mean", "Std", "Min", "Max"])
df_nuc.to_csv('Data/RAF_nuc_norms.csv')

df_fro= pd.DataFrame(raf_fro, columns = ["n", "f", "N", "Mean", "Std", "Min", "Max"])
df_fro.to_csv('Data/RAF_fro_norms.csv')
