
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

def exp_MLE(data):
    return len(data)/sum(data)

failure = glob.glob("Data/*-Non-RAF_Perturbations.csv")
#print(failure)



non_raf = []
for file in failure:
    run = file.split("-")[0][5:]
    #print(run)
    n = float(run.split("_")[0])
    f = float(run.split("_")[1])
    #print(file)
    data = pd.to_numeric(pd.read_csv("{}".format(file)).iloc[:,0])
    #print(data)
    non_raf.append([n,f, exp_MLE(data)])

df = pd.DataFrame(non_raf, columns = ["n", "f", "Exp MLE"])
print(df)
df.to_csv('Data/Non-RAF_Perturbations_MLE.csv')


success = glob.glob("Data/*[0-9]-RAF_Perturbations.csv")
#print(success)

raf = []
for file in success:
    run = file.split("-")[0][5:]
    #print(run)
    n = float(run.split("_")[0])
    f = float(run.split("_")[1])
    data = pd.to_numeric(pd.read_csv("{}".format(file)).iloc[:,0])
    #print(data)
    raf.append([n,f, exp_MLE(data)])

df = pd.DataFrame(raf, columns = ["n", "f", "Exp MLE"])
print("---")
print(df)
df.to_csv('Data/RAF_Perturbations_MLE.csv')
