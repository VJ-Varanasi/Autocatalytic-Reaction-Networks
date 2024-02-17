
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
import math



def len_R (n, k = 2):
    # Set of Molecules
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

len_r = {2: 8, 3:36, 4: 128, 5: 376, 6: 1004, 7:2528, 8 : 6096, 9:14260, 10:32668}
counts = {2: 0.35, 3: 0.25, 4: 0.45, 5: 1.12, 6:1.3, 7:1.4, 8:1.4, 10:1.4}

def get_files(string):
    # all_files = glob.glob("Data/*{}}".format(string))
    out_files = []
    # sizes = np.unique([file.split("-")[0][5:].split("_")[0] for file in all_files])
    for size in counts:
        out_files.append("Data/{}_{}{}".format(size, counts[size], string))

    
    return out_files



#failure_files = glob.glob("Data/*-Non-RAF_Perturbations.csv")
failure = get_files("-Non-RAF_Perturbations.csv")


k = len(failure)
n =math.ceil(np.sqrt(k))
fig, axs = plt.subplots(n,n, figsize =(10, 7))
fig.suptitle('Scaled Non-RAF Perturbations', fontsize=16)
fig.tight_layout()


for ind in range(k):
    j = ind % n
    i = ind // n 
    file=  failure[ind]

    run = file.split("-")[0][5:]
    #print(file)
    data = pd.to_numeric(pd.read_csv("{}".format(file)).iloc[:,0])
    
    run_n = float(run.split("_")[0])
    run_f = float(run.split("_")[1])
    if run_n in len_r:
        R = len_r[run_n]
    else:
        R = len_R(run_n)
        len_r[run_n] = R

    
    data = data / R
    
    

    axs[i,j].hist(data, color = "orange", label = "N = {}".format(len(data)))
    axs[i,j].legend()
    axs[i,j].title.set_text("n = {}, f= {}".format(run_n, run_f))


plt.savefig("Images/Scaled-Non-RAF_Perturbations_Hist.png")




#print(success)
#success = glob.glob("Data/*[0-9]-RAF_Perturbations.csv")

success = get_files("-RAF_Perturbations.csv")

k = len(success)
n = math.ceil(np.sqrt(k))
fig, axs = plt.subplots(n,n, figsize =(10, 7), constrained_layout = True)
fig.suptitle('Scaled RAF Perturbations', fontsize=16)
fig.tight_layout()

for ind in range(k):
    j = ind % n
    i = ind // n 
    file=  success[ind]

    run = file.split("-")[0][5:]

    data = pd.to_numeric(pd.read_csv("{}".format(file)).iloc[:,0])
    
    run_n = float(run.split("_")[0])
    run_f = float(run.split("_")[1])
    if run_n in len_r:
        R = len_r[run_n]
    else:
        R = len_R(run_n)
        len_r[run_n] = R

    data = data / R
    

    axs[i,j].hist(data, color = "green", label = "N = {}".format(len(data)))
    axs[i,j].legend()
    axs[i,j].title.set_text("n = {}, f= {}".format(run_n, run_f))


plt.savefig("Images/Scaled-RAF_Perturbations_Hist.png")


# fig, ax = plt.subplots(figsize =(10, 7))
# ax.hist(success_counts, color = "green", label = "Perturbations of RAF ({})".format(RAFs))
# plt.title("Perturbations of RAF ({})".format(RAFs))
# plt.show()
# plt.savefig("Images/{}-{}-RAF_Perturbations_Hist_{}.png".format(n,f,RAFs))

# fig, ax = plt.subplots(figsize =(10, 7))
# ax.hist(failure_counts, color = "orange", label= "Perturbations of Non-RAF ({})".format(N-RAFs))


# pd.DataFrame(success_counts).to_csv('Data/{}_{}-RAF_Perturbations.csv'.format(n, f), mode='a', index=False)
# pd.DataFrame(failure_counts).to_csv('Data/{}_{}-Non-RAF_Perturbations.csv'.format(n, f), mode='a', index=False)

# plt.title("Perturbations of Non-RAF ({})".format(N-RAFs))
# plt.show()
# plt.savefig("Images/{}-{}-Non-RAF_Perturbations_Hist_{}.png".format(n,f,N-RAFs))