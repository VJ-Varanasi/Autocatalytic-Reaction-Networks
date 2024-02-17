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
import json
import pickle
from utils import *


def add_catalyst(F, R, C):
    add_new = 0
    while add_new == 0:
        reaction = random.sample(list(R.keys()),1)[0]
        if reaction %2 == 1:
            k = 1
        else:
            k = -1
        catalyst = random.sample(F,1)[0]
        
        if reaction in C:
            if catalyst not in C[reaction]:
                C[reaction].append(catalyst)
                C[reaction + k].append(catalyst)
                add_new = 1

        else:
            C[reaction] = [catalyst]
            C[reaction+k] = [catalyst]
            add_new = 1
    return(C)

def remove_catalyst(C, F):
    #print(np.unique([num for sublist in list(C.values()) for num in sublist]))
    food_catalysts = list(set(np.unique([num for sublist in list(C.values()) for num in sublist])) & set(F))
    #print("FC: {}".format(food_catalysts))
    if len(food_catalysts) > 0:
        catalyst = random.sample(food_catalysts, 1)[0]

        for i in list(C.keys()):
            
            if catalyst in C[i]:
                if i %2 == 1:
                    k = 1
                else:
                    k = -1

                if len(C[i]) > 1:
                    C[i].remove(catalyst)
                    C[i+k].remove(catalyst)
                else:
                    del C[i]
                    del C[i+k]
                
                break
 
    return C

def get_M (X,R,C,X_dict):

    M1 = np.zeros(len(X) * len(R))
    M1= M1.reshape((len(X), len(R)))
    
    for reaction in C:
        for c in C[reaction]:
            i = X_dict[c]
            j = reaction -1
            M1[i,j] = 1

    return M1

def sparse_gen (N, f, n, t=2):
    success = []
    failure = []

    if n == 2:
        t = 1

    for j in range(N):
        X,F,R = create_XFR(n, t)
        # X_dict = {}
        # for i in range(len(X)):
        #     X_dict[X[i]] = i
        p = f/len(R)
        C = create_catalysts(X, len(R),p)
        #print(C)
    
        raf = RAF(X, F, R, C)
        if raf == 1:
            success.append(C)
        else:
            failure.append(C)
    
    print(len(success)/N)
    
    success_name = 'Data/Dict-Matrix-{}-{}-RAF.txt'.format(n,f)
    #print(success)
    with open(success_name,'a+') as file:  
        #for i in success:
        for C in success:
            file.write(str(C))
            file.write("\n")
        #json.dump(success, file)
        #pickle.dump(success, file)
  
   
     
    failure_name = 'Data/Dict-Matrix-{}-{}-Non-RAF.txt'.format(n,f)
    with open(failure_name,'a+') as file:  
        for C in failure:
            file.write(str(C))
            file.write("\n")
        #json.dump(failure, file)
            #file.write(pickle.dumps(i))
    

    return


Ns = ast.literal_eval(sys.argv[1])
ns = ast.literal_eval(sys.argv[2])
fs = ast.literal_eval(sys.argv[3])

# for i in range(len(Ns)):
#     sparse_gen(Ns[i], fs[i], ns[i])


#             vals = vals + [(Ns[i],fs[i],ns[i])]




pool = multiprocess.Pool(processes=int(os.getenv('SLURM_CPUS_ON_NODE')))

if __name__ == '__main__':
    with pool as p:
        vals = []
        
        for i in range(len(Ns)):
            vals = vals + [(Ns[i],fs[i],ns[i])]
            
        p.starmap(sparse_gen, vals)