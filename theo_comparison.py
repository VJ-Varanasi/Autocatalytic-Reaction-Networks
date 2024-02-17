
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


def create_XFR(n, t = 2,k = 2):

    if n ==2:
        t = 1

    X = []
    F = []
    alphabet = string.ascii_uppercase[0:k]
    for i in range(1, n+1):    
        vals = [''.join(m) for m in itertools.product(alphabet, repeat=i)]
        X = X + vals
        if i <= t:
            F = F + vals
    
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
    
    return(X,F, R)

def create_catalysts (X, react_count, p):
    C = {}
    for i in X:
        for j in range(1, react_count):
            if np.random.random(1)[0] < p:
                if j%2 == 1:
                    k = 1
                else:
                    k = -1

                if j in C.keys():
                    if i not in C[j]:
                        C[j].append(i)
                else:
                    C[j]= [i]

                if j+k in C.keys():
                    if i not in C[j+k]:
                        C[j+k].append(i)
                else:
                    C[j+k]= [i]
                
    return (C)

def closure(F, R):
    no_change = 0
    X = list(F)

    while no_change !=1:
        no_change = 1

        for i in list(R.values()):
            sufficient = 1
            for j in i[0]:
                if j not in X:
                    sufficient = 0
            if sufficient == 1:
                for k in i[1]:
                    if k not in X:
                        X.append(k)
                        no_change = 0
    return(X)


def Rsupp (R):
    supp = []
    for i in list(R.values()):
        cands = i[0] + i[1]
        for j in cands:
            if j not in supp:
                supp.append(j)
    return(supp)

def reduceR(R, C):

    catalyzed = list(C.keys())
    uncat_R = list(set(R.keys()) - set(catalyzed))


    for i in uncat_R:
        del R[i]
    
    
    no_change = 0
    while no_change != 1:
        no_change = 1

        suppR= Rsupp(R)
        Rs = list(C.keys())

        for i in Rs:
            for j in C[i]:
                if j not in suppR:
                    C[i] = C[i].remove(j)
                
                if not C[i]:
                    # print("DELETE")
                    del C[i]
                    if i in R:
                        del R[i]
                    no_change = 0
                    break
         
    return(R,C)

def reduceToF(F, R):
    W = closure(F,R)
    r_num= list(R.keys())
    
    for i in r_num:
        remove = 0
        for j in R[i][0]:
            if j not in W:
                remove = 1
                break
        if remove == 1:
            del R[i]

    
    return R


def RAF(X,F,R,C):
    X_old = X.copy()
    R_old = copy.deepcopy(R)

    i = 0 
    change = 0
    while change != 1:
        R, C_old = reduceR(R, C)
        X = closure(F,R)
        R = reduceToF(F,R)
        i= i+1

        if R != False and X != False:
            if X_old == X and R_old == R:
                change = 1
            else:
                R_old = copy.deepcopy(R)
                X_old = X.copy()
        else:
            break
    
    if not R:
        return 0
    else:
        #print(X_old)
        #print(R_old)
        #graph(X_old, F, R_old,C_old)
        return 1


def theo(f,n):
    if n not in len_r:
        Rn = len_R(n)
    else:
        Rn  = len_r[n]
    
    return 1- (1 - f/Rn)**240

def comp(N,f,n):
    raf_count = 0
    
    for i in range(N):
        X,F,R = create_XFR(n)

        p = f/len(R)
        #print(p)
        C = create_catalysts(X, len(R),p)
        raf_count += RAF(X,F,R,C)
    
    return raf_count/N



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


Ns = ast.literal_eval(sys.argv[1])
ns = ast.literal_eval(sys.argv[2])
fs = ast.literal_eval(sys.argv[3])


pool = multiprocess.Pool(processes=int(os.getenv('SLURM_CPUS_ON_NODE')))

if __name__ == '__main__':
    with pool as p:
        vals = []
        
        for i in range(len(Ns)):
            vals = vals + [(Ns[i],fs[i],ns[i])]
            
        p.starmap(comp, vals)