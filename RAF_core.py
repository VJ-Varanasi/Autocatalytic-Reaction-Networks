


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
import json
import pickle
import ast

#Convert data to matrix
#Find raf core
#convert core to matrix
#find Fro distance

def get_M (X,R,C):

    M1 = np.zeros(len(X) * len(R))
    M1= M1.reshape((len(X), len(R)))

    X_dict = {}
    for i in range(len(X)):
        X_dict[X[i]] = i
    
    for reaction in C:
        for c in C[reaction]:
            i = X_dict[c]
            j = int(reaction) -1
            M1[i,j] = 1

    return M1


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

def RAF_C(X,F,R,C):
    
    

    X_prev = copy.deepcopy(X)
    R_prev = copy.deepcopy(R)
    C_prev = copy.deepcopy(C)

    i = 0 
    change = 0
    while change != 1:
        R_prev, C_curr = reduceR(R_prev, C_prev)
        X_curr = closure(F,R_prev)
        R_curr = reduceToF(F,R_prev)
        i= i+1

        if R_curr != False and X_curr != False:
            if X_curr == X_prev and R_curr == R_prev:
                change = 1
            else:
                R_prev = copy.deepcopy(R_curr)
                X_prev = copy.deepcopy(X_curr)
        else:
            break
    
    return C_curr

def norm_dist(C, n):

    X,F,R = create_XFR(n)

    mat_i = get_M(X,R,C)
    core = RAF_C(X,F,R,C)
    #print(core)
    #print(X)
    #print(R)
    mat_core = get_M(X,R,core)

    return linalg.norm(mat_i - mat_core , "nuc")
   


success = glob.glob("Data/Dict-Matrix-*[0-9]-RAF.txt")


raf_nuc = []
for file in success:
    norms = []
    n = int(file.split("-")[2])
    f = float(file.split("-")[3])
    #print(file)
    #data = pd.to_numeric(pd.read_csv("{}".format(file)).iloc[:,0])
    #print(file)
    file1 = open(file, 'r')
    Lines = file1.readlines()

    for line in Lines:

        
        #C = line[1:][:-2]
        #print(line)
        #print(type(line))
        #print("--")
        C= ast.literal_eval(line)
        #print(C)
        #print(type(C))
        #print("")
        norms.append(norm_dist(C,n))
        


    #with open(file, "rb") as f_dict:
        # s = "[" + f_dict.read() + "]"
        # print(s)
        # dict_list = ast.literal_eval(str(s))
        #print(s)
        #dict_list = ast.literal_eval(s)
    #dict_list = json.load(yeet)
    
    #print(dict_list)

    #result = [ast.literal_eval('{%s}' % item[1:-1]) for item in dict_list]
    #print(ast.literal_eval(dict_list))
    # for C in yeet.columns:
    #         #     norms.append(norm_dist(C,n))
    



    raf_nuc.append([n, f, len(norms), np.mean(norms), np.std(norms), np.min(norms), np.max(norms)])


    #print(data)
    # nuc, fro = norm_dist(mat)
    # raf_nuc.append([n,f] + nuc)
    # raf_fro.append([n,f] + fro)

df_nuc= pd.DataFrame(raf_nuc, columns = ["n", "f", "N", "Mean", "Std", "Min", "Max"])
df_nuc.to_csv('Data/RAF_core_nuc_norms.csv', index=False)

# df_fro= pd.DataFrame(raf_fro, columns = ["n", "f", "N", "Mean", "Std", "Min", "Max"])
# df_fro.to_csv('Data/RAF_fro_norms.csv')
