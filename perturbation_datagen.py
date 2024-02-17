
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
        R, C = reduceR(R, C)
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
        return 1

def add_catalyst(X, R, C):
    add_new = 0
    while add_new == 0:
        reaction = random.sample(list(R.keys()),1)[0]
        if reaction %2 == 1:
            k = 1
        else:
            k = -1
        catalyst = random.sample(X,1)[0]
        
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

def remove_catalyst(C):
    reaction = random.sample(list(C), 1)[0]
    if reaction %2 == 1:
        k = 1
    else:
        k = -1
    
    if len(C[reaction]) == 1:
        del C[reaction]
        del C[reaction+ k]
  
    else:
      
        catalyst = random.sample(C[reaction], 1)[0]
        
        C[reaction].remove(catalyst)
        
        C[reaction+k].remove(catalyst)
       
    return(C)


def stability_test_count (N, f, n, t=2):
    RAFs = 0
    success_counts = []
    failure_counts = []

    if n == 2:
        t = 1


    for j in range(N):
        X,F,R = create_XFR(n, t)
        p = f/len(R)
        C = create_catalysts(X, len(R),p)
        

        raf = RAF(X.copy(), F.copy(), R.copy(), C)
        RAFs += raf
        if raf ==1:
            count = 0
            while raf == 1:
                #print("C:{}".format(C))
                
                new_C = remove_catalyst(C)
                #print("New C:{}".format(new_C))
                if new_C:
                    raf = RAF(X.copy(),F.copy(),R.copy(),new_C)
                    count += 1
                    C = new_C.copy()
                    #success_stability += RAF(X.copy(),F.copy(),R.copy(),new_C
                else:
                    count += 1
                    break
            success_counts.append(count)
        else:
            count = 0
            while raf != 1:
                #print("C:{}".format(C))
                new_C = add_catalyst(X.copy(),R.copy(), C)
                #print("New C:{}".format(new_C))
        
                raf = RAF(X.copy(),F.copy(),R.copy(),new_C)
                C = new_C.copy()
                count += 1
                    #success_stability += RAF(X.copy(),F.copy(),R.copy(),new_C)
            failure_counts.append(count)
            #print("R:{}".format(R))
            #new_C = add_catalyst(X.copy(),R.copy(), C)
            #failure_stability += RAF(X,F,R,new_C)

    pd.DataFrame(success_counts).to_csv('Data/{}_{}-RAF_Perturbations.csv'.format(n, f), mode='a', index=False)
    pd.DataFrame(failure_counts).to_csv('Data/{}_{}-Non-RAF_Perturbations.csv'.format(n, f), mode='a', index=False)


    return



N = ast.literal_eval(sys.argv[1])
n_range = ast.literal_eval(sys.argv[2])
f_range = ast.literal_eval(sys.argv[3])

fs = np.linspace(f_range[0], f_range[1], f_range[2])
ns = np.arange(n_range[0], n_range[1]+1, 1)

points = list(itertools.product(fs, ns))
#print(points)


pool = multiprocess.Pool(processes=int(os.getenv('SLURM_CPUS_ON_NODE')))

if __name__ == '__main__':
    with pool as p:
        vals = []
        
        for point in points:
            vals = vals + [(N,point[0],point[1])]
            
        p.starmap(stability_test_count, vals)
            
