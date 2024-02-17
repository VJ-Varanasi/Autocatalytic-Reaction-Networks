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
from numpy import linalg 
import glob


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
    X_prev = X.copy()
    R_prev = copy.deepcopy(R)
    C_prev = copy.deepcopy(C)

    i = 0 
    change = 0
    while change != 1:
        R_new, C_new = reduceR(copy.deepcopy(R_prev), copy.deepcopy(C_prev))
        X_new = closure(F,R_new)
        R_new = reduceToF(F,R_new)
        i= i+1

        if R_new != False and X_new != False:
            if X_prev == X_new and R_prev == R_new:
                change = 1
            else:
                R_prev = copy.deepcopy(R_new)
                X_prev = X_new.copy()
                C_prev= copy.deepcopy(C_new)
        else:
            break
    
    if not R:
        return 0 #, X_new
    
    else:
        #print(X_old)
        #print(R_old)
        #graph(X_old, F, R_old,C_old)
        return 1 #, X_new




# def test(N,f,n):
#     raf_count = 0
#     non_first = 0
#     for i in range(N):
#         X,F,R = create_XFR(n)
#         first_order = {}
#         for i in F:
#             for j in F:
#                 add = i + j
#                 if add not in F:
#                     if add not in first_order:
#                         first_order[add] = 1
#                     else:
#                         first_order[add] +=1


#         p = f/len(R)
#         #print(p)
#         C = create_catalysts(X, len(R),p)
#         raf = RAF(X,F,R,C)
#         if raf == 1:
#             raf_count +=1
#             gen = [i for i in X_old if i not in F]
#             #print(gen)
#             if len(gen) != 0:
#                 first_order_gen = [j for j in gen if j in list(first_order.keys())]
#                 if len(first_order_gen) == 0:
#                     print(C)
#                     #graph(X_old, F,R,C)
#                     non_first +=1

#     print(non_first)
#     print(raf_count)
#     print("--")
#     return 



def new_RAF (F,R,C):

    first_order = {}
    for i in F:
        for j in F:
            add = i + j
            if add not in F:
                if add not in first_order:
                    first_order[add] = 1
                else:
                    first_order[add] +=1

    react_list = []
    for reaction in R:
        if reaction % 2 == 1:
            if R[reaction][1][0] in F:
                react_list.append(reaction)
                react_list.append(reaction +1)
            elif R[reaction][1][0] in first_order.keys():
                react_list.append(reaction)
                react_list.append(reaction +1)

    count = 0
    
    for reaction in react_list:
        if reaction in C:
            if len(list(set(C[reaction]) & set(F))) != 0:
               #print(reaction)
                count += len(list(set(C[reaction]) & set(F))) 
    
    return count
 


failure = glob.glob("Data/Dict-Matrix-*-Non-RAF.txt")

df_data = []

for file in failure:
    #Cs = []
    n = int(file.split("-")[2])
    f = float(file.split("-")[3])

    file1 = open(file, 'r')
    Lines = file1.readlines()
    
    Cs = [ast.literal_eval(line) for line in Lines]

    X,F,R = create_XFR(n)

    # if (n,f) not in failure_dict:
    #     failure_dict[(n,f)] = [ast.literal_eval(line) for line in Lines] 
    # else:
    #     failure_dict[(n,f)].append([ast.literal_eval(line) for line in Lines])

    count = 0
    for c in Cs:
        if new_RAF(F,R,c) != 0:
            count += 1

    df_data.append([n,f, count])  

print("NON-RAF")
non_raf = pd.DataFrame(df_data, columns = ["n", "f", "Non-First-Order Non-RAF"])
print(non_raf)


    

success = glob.glob("Data/Dict-Matrix-*[0-9]-RAF.txt")

df_data= []
for file in success:
    n = int(file.split("-")[2])
    f = float(file.split("-")[3])

    file1 = open(file, 'r')
    Lines = file1.readlines()
    
    Cs = [ast.literal_eval(line) for line in Lines]

    X,F,R = create_XFR(n)

    # if (n,f) not in failure_dict:
    #     failure_dict[(n,f)] = [ast.literal_eval(line) for line in Lines] 
    # else:
    #     failure_dict[(n,f)].append([ast.literal_eval(line) for line in Lines])

    count = 0
    for c in Cs:
        if new_RAF(F,R,c) == 0:
            count += 1

    df_data.append([n,f, count])  

print("RAF")
raf = pd.DataFrame(df_data, columns = ["n", "f", "Non-First-Order RAF Set"])
print(raf)

# Ns = ast.literal_eval(sys.argv[1])
# ns = ast.literal_eval(sys.argv[2])
# fs = ast.literal_eval(sys.argv[3])


# pool = multiprocess.Pool(processes=int(os.getenv('SLURM_CPUS_ON_NODE')))

# if __name__ == '__main__':
#     with pool as p:
#         vals = []
        
#         for i in range(len(Ns)):
#             vals = vals + [(Ns[i],fs[i],ns[i])]
            
#         p.starmap(test, vals)