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


def create_XFR(n, t,k = 2):

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

def graph_anlaysis (X, F, R, C):
    DG = nx.MultiDiGraph()

    for node in X:
        if node in F:
            DG.add_node(node, ncolor = "orange")
        else:
            DG.add_node(node, ncolor = 'lawngreen')

    for node in R: 
        text = "R{}".format(node)
        DG.add_node(text, ncolor = 'salmon')
        for i in range(len(R[node][0])):
            DG.add_edge(R[node][0][i], text, color = "steelblue")
        for i in range(len(R[node][1])):
            DG.add_edge(text,R[node][1][i], color ="steelblue")
    
    bet_cet_dict = nx.betweenness_centrality(DG)
    [bet_cet_dict[node] for node in bet_cet_dict if node in F ]
    bet_cet_init = np.mean([bet_cet_dict[node] for node in bet_cet_dict if node in F ])

    if nx.is_strongly_connected(DG):
        init_diameter = nx.diameter(DG)
    else:
        init_diameter = 0

    two_cycle_init = len(list(nx.simple_cycles(DG, 2)))
    par_edge_init = len(list(nx.simple_cycles(DG, 1)))
    
    
    food_edges = 0
    for j in C:
        for k in C[j]:
            if k in X and j in list(R.keys()):
                text = "R{}".format(j)
                DG.add_edge(k, text, color ="fuchsia")
                if k in F:
                    food_edges = food_edges + 1



    bet_cet_dict = nx.betweenness_centrality(DG)
    [bet_cet_dict[node] for node in bet_cet_dict if node in F ]
    bet_cet_final = np.mean([bet_cet_dict[node] for node in bet_cet_dict if node in F ])

    if nx.is_strongly_connected(DG):
        final_diameter = nx.diameter(DG)
    else:
        final_diameter = 0
    
    par_edge_final = len(list(nx.simple_cycles(DG, 1)))

    two_cycle_final = len(list(nx.simple_cycles(DG, 2)))




    ## Connected Components
    out_dict = {}
    
    #Model Params
    out_dict['Food Size'] = len(F)
    out_dict['nodes'] = DG.number_of_nodes()
    out_dict['Non-Food Catalysts'] = DG.number_of_edges() - food_edges
    out_dict['Food Catalysts'] = food_edges
    

    #out_dict['Change in Simple Cycles'] = simple_cycle_final - simple_cycle_init

    
    out_dict['Change in Diameter'] = final_diameter - init_diameter

    out_dict['Two Cycle'] = two_cycle_final - two_cycle_init
    out_dict['Parallel Edges'] = par_edge_final - par_edge_init

    out_dict ['Change in Bet_Centrality'] = bet_cet_final - bet_cet_init 

    #out_dict['coloring'] = max(nx.coloring.greedy_color(DG).values())+1


    # edge_colors = []

    # for (u,v,attrib_dict) in list(DG.edges.data()):
    #     edge_colors.append(attrib_dict['color'])
    
    # node_colors = list(nx.get_node_attributes(DG, "ncolor").values())

    # nx.draw(DG, node_color= node_colors, edge_color = edge_colors, with_labels=True, font_weight='bold', node_size = 750, pos=nx.circular_layout(DG), connectionstyle='arc3, rad = 0.1')
    # plt.show()
    return out_dict

# RAF
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

def data_gen_helper(n, f_max, f_min, t, p_count):
    output = []
    if n == 2:
        t = 1
    X,F,R = create_XFR(n, t)
    p_min = f_min / len(R)
    p_max = f_max / len(R)
    ps = np.linspace(p_min,p_max, p_count)
    for p in tqdm(ps):
        C = create_catalysts(X, len(R), p)
        analysis = graph_anlaysis(X,F,R,C)
        
        raf = RAF(X.copy(),F.copy(),R.copy(),C)
        output.append([raf, n, p] + list(analysis.values()))
        columns = ["RAF", "n", "p"] + list(analysis.keys())

    return(output, columns)



f_min = float(ast.literal_eval(sys.argv[1]))
f_max = float(ast.literal_eval(sys.argv[2]))
Ns = ast.literal_eval(sys.argv[3])
ns = ast.literal_eval(sys.argv[4])
t = int(sys.argv[5])
p_count = int(sys.argv[6])

pool = multiprocess.Pool(processes=int(os.getenv('SLURM_CPUS_ON_NODE')))

if __name__ == '__main__':
    with pool as p:
        output = []
        vals = []
        
        for i in range(len(Ns)):
            for j in range(Ns[i]):
                vals = vals + [(ns[i], f_max,f_min, t, p_count)]
            
        for result in p.starmap(data_gen_helper, vals):
            columns = result[1]
            output = output + result[0]
    
    df = pd.DataFrame(output, columns = columns)
    df.to_csv('Data/Kauffman_Data_f{}_{}.csv'.format(f_min, f_max), mode='a', index=False)
        
