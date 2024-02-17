# Read Data In
# Creat Vector from Data
# Run PCA
# Run Logistic Regression
# Plot data + Logistic Regression Boundary

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
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from matplotlib.colors import ListedColormap

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

def get_M (X,R,C):
    M1 = np.zeros(len(X) * len(R))
    M1= M1.reshape((len(X), len(R)))

    X_dict = {}
    for i in range(len(X)):
        X_dict[X[i]] = i
    
    for reaction in C:
        for c in C[int(reaction)]:
            i = X_dict[c]
            j = int(reaction) -1
            M1[i,j] = 1

    return M1


def plot_LDA(n, ax, success,failure):
    # Sample labeled matrices (you should replace this with your data)
    # Each matrix is assumed to be a 2x2 matrix for simplicity
    X,F,R = create_XFR(n)
    matricies = [ get_M(X,R, C) for C in success] + [get_M(X,R, C) for C in failure]
   

    # Create labels for the matrices
    labels = [1]*len(success) + [0]*len(failure)  

    
    # Convert matrices to vectors
    vectors = [matrix.flatten() for matrix in matricies]

    # Perform Linear Discriminant Analysis (LDA)
    pca = PCA(n_components=2)  # Reduce to 2 dimensions
    X_pca = pca.fit_transform(vectors)


    # classifier = LogisticRegression()
    # classifier.fit(X_pca, labels)
    # # Plot the decision boundary

    # h = .02  # Step size in the mesh
    # x_min, x_max = X_pca[:, 0].min() - 1, X_pca[:, 0].max() + 1
    # y_min, y_max = X_pca[:, 1].min() - 1, X_pca[:, 1].max() + 1
    # xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))


    # # Predict the class for each point in the meshgrid
    # Z = classifier.predict(np.c_[xx.ravel(), yy.ravel()])
    # Z = Z.reshape(xx.shape)

    # plt.figure(figsize=(8, 6))
    # plt.contourf(xx, yy, Z, cmap=cmap_background, alpha=0.4)
    # plt.scatter(X_pca[:, 0], X_pca[:, 1], c=labels, cmap=cmap_points, edgecolor='k', s=40)
    # plt.title('Decision Boundary using PCA')
    # plt.xlabel('Principal Component 1')
    # plt.ylabel('Principal Component 2')

    # Create a color map for the plot
    # cmap_background = ListedColormap(['#FFAAAA', '#AAAAFF'])
    # cmap_points = ListedColormap(['#FF0000', '#0000FF'])

    
    colors = ['r', 'g']
    target_names = ['RAF', 'Non-RAF']

    for color, i, target_name in zip(colors, [0, 1], target_names):
        ax.scatter(X_pca[labels == i, 0], X_pca[labels == i, 1], color=color, alpha=0.8, label=target_name)


    ax.title('PCA of n = {}'.format(n))
    ax.xlabel('PCA Component 1')
    ax.ylabel('PCA Component 2')

    # Plot the decision boundary
    # coeff = pca.coef_[0]
    # intercept = lda.intercept_

    # xx = np.linspace(min(X_pca[:, 0]), max(X_pca[:, 0]), 100)
    # yy = -(coeff[0] / coeff[1]) * xx - (intercept[0] / coeff[1])

    # ax.plot(xx, yy, color='k', linestyle='--')


    return 

# success = glob.glob("Data/Dict-Matrix-*[0-9]-RAF.txt")

# nf_pairs = []

# success_dict = {}

# for file in success:
#     Cs = []
#     n = int(file.split("-")[2])
#     f = float(file.split("-")[3])

#     file1 = open(file, 'r')
#     Lines = file1.readlines()
    

#     if (n,f) not in success_dict:
#         success_dict[(n,f)] = [ast.literal_eval(line) for line in Lines] 
#     else:
#         success_dict[(n,f)].append([ast.literal_eval(line) for line in Lines])
    



failure = glob.glob("Data/Dict-Matrix-*-Non-RAF.txt")
failure_dict= {}

for file in failure:
    Cs = []
    n = int(file.split("-")[2])
    f = float(file.split("-")[3])

    file1 = open(file, 'r')
    Lines = file1.readlines()
    

    if (n,f) not in failure_dict:
        failure_dict[(n,f)] = [ast.literal_eval(line) for line in Lines] 
    else:
        failure_dict[(n,f)].append([ast.literal_eval(line) for line in Lines])


#keys1 = set(success_dict.keys())
#keys2 = set(failure_dict.keys())

# Find the intersection of keys (common keys)
#common_keys = keys1.intersection(keys2)
                                   

# N = len(common_keys)
# rows = int(n**0.5)  # Calculate the square root to get the number of rows
# cols = (N + rows - 1) // rows  # Calculate the number of columns

# # Create the subplot layout
# fig, axes = plt.subplots(rows, cols, figsize=(10, 10))

# # You can iterate through your list of keys and axes to create subplots
# for i, key in enumerate(common_keys):
#     ax = axes[i // cols, i % cols]
#     n = key[0]
#     success = success_dict[key]
#     failure = failure_dict[key]
#     plot_LDA(n, ax, success, failure)
#     # Add your plot or content here

# # Remove any empty subplots if the number of keys is not a perfect square
# for i in range(n, rows * cols):
#     fig.delaxes(axes.flatten()[i])

# plt.tight_layout()
# plt.savefig("Images/LDA")

n = ast.literal_eval(sys.argv[1])
key = set([i for i in list(failure_dict.keys()) if i[0] == n]).pop()

X,F,R = create_XFR(n)

#success = success_dict[key]
failure = failure_dict[key]
#matricies = [ get_M(X,R, C) for C in success] + [get_M(X,R, C) for C in failure]
matricies =  [get_M(X,R, C) for C in failure]

# Create labels for the matrices
#labels = [1]*len(success) + [0]*len(failure)  
labels = [0]*len(failure)  
print("Data Size")
print(n)
print(len(failure))
#print("---")

# Convert matrices to vectors
vectors = [matrix.flatten() for matrix in matricies]

data =np.array(vectors)
#print("Data")
#print(data.shape)


# Perform Linear Discriminant Analysis (LDA)
n_components = 40

# Perform PCA
pca = PCA(n_components=n_components)
principal_components = pca.fit_transform(data)
#print("PCA")
#print(principal_components.shape)

print("")
corr_mat = np.zeros((data.shape[1], n_components))

zeros = []
for i in range(data.shape[1]):
    for j in range(n_components):
       #print(i)
       #print(data[:,i])
       if sum(data[:,i]) != 0:
        corr_mat[i,j] = np.corrcoef(data[:,i], principal_components[:,j])[0,1]
       else: 
        zeros.append(i)
        corr_mat[i,j] = - 100
           



# Calculate the correlation of the principal components to each row of the original data
#correlations = np.corrcoef(data, principal_components, rowvar=False)

#print(data[:,34])

# Print the correlations
#print("Correlation of Principal Components to Each Row of the Original Data:")
max_corr = np.argmax(corr_mat, axis= 0)
print(max_corr)

for i in range(n_components):
    print("Reaction {}: {}".format(max_corr[i], corr_mat[max_corr[i], i]))

print("")
print(np.unique(zeros))



#print(correlations)

