from utils import *
from numpy import linalg as LA

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


success = glob.glob("Data/Dict-Matrix-*[0-9]-RAF.txt")

nf_pairs = []

success_dict = {}

for file in success:
    Cs = []
    n = int(file.split("-")[2])
    f = float(file.split("-")[3])

    file1 = open(file, 'r')
    Lines = file1.readlines()
    

    if (n,f) not in success_dict:
        success_dict[(n,f)] = [ast.literal_eval(line) for line in Lines] 
    else:
        success_dict[(n,f)].append([ast.literal_eval(line) for line in Lines])
    



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


keys1 = set(success_dict.keys())
keys2 = set(failure_dict.keys())

#Find the intersection of keys (common keys)
common_keys = keys1.intersection(keys2)
print(common_keys)                                   

#key = set([i for i in common_keys if i[0] == n]).pop()

for key in common_keys:
    success = success_dict[key]
    failure = failure_dict[key]
    n = key[0]
    f = key[1]

    X,F,R = create_XFR(n)

    matricies = [get_M(X,R, C) for C in success] +[get_M(X,R, C) for C in failure]
    s = len(matricies)

    dist = np.zeros((s, s))

    for i in range(s):
        for j in range(i+1, s):
            dist[i,j] = LA.norm(matricies[i] - matricies[j], "fro")
    

    k = 0
    dist_dict ={}

    while k < np.max(dist) +1:
    #while k < 5:
        #clust = np.argwhere((dist < k) & (dist > 0))
        dist_dict[k] = np.argwhere((dist < k) & (dist > 0))
        k+=1
    
    plt.figure(figsize=(15, 9))
    plt.title("({}, {}) Hierarchical Clustering Dendrogram".format(n,f))
    plt.plot(list(dist_dict.keys()), [len(dist_dict[i]) for i in list(dist_dict.keys())])
    plt.xlabel("Distance Frobenius Norm")
    plt.xlabel("Count")
    plt.savefig("Images/AggClusterDendrogram_({},{})".format(n,f))





    #vectors = [matrix.flatten() for matrix in matricies]



    # data =np.array(vectors)
  

    # n_clusters = 2  # Specify the number of clusters
    # agg_cluster = AgglomerativeClustering(n_clusters=n_clusters, linkage='single', compute_distances = True)
    # #labels = agg_cluster.fit_predict(data)
    # #print(labels)
    # #print(agg_cluster.get_params)

    # # setting distance_threshold=0 ensures we compute the full tree.
    # model = agg_cluster.fit(data)

    # plt.title("({}, {}) Hierarchical Clustering Dendrogram".format(n,f))
    # # plot the top three levels of the dendrogram
    # plot_dendrogram(model, truncate_mode="level", p=3)
    # plt.xlabel("Number of points in node (or index of point if no parenthesis)")
    # plt.savefig("Images/AggClusterDendrogram_({},{})".format(n,f))



# plt.figure(figsize=(15, 9))
# # Plot the dendrogram to visualize the hierarchy
# distance_matrix = squareform(pdist(data))
# linked = linkage(distance_matrix, method='single')

# # Customize the dendrogram
# dendrogram(linked, orientation='top', distance_sort='descending', labels=labels,
#            leaf_font_size=12, leaf_rotation=45, above_threshold_color='gray', color_threshold=25)

# plt.title('Dendrogram')
# plt.xticks(np.arange(0, len(data), step=20)) 
# plt.savefig("Images/AggClusterDendrogram")



# # Plot the data points with cluster labels
# plt.scatter(data[:, 0], data[:, 1], c=labels)
# plt.title('Agglomerative Hierarchical Clustering Result')
# plt.xlabel('Feature 1')
# plt.ylabel('Feature 2')

# legend_labels = list(set(labels))
# for label in legend_labels:
#     cluster_points = data[labels == label]
#     plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f'Cluster {label}', cmap='viridis')
# plt.legend()

# plt.savefig("Images/AggClusterScatter")