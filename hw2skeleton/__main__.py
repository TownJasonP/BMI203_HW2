import sys
import itertools
from .io import read_active_sites,write_clustering, write_mult_clusterings
from .cluster import *#compute_jaccard_similarity,cluster_by_partitioning, cluster_hierarchically, initialize_k_clusters, update_clusters
import matplotlib.pyplot as plt
import numpy as np

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites, 3)
    print(quality_index(clustering))
    write_clustering(sys.argv[3], clustering)
    #test_cluster_number(cluster_by_partitioning, active_sites)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clustering = cluster_hierarchically(active_sites, 3)
    print(clustering)
    print(quality_index(clustering))
    write_clustering(sys.argv[3], clustering)


if sys.argv[1][0:2] == '-R':
    print("Clustering using random grouping method (control)")
    clustering = cluster_randomly(active_sites, 3)
    print(clustering)
    print(quality_index(clustering))
    write_clustering(sys.argv[3], clustering)
    #test_cluster_number(cluster_randomly, active_sites)

# compare clustering results
def compute_overlap(a,b):
    # ONLY VALID IF MULTIPLES ARE NOT EXPECTED WITHIN a OR WITHIN b!!!
    count = 0
    for i in a:
        if i in b:
            count = count + 1
    return(count)

def compare_cluster_overlap(n):
    k_clusters = cluster_by_partitioning(active_sites, n)
    h_clusters = cluster_hierarchically(active_sites, n)

    permutations = list(itertools.permutations(range(n)))

    overlaps = []
    for i in permutations:
        perm_overlap = []
        for j in range(n):
            perm_overlap.append(compute_overlap(k_clusters[j],h_clusters[i[j]]))
        overlaps.append(np.sum(perm_overlap))
    return float(max(overlaps))

def compare_cluster_overlap_rand(n):
    k_clusters = cluster_randomly(active_sites, n)
    h_clusters = cluster_randomly(active_sites, n)

    permutations = list(itertools.permutations(range(n)))

    overlaps = []
    for i in permutations:
        perm_overlap = []
        for j in range(n):
            perm_overlap.append(compute_overlap(k_clusters[j],h_clusters[i[j]]))
        overlaps.append(np.sum(perm_overlap))
    return float(max(overlaps))

plt.figure(facecolor = 'white')
x = []
y = []
y_c = []
for j in range(5): # number sims to run
    print('Simulation # ' + str(j+1) + ' of 5')
    for i in range(2,8): # range of cluster numbers to test
        x.append(i)
        y.append(compare_cluster_overlap(i)/len(active_sites))
        y_c.append(compare_cluster_overlap_rand(i)/len(active_sites))

plt.ylabel('Overlap of Cluster Contents')
plt.xlabel('Number of Clusters')
plt.scatter(x,y, marker = 'o', color = 'blue', alpha = 0.4, label = 'H, P')
plt.scatter(x,y_c, marker = 'o', color = 'green', alpha = 0.4, label = 'R, R')
plt.ylim(0,1)
plt.xlim(0,10)
plt.legend()
plt.show()

# Compare clustering quality across methods and across cluster number
plt.figure(facecolor = 'white')
for i in [cluster_hierarchically, cluster_by_partitioning, cluster_randomly]:
    methods = {'cluster_randomly':'Random', 'cluster_hierarchically':'Hierarchical', 'cluster_by_partitioning':'Partition'}
    color_dict = {'cluster_randomly':'green', 'cluster_hierarchically':'red', 'cluster_by_partitioning':'blue'}

    k, q = test_cluster_number(i, active_sites, 5)
    plt.scatter(k,q, alpha = 0.4, marker = 'o', color = color_dict[i.__name__], label = methods[i.__name__])

plt.xlim(0,20)
plt.xlabel('Number of Clusters')
plt.ylabel('Quality Index')
plt.grid()
plt.legend()
plt.savefig('cluster_quality.png')
