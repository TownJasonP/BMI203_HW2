import sys
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
    test_cluster_number(cluster_by_partitioning, active_sites)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clustering = cluster_hierarchically(active_sites, 3)
    print(clustering)
    print(quality_index(clustering))
    write_clustering(sys.argv[3], clustering)
    test_cluster_number(cluster_hierarchically, active_sites)

if sys.argv[1][0:2] == '-R':
    print("Clustering using random method (control)")
    clustering = cluster_randomly(active_sites, 3)
    print(clustering)
    print(quality_index(clustering))
    write_clustering(sys.argv[3], clustering)
    test_cluster_number(cluster_randomly, active_sites)
