from .utils import Atom, Residue, ActiveSite
import numpy as np
import collections
import matplotlib.pyplot as plt

def flatten(x):
    '''
    flatten a nested list into a simple list
    (Why is this not a native python function?!)
    '''
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def compute_jaccard_similarity(site_a, site_b):
    """
    Compute the Jaccard similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    a = [i.type for i in site_a.residues]
    b = [i.type for i in site_b.residues]

    similarity = 0.0
    intersection = len(set(a) & set(b))
    union = len(set(a) | set(b))
    similarity = float(intersection)/float(union)
    return similarity

def initialize_k_clusters(data,k):
    '''
    Randomly choose starting points and assign all data to nearest starting
    points

    Input: list of data to be clustered (data), number of clusters (k)
    Output: list of k lists of data with each inner list representing a
    cluster
    '''

    n = len(data)
    clusters = [[] for i in range(k)]

    # randomly pick starting points
    center_indices = np.random.choice(range(n), size = k, replace = False)
    centers = [data[i] for i in center_indices]

    # assign all data points to nearest center
    for data_point in data:
        similarity = []
        if data_point in centers: # make sure centers are in their respective clusters
            assign_to = centers.index(data_point)
        else:
            for center in centers:
                similarity.append(compute_jaccard_similarity(data_point, center))
            assign_to = similarity.index(max(similarity))
        clusters[assign_to].append(data_point)

    return(clusters)

def update_clusters(clusters):
    '''
    Find new centers of clusters and reassign data to nearest new center

    Input: list of k lists representing clustered data
    Output: updated list of k lists representing clustered data (hopefully with
    more representative centers)
    '''

    # find new centers (largest average similarity to all other elements within a cluster)
    new_centers = []
    for cluster in clusters:
        similarities = []
        for element in list(cluster):
            sub_similarities = []
            for all_other_elements in cluster:
                x = compute_jaccard_similarity(element,all_other_elements)
                sub_similarities.append(x)
            similarities.append(np.average(sub_similarities)) # can change average to other metrics for different linkage types?

        max_sim = max(similarities)
        central_element = cluster[similarities.index(max_sim)]
        new_centers.append(central_element)


    # assign all data points to new centers
    k = len(clusters)
    data = [val for sublist in clusters for val in sublist] # flatten the list of lists
    new_clusters = [[] for i in range(k)]
    for data_point in data:
        if data_point in new_centers: # make sure centers are in their respective clusters
            assign_to = new_centers.index(data_point)
        else:
            similarity = []
            for center in new_centers:
                similarity.append(compute_jaccard_similarity(data_point, center))
            assign_to = similarity.index(max(similarity))

        new_clusters[assign_to].append(data_point)

    return(new_clusters)

def cluster_by_partitioning(active_sites, k):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # randomly initialize clusters
    clusters = initialize_k_clusters(active_sites, k)

    # iteratively update clusters until convergence
    while clusters != update_clusters(clusters) or clusters != update_clusters(update_clusters(clusters)):
        clusters = update_clusters(clusters)

    return clusters


def cluster_hierarchically(active_sites, k):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # calculate initial pairwise distances (list of n^2-n (active_site_i-active_site_j, distance) pairs)
    pairwise_distances = []
    n = len(active_sites)
    for i in range(n):
        for j in range(i+1,n):
            label = [active_sites[i],active_sites[j]]
            distance = compute_jaccard_similarity(active_sites[i],active_sites[j])
            pairwise_distances.append((label,distance))

    # find minimum distance in pairwise distances and group those objects
    max_value = max([i[1] for i in pairwise_distances])
    index = [i[1] for i in pairwise_distances].index(max_value)
    to_group = pairwise_distances[index][0]

    new_clusters = []
    for i in active_sites:
        if i not in to_group:
            new_clusters.append(i)
    new_clusters.append(to_group)

    # iteratively cluster the closest 'objects' in the cluster
    while len(new_clusters) > k:
        # calculate new pairwise distances:
        pairwise_distances = []
        n = len(new_clusters)
        for i in range(n):
            for j in range(i+1,n):
                label = [new_clusters[i],new_clusters[j]]
                flattened_cluster_i = flatten(new_clusters[i])
                flattened_cluster_j = flatten(new_clusters[j])
                distances = []
                for site_i in flattened_cluster_i:
                    for site_j in flattened_cluster_j:
                        distances.append(compute_jaccard_similarity(site_i,site_j))

                distance = (np.average(distances))
                pairwise_distances.append((label,distance))

        # find minimum distance in pairwise distances and group those objects
        max_value = max([i[1] for i in pairwise_distances])
        index = [i[1] for i in pairwise_distances].index(max_value)
        to_group = pairwise_distances[index][0]

        updated = []
        for i in new_clusters:
            if i not in to_group:
                updated.append(i)
        updated.append(to_group)

        new_clusters = updated

    output = []
    for i in new_clusters:
        output.append(flatten(i))

    return output

def cluster_randomly(active_sites, k):
    '''
    Generate randomly clustered data as a control measure

    input: list of active sites to put into k clusters
    output: list of lists containing active sites representing k clusters
    '''
    clusters = [[]for i in range(k)]
    for a_s in active_sites:
        ind = np.random.randint(0,k)
        clusters[ind].append(a_s)
    return(clusters)



def quality_index(clustering_result):
    '''
    Returns a weighted average of ratios of average intra-cluster similarity to average
    extra-cluster similarity

    input: list of clusters
    output: float representing a 'quality index' of clustering
    '''
    all_data = flatten(clustering_result)

    ratios = []
    sizes = []

    for cluster in clustering_result:
        intra_cluster_sim = []
        outside_cluster_sim = []

        if len(cluster) > 1:
            for j in enumerate(cluster):
                for i in enumerate(cluster):
                    if i[0] > j[0]:
                        intra_cluster_sim.append(compute_jaccard_similarity(i[1],j[1]))
            for j in cluster:
                for i in all_data:
                    if i not in cluster:
                        outside_cluster_sim.append(compute_jaccard_similarity(i,j))

            ins = np.average(intra_cluster_sim)
            out = np.average(outside_cluster_sim)
            #print(out)
            ratios.append(ins/out)
            sizes.append(len(cluster))

    weighted_avg_ratio = np.sum(np.multiply(ratios,sizes))/len(all_data)

    return(weighted_avg_ratio)

def test_cluster_number(clustering_method, data):
    '''
    produce an elbow plot to help determine ideal number of clusters

    input: clustering algorithm and data to be clustered
    output: scatterplot showing quality index as a function of cluster size
    '''

    k = []
    quality = []
    for j in range(20): # Repeat 5 times
        for i in range(2,20): # Check different numbers of clusters
            k.append(i)
            clusters = clustering_method(data,i)
            q = quality_index(clusters)
            quality.append(q)

    plt.figure(facecolor = 'white')
    plt.scatter(k,quality, color = 'blue', alpha = 0.4, marker = 'o')
    plt.ylim(0,max(quality))
    plt.xlim(0,21)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Quality Index')
    plt.grid()
    plt.show()
