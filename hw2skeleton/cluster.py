from .utils import Atom, Residue, ActiveSite
import numpy as np

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


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
