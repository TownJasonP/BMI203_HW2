from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_similarity():

    activesite_a = ['A','B','C']
    activesite_b = ['A','B','D']

    # update this assertion
    assert cluster.compute_jaccard_similarity(activesite_a, activesite_b) == 0.5

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_by_partitioning(active_sites) == []

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites) == []
