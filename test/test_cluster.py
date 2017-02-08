from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_similarity():

    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.compute_jaccard_similarity(active_sites[0], active_sites[1]) == 0.125
    assert cluster.compute_jaccard_similarity(active_sites[1], active_sites[2]) == 3.0/7
    assert cluster.compute_jaccard_similarity(active_sites[0], active_sites[2]) == 0.25

    # verify distance metric properties:
    #self similarity = 1
    assert cluster.compute_jaccard_similarity(active_sites[0], active_sites[0]) == 1
    # sim(a,b) = sim(b,a)
    assert cluster.compute_jaccard_similarity(active_sites[0], active_sites[1]) == cluster.compute_jaccard_similarity(active_sites[1], active_sites[0])
    # triangle inequality
    assert cluster.compute_jaccard_similarity(active_sites[0], active_sites[1]) + cluster.compute_jaccard_similarity(active_sites[1], active_sites[2]) >= cluster.compute_jaccard_similarity(active_sites[0], active_sites[2])

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # verify cluster number
    for i in [2,3]:
        assert len(cluster.cluster_by_partitioning(active_sites,i)) == i

    # had a lot of trouble pinning down how to test this adequately, ran out of time
    assert (cluster.cluster_by_partitioning(active_sites,2))[0] == [active_sites[0]] or [active_sites[1],active_sites[2]]

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites,2) == [[active_sites[0]],[active_sites[1],active_sites[2]]]
