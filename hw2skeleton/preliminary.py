import numpy as np
from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites,write_clustering, write_mult_clusterings
from .cluster import *

datapath = '/Users/student/Documents/Algorithms/BMI203_HW2/data'
active_sites = read_active_sites(datapath)

lengths = []
for i in active_sites:
    lengths.append(len(i.residues))

# plot histogram of number of residues in active sites
plt.hist(lengths, bins = np.linspace(0,50,21))
plt.xlabel('Residue Number')
plt.ylabel('Frequency')
plt.show()

# plot histogram of types of residues in active sites
residues_str = []
for i in active_sites:
    for j in i.residues:
        residues_str.append(j.type)
unique = list(set(residues_str))

counts = np.zeros(len(unique))
for i in enumerate(unique):
    for j in residues:
        if i[1] == j:
            counts[i[0]] += 1

plt.bar(left = range(len(counts)), width = 1, height = counts)
plt.xticks([i+0.5 for i in range(len(counts))],unique, rotation = 'vertical')
plt.xlim(0,17)
plt.show()

probs = [i/sum(counts) for i in counts]

# plot co-occurence of amino acids in active sites
residues = []
for i in active_sites:
    residues.append(i.residues)

result = np.zeros((len(unique),len(unique)))
for i in enumerate(unique):
    for j in enumerate(unique):
        for k in residues:
            if i[1] in [aa.type for aa in k]:
                new_list = [aa.type for aa in k]
                new_list.remove(i[1])
                if j[1] in new_list:
                    result[i[0]][j[0]] += 1

total = np.sum(np.sum(result))

expected_result = np.zeros((len(unique),len(unique)))
for i in enumerate(probs):
    for j in enumerate(probs):
        expected_result[i[0]][j[0]] = i[1]*j[1]

plt.imshow(result/total, interpolation = 'nearest', cmap = 'YlOrRd', vmin = 0, vmax = 0.05)
plt.xticks(range(17),unique, rotation = 'vertical')
plt.yticks(range(17),unique)
plt.colorbar()
plt.title('Observed Co-occurence')
plt.show()

plt.imshow(expected_result, interpolation = 'nearest', cmap = 'YlOrRd', vmin = 0, vmax = 0.05)
plt.xticks(range(17),unique, rotation = 'vertical')
plt.yticks(range(17),unique)
plt.colorbar()
plt.title('Expected Co-occurence')
plt.show()

# find 'expected' Jaccard similarity for randomly generated active sites
n = 50
matrix = np.zeros((n,n))

for sim_i in range(1,n):
    for sim_j in range(1,n):
        rand_data_a = np.random.choice(unique, p = probs, size = (100,sim_i))
        rand_data_b = np.random.choice(unique, p = probs, size = (100,sim_j))

        sims = []
        for i in enumerate(rand_data_a):
            for j in enumerate(rand_data_b):
                sims.append(compute_jaccard_similarity(i[1],j[1]))
        matrix[sim_i][sim_j] = np.average(sims)
