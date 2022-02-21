import numpy as np
import matplotlib.pyplot as plt

import networkx as nx
from networkx.algorithms import community

import scipy

from helpers import make_dirs, get_path, set_rc_params
from helpers import load_sequence_data, iter_sequence_data

import os

set_rc_params()


def filter_missing_RM_systems(A_ij_list) :

    # Determine RM systems that are missing
    RMs_IDs = np.argwhere(np.vstack(A_ij_list).sum(axis = 0)).flatten()

    # Remove missing RM systems
    return [A_ij_s[:, RMs_IDs] for A_ij_s in A_ij_list]




# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_S6')
make_dirs(fig_path)

# Prepare figure
fig = plt.figure(figsize=(5, 20))
ax = plt.gca()

A_ij, s, data_taxonomy, data_motifs, genus_id_to_name = load_sequence_data()

# Determine the number of RM systems in the data set
n_RM = A_ij.shape[1]


# Determine the presence-abence of the genera to compare
A_ij_truncated = []
names = []
for strain_id, _, Bs, genus in iter_sequence_data(A_ij, s, data_taxonomy) :

    A_ij_s = np.zeros((len(Bs), n_RM))
    for b_id, b in enumerate(Bs) :
        for r in b :
            A_ij_s[b_id, r] = 1

    A_ij_truncated.append(A_ij_s)
    names.append(genus)


A_ij_filtered = filter_missing_RM_systems(A_ij_truncated)

# Determine the genera where the frequency of each RM system is highest
RMs_freq = np.vstack([A_ij_s.mean(axis=0) for A_ij_s in A_ij_filtered])
dominant_genus = RMs_freq.argmax(axis=0)

# Reorder the matrix to maximize network modularity
A_ij_reordered = [A_ij_s[:, np.argsort(dominant_genus)] for A_ij_s in A_ij_filtered]

# Output example of occurence
for example_ID in np.argwhere(dominant_genus == 0).flatten() :
    print(f'RM system motif {data_motifs.loc[example_ID, "R-M motif"]}\t is most frequently found in {names[dominant_genus[example_ID]]} where it is found in {100*RMs_freq[dominant_genus[example_ID], example_ID]:.1f}\t percent of strains')


# Produce the figure
offset = 0
for genus_id, A_ij_s in enumerate(A_ij_reordered) :

    color = plt.cm.tab10(genus_id % 10)
    (I, J, _) = scipy.sparse.find(A_ij_s)

    plt.scatter(J, I+offset, marker = 's', color = color, s = 4, edgecolors = 'none')

    offset += A_ij_s.shape[0]


n_Bs = np.array(list(map(len, A_ij_truncated)))
ytick_pos = np.cumsum(n_Bs) - n_Bs/2

# Adjust the axes
ax.set_aspect('equal', 'box')

d = 1
ax.set_xlim(-d, A_ij_filtered[0].shape[1]+d)
ax.set_ylim(-d, offset+d)

ax.set_xticks([])
ax.set_xlabel('RM systems')

ax.set_yticks(ytick_pos)
ax.set_yticklabels(names, fontstyle = 'italic', verticalalignment = 'center')
ax.invert_yaxis()

# Finalize figures
fig.savefig(os.path.join(fig_path, 'figS6.png'), bbox_inches = 'tight', dpi = 300)



# Report the network modularity
strains_to_filter = [['Nothing'], ['Helicobacter']]

for strain_filter in strains_to_filter :

    # Determine IDs to filter
    filter_ids = np.argwhere(np.isin(names, strain_filter)).flatten()
    A_ij_filtered = filter_missing_RM_systems([A_ij_s for genus_id, A_ij_s in enumerate(A_ij_truncated) if genus_id not in filter_ids])

    n_Bs = np.array(list(map(len, A_ij_filtered)))

    # Report the network modularity
    G = nx.Graph()

    strain_ids = np.arange(sum(n_Bs))
    RMs_ids    = [str(r) for r in range(A_ij_filtered[0].shape[1])]

    # Determine the genera where the frequency of each RM system is highest
    RMs_freq = np.vstack([A_ij_s.mean(axis=0) for A_ij_s in A_ij_filtered])
    dominant_genus = RMs_freq.argmax(axis=0)

    # Add nodes with the node attribute "bipartite"
    G.add_nodes_from(strain_ids, bipartite=0)  # Stain nodes
    G.add_nodes_from(RMs_ids,    bipartite=1)  # RM nodes

    for strain_id, row in enumerate(np.vstack(A_ij_filtered)) :

        if strain_id in filter_ids :
            continue

        for RMs_id in np.argwhere(row).flatten() :
            G.add_edge(strain_id, str(RMs_id))

    # Partition strains by genus
    partition_strain_nodes = [set(np.arange(N) + offset) for N, offset in zip(n_Bs, np.cumsum([0] + n_Bs[:-1].tolist()))]

    # Partition RM systems by dominant genus
    partition_RM_nodes     = [set([RMs_id for g, RMs_id in zip(dominant_genus, RMs_ids) if g == genus]) for genus in range(len(n_Bs))]

    # Check modularity of partiitons
    partitions = [strains.union(RMs) for strains, RMs in zip(partition_strain_nodes, partition_RM_nodes)]

    str_strains = ', '.join(strain_filter)
    print(f'Modularity using dominant genus partitions (excluding {str_strains}): {community.modularity(G, partitions):.2f}')



