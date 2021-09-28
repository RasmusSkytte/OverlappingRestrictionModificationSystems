import numpy as np
import matplotlib.pyplot as plt

from helpers import make_dirs, get_path, set_rc_params
from helpers import load_sequence_data, iter_sequence_data, generate_random_sample, compute_relative_overlaps
from helpers import plot_bipartite_network, plot_network_histogram, compute_network, conserve_RM_degree

import os

set_rc_params()

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_3')
make_dirs(fig_path)

lpad = 5   # Distance between x/y labels and the axes
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8, 8))
fig.subplots_adjust(wspace=0.25, hspace=0.1)

A_ij, s, data_taxonomy, _, genus_id_to_name = load_sequence_data()

# Determine the genera to compare
average_RM_abundence = np.array([])
for strain_id, A_ij_s, Bs, genus in iter_sequence_data(A_ij, s, data_taxonomy) :
    average_RM_abundence = np.append(average_RM_abundence, np.mean(A_ij_s.sum(axis=1)))

sort_ind = np.argsort(average_RM_abundence)

q_cuts = np.quantile(average_RM_abundence[sort_ind], (0.25, 0.50, 0.75))
q_inds = np.searchsorted(average_RM_abundence[sort_ind], q_cuts)

# Loop over subsets to generate figures
bins = np.linspace(0, 1, 6)
for d, (strain_id, A_ij_s, Bs, genus) in enumerate(iter_sequence_data(A_ij, s, data_taxonomy)) :

    # Only plot if index is in q_inds
    if d not in sort_ind[q_inds] :
        continue

    # Determine the plotting index
    k = np.where(sort_ind[q_inds] == d)[0][0]

    # Select the plotting axes
    ax = axes[k]

    # Get the matrices for the networks
    _, A_RMs_s, A_Bs_s = compute_network(Bs)

    # Get the matrices for the random networks
    B_rand = generate_random_sample(A_RMs_s, A_ij_s, conserve_RM_degree = conserve_RM_degree())
    A_ij_rand, A_RMs_rand, _ = compute_network(B_rand)

    # Add the plots
    scaling = 0.4
    plot_bipartite_network(Bs,     A_ij_s,    A_RMs_s,     genus,       scaling=scaling, ax=ax[0])
    plot_bipartite_network(B_rand, A_ij_rand, A_RMs_rand, 'Null model', scaling=scaling, ax=ax[2])

    plot_network_histogram(compute_relative_overlaps(A_ij_s,    Bs),     bins=bins,               ax=ax[1], xlabel=None, normalized=True, color=plt.cm.Dark2(0))
    plot_network_histogram(compute_relative_overlaps(A_ij_rand, B_rand), bins=bins, hatch='///',  ax=ax[1], xlabel=None, normalized=True, color='none')

    # Adjust x-position of histogram
    pos = ax[1].get_position()
    pos.x0 += 0.04
    pos.x1 += 0.04
    ax[1].set_position(pos)

    ax[1].set_xlim(0, 1)
    ax[1].set_ylim(0, 0.2)
    ax[1].set_xticks([])

    ax[1].set_ylabel('pmf.', labelpad=5)

    # Adjust the color
    for axis in ['top','bottom','left','right']:
        ax[0].spines[axis].set_linewidth(2)
        ax[0].spines[axis].set_color(plt.cm.Dark2(0))

axes[-1][1].set_xticks(bins)
axes[-1][1].set_xlabel('$I\;/\;U$')

# Finalize figures
fig.savefig(os.path.join(fig_path, 'fig3.png'), bbox_inches = 'tight')