import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from helpers import make_dirs, get_path, set_rc_params
from helpers import load_sequence_data, iter_sequence_data, compute_relative_overlaps
from helpers import plot_bipartite_network, plot_network_histogram, compute_network

import os

set_rc_params()

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_S2')
make_dirs(fig_path)

lpad = 5   # Distance between x/y labels and the axes
A_ij, s, data_taxonomy, _, genus_id_to_name = load_sequence_data()


# Count genera
genera = len(list(iter_sequence_data(A_ij, s, data_taxonomy, tqdm_disable=True)))

n_cols = 1
n_rows = 6
n_rows_first = 4 # First plots should have fewer rows
n_plots = n_cols * n_rows

# Loop over subsets to generate figures
bins = np.linspace(0, 1, 6)
panel_id = 0
plot_id = 0
ax_id = 0
for _, A_ij_s, Bs, genus in iter_sequence_data(A_ij, s, data_taxonomy) :

    # Check if new panel should be made
    if plot_id % n_plots == 0 or (panel_id == 0 and plot_id == n_cols * n_rows_first):
        if plot_id > 0 :

            # Delete unused axes
            if ax_id < len(axes) :
                for ax in axes[ax_id:] :
                    ax.remove()

            fig.savefig(os.path.join(fig_path, f'figS2_{panel_id}.png'), bbox_inches = 'tight')
            plt.subplots_adjust(wspace=0.25)
            panel_id += 1
            plot_id = 0
            ax_id   = 0

        fig, axes = plt.subplots(nrows=n_rows, ncols=4*n_cols, figsize=(n_cols * 16, n_rows * 8/2))
        plt.subplots_adjust(wspace=0.25)
        axes = axes.flatten()

    # Get the matrices for the networks
    _, A_RMs_s, A_Bs_s = compute_network(Bs)

    # Add the plots
    plot_bipartite_network(Bs, A_ij_s, A_RMs_s, genus, scaling=0.5, ax=axes[ax_id])

    clip = 18
    data = A_ij_s.sum(axis=1)
    plot_network_histogram(data, bins=np.arange(clip+2)-0.5, clip=clip, labelpad=5, color=plt.cm.Pastel1(2), ylabel=None, ax=axes[ax_id+1])

    axes[ax_id+1].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ticks_loc = axes[ax_id+1].get_xticks().tolist()
    ticks_labels = ['{:,.0f}'.format(x) for x in ticks_loc]
    ticks_labels[-1] = ' ' + ticks_labels[-1] + '+'
    axes[ax_id+1].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    axes[ax_id+1].xaxis.set_ticklabels(ticks_labels)


    data = A_ij_s.sum(axis=0)
    data = data[data > 0]
    plot_network_histogram(data, bins=np.arange(clip+2)-0.5, clip=clip, xlabel='# bacteria with RM system', labelpad=5, color=plt.cm.Pastel1(3), ylabel=None, ax=axes[ax_id+2])

    axes[ax_id+2].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ticks_loc = axes[ax_id+2].get_xticks().tolist()
    ticks_labels = ['{:,.0f}'.format(x) for x in ticks_loc]
    ticks_labels[-1] = ' ' + ticks_labels[-1] + '+'
    axes[ax_id+2].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    axes[ax_id+2].xaxis.set_ticklabels(ticks_labels)


    data = compute_relative_overlaps(A_ij_s, Bs)
    plot_network_histogram(data, bins=np.linspace(0, 1, 11), clip=clip, xlabel='# intersection / # union', labelpad=5, color=plt.cm.Pastel1(4), ylabel=None, ax=axes[ax_id+3])
    axes[ax_id+3].set_xlim(0, 1)
    axes[ax_id+3].set_xticks(np.linspace(0, 1, 5))


    # Increment counters
    plot_id += 1
    ax_id   += 4


# Delete unused axes
if ax_id < len(axes) :
    for ax in axes[ax_id:] :
        ax.remove()

# Save final figure
fig.savefig(os.path.join(fig_path, f'figS2_{panel_id}.png'), bbox_inches = 'tight')