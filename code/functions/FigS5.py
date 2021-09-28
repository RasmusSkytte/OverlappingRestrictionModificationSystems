import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from helpers import make_dirs, get_path, set_rc_params
from helpers import load_sequence_data, iter_sequence_data, compute_relative_overlaps
from helpers import add_figure_labels

import os

set_rc_params()

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_S5')
make_dirs(fig_path)


A_ij, s, data_taxonomy, data_motifs, genus_id_to_name = load_sequence_data()


df_overlap = pd.DataFrame({'genus' : [], 'genus_id' : [], 'nRM_i' : [], 'nRM_j' : [], 'o_ij' : []})

# Loop over subsets to generate figures
for genus_id, A_ij_s, Bs, genus in iter_sequence_data(A_ij, s, data_taxonomy) :

    # Get the overlap for the subset
    o_ij_s = compute_relative_overlaps(A_ij_s, Bs, ind='all')

    # Exclude diagnoal elements
    triu_ind_x, triu_ind_y = np.triu_indices(len(Bs), k=1)
    tril_ind_x, tril_ind_y = np.tril_indices(len(Bs), k=-1)
    ind = (np.concatenate((triu_ind_x, tril_ind_x)), np.concatenate((triu_ind_y, tril_ind_y)))

    # Store the number of RM systems in each member of the pair
    nRM_i = np.array([list(map(len, Bs))] * len(Bs))
    nRM_j = nRM_i.T

    df_overlap_s = pd.DataFrame({
                        'genus' : genus,
                        'genus_id' : genus_id,
                        'nRM_i' : nRM_i[ind],
                        'nRM_j' : nRM_j[ind],
                        'o_ij' : o_ij_s[ind]})

    df_overlap = df_overlap.append(df_overlap_s)


fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 6))
plt.subplots_adjust(wspace=0.3)

for ax_id, strains_to_remove in enumerate(['', 'Helicobacter']) :

    # Filter out strains
    I = np.isin(df_overlap.genus, strains_to_remove)
    df = df_overlap.loc[~I]

    # Take the average over the different genera
    df = df.groupby(['nRM_i', 'nRM_j']).o_ij.mean().reset_index()

    # Convert to 2D representation
    df = pd.pivot_table(df, index='nRM_i', columns='nRM_j', values='o_ij')
    df.fillna(0, inplace=True)
    index = pd.Index(np.arange(df.index.max()+1))
    df = df.reindex(index=index, columns=index, fill_value=0)
    nRM = np.arange(-0.5, df.index.max()+1.5)

    h = axes[ax_id].pcolor(nRM, nRM, df.to_numpy(), cmap='Purples', vmin=0, vmax=0.6)

    if ax_id == len(axes) - 1 :
        cax = inset_axes(axes[-1], width='5%', height='95%', loc='center right')
        cbar = fig.colorbar(h, cax=cax)
        cbar.ax.yaxis.set_ticks_position('left')

    axes[ax_id].set_xticks(np.arange(19, step=3))
    axes[ax_id].set_yticks(np.arange(19, step=3))
    axes[ax_id].set_xlim(-0.5, 18.5)
    axes[ax_id].set_ylim(-0.5, 18.5)
    axes[ax_id].set_aspect('equal', adjustable='box')
    axes[ax_id].set_xlabel('#RM$_i$')
    axes[ax_id].set_ylabel('#RM$_j$')



# Add the figure labels
add_figure_labels(['A', 'B'], axes, dx=-0.062, dy=0.013)

fig.savefig(os.path.join(fig_path, 'figS5.png'), bbox_inches = 'tight')