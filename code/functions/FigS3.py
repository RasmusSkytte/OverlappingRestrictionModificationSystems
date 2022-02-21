import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from helpers import make_dirs, get_path, set_rc_params
from helpers import load_reference_data, iter_reference_data, compute_relative_overlaps
from helpers import add_figure_labels

import os

set_rc_params()

# Prepare figure folders
fig_path = get_path()
make_dirs(os.path.join(fig_path, 'Figure_S3'))

A_ij, s, genus_id_to_name = load_reference_data()

# Prepare figure
lpad = 5   # Distance between x/y labels and the axes
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(16, 6))
fig.subplots_adjust(wspace=0.3)
axes = axes.flatten()


n_max = 18
# Loop over subsets to generate figures
for strain_id, A_ij_s, Bs in iter_reference_data(A_ij, s) :

    # Get the distribution of number of RM systems
    data = np.array(np.sum(A_ij_s > 0, axis=1))
    axes[0].hist(data, bins=np.arange(start=-0.5, stop=np.max(data)+1.5), color=plt.cm.tab10(strain_id), linewidth=plt.rcParams['axes.linewidth'], edgecolor='k', density=True, alpha = 0.5)

    data = np.array(np.sum(A_ij_s > 0, axis=0))
    data = np.clip(data, 0, n_max)
    data = data[data > 0]
    axes[1].hist(data, bins=np.arange(start=-0.5, stop=np.max(data)+1.5), color=plt.cm.tab10(strain_id), linewidth=plt.rcParams['axes.linewidth'], edgecolor='k', density=True, alpha = 0.5, label=genus_id_to_name[strain_id])

    # Get the overlap for the subset
    data = compute_relative_overlaps(A_ij_s, Bs)
    weights = None
    weights = np.ones_like(data)/len(data)
    axes[2].hist(data, bins=np.linspace(0, 1, 11), weights=weights,       color=plt.cm.tab10(strain_id), linewidth=plt.rcParams['axes.linewidth'], edgecolor='k',               alpha = 0.5)


# Plot the RM per B distribution
axes[0].set_xticks(np.arange(0, n_max+2, step=2))
axes[0].set_xlim(-0.5, n_max+0.5)
axes[0].set_ylim(0, 0.4)
axes[0].set_xlabel('# RM per B', labelpad=lpad)
axes[0].set_ylabel('pmf.',       labelpad=lpad)


# Plot the B per RM distribution
axes[1].set_xticks(np.arange(0, n_max+2, step=2))
axes[1].set_xlim(-0.5, n_max+0.5)
axes[1].set_ylim(0, 0.5)
ticks_loc = axes[1].get_xticks().tolist()
ticks_labels = ['{:,.0f}'.format(x) for x in ticks_loc]
ticks_labels[-1] = ' ' + ticks_labels[-1] + '+'
axes[1].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes[1].xaxis.set_ticklabels(ticks_labels)
axes[1].set_xlabel('# B per RM', labelpad=lpad)
axes[1].legend().get_frame().set_edgecolor('k')


# Plot the histogram of overlap
axes[2].set_xlim(0, 1)
axes[2].set_ylim(0, 0.35)
axes[2].set_xlabel('# Intersection / # Union', labelpad=lpad)


# Add the figure labels
add_figure_labels(['A', 'B', 'C'], axes, dx=-0.04, dy=0.025)

fig.savefig(os.path.join(fig_path, 'Figure_S3', 'figS3.png'), bbox_inches = 'tight')