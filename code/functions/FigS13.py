import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import scipy.io

import os

from helpers import make_dirs, get_path, set_rc_params, add_figure_labels
from helpers import iter_simulation_data, compute_network


set_rc_params()

n_open_triplets = []
n_closed_triplets = []
Ks = []
iterations = []

# Loop over subsets to generate figures
for K, seed_id, log_iterations in iter_simulation_data() :

    # Skip early data points
    if log_iterations < 3 :
        continue

    # Load the data
    lname = os.path.join(get_path('data'), 'Fig6', f'RM_{K}_seed_{seed_id}.mat')
    if os.path.exists(lname) :
        data = scipy.io.loadmat(lname)

    else :
        raise ValueError(f'{lname} is missing!')


    # Get data from the run
    Bs = [b[0, :].tolist() for b in data['B_samples'][0, :][log_iterations].flatten().tolist()] # wtf scipy?

    A_ij, A_RMs, A_Bs = compute_network(Bs)

    A = (A_Bs > 0).astype(int) - np.eye(len(A_Bs))

    A2 = np.dot(A, A)
    A3 = np.dot(A2, A)

    tA2 = np.trace(A2)
    tA3 = np.trace(A3)

    closed_triplets = tA3 / (2*3) # 2 paths, 3 members in the triplet
    open_triplets = (A2.sum() - tA2 - tA3) / 2

    n_open_triplets.append(open_triplets)
    n_closed_triplets.append(closed_triplets)
    Ks.append(K)
    iterations.append(log_iterations)


# Combine results
df = pd.DataFrame({'K' : Ks, 'iterations' : iterations, 'closed_triplets' : n_closed_triplets, 'open_triplets' : n_open_triplets})
df['triplet_ratio'] = df.closed_triplets / (df.closed_triplets + df.open_triplets)

mean = df.groupby(['K', 'iterations']).mean()
std =  df.groupby(['K', 'iterations']).std()


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(8, 3), sharex=True)
plt.subplots_adjust(wspace = 0.45)

for K in np.unique(Ks) :
    n = np.sum(df.K == K)
    m = mean.loc[pd.IndexSlice[K, :]]
    s = std.loc[pd.IndexSlice[K, :]]

    axes[0].errorbar(np.power(10, m.index), m.open_triplets,   s.open_triplets,   fmt = '.', capsize = 5, markeredgewidth = 2)
    axes[1].errorbar(np.power(10, m.index), m.closed_triplets, s.closed_triplets, fmt = '.', capsize = 5, markeredgewidth = 2)
    axes[2].errorbar(np.power(10, m.index), m.triplet_ratio,   s.triplet_ratio,   fmt = '.', capsize = 5, markeredgewidth = 2, label=f'K = {K}')


for ax in axes :
    ax.set_xscale('log')
    ax.set_xlabel('Time (# Additions)')
    ax.set_xticks(np.power(10, np.arange(3, 7)))

for ax in axes[:2] :
    ax.set_yscale('log')
    ax.set_ylim(1e-1, 1e7)


axes[0].set_ylabel('# hierarchical triplets', labelpad = -4)
axes[1].set_ylabel('# looped triplets', labelpad = -4)

axes[2].set_ylim(0, 1)
axes[2].set_ylabel('Fraction looped triplets')
axes[2].legend(bbox_to_anchor=(1.05, 1.04), loc='upper left').get_frame().set_edgecolor('k')


# Add the figure labels
add_figure_labels(['A', 'B', 'C'], axes, dx=-0.067, dy=0.05)


# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_S13')
make_dirs(fig_path)

fig.savefig(os.path.join(fig_path, 'figS13.png'), bbox_inches = 'tight')