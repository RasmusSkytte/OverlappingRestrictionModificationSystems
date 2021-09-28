import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from helpers import make_dirs, flatten_list, get_path, set_rc_params, get_metrics_of_real_network
from helpers import load_sequence_data, iter_sequence_data, compute_relative_overlaps
from helpers import add_figure_labels

import os
import inspect

set_rc_params()

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_2')
make_dirs(fig_path)


A_ij, s, data_taxonomy, data_motifs, genus_id_to_name = load_sequence_data()

# Allocate arrays
nRM_per_B = []  # B : bacterium
nB_per_RM = []

df_overlap = pd.DataFrame({'genus' : [], 'genus_id' : [], 'o_ij' : []})

max_o_ij = []
nSamples = 0
nFiltered = 0
nViableSamples = 0
nUnviable = len(np.unique(s))
nViable = 0
names = []
found_RMs = []
found_motifs = np.zeros(np.size(A_ij, 1))
found_genus_ids = []

hist_strain_counts = []
hist_avg_RM_per_B  = []


# Count the total filtered sample
for _, A_ij_s, _, _ in iter_sequence_data(A_ij, s, data_taxonomy, treshold=0, tqdm_disable=True) :
    nSamples += len(A_ij_s)

# Loop over subsets to generate figures
for genus_id, A_ij_s, Bs, genus in iter_sequence_data(A_ij, s, data_taxonomy) :

    # Count the samples
    nViableSamples += np.size(A_ij_s, 0)

    # Note the ID of the found RM systems
    found_RMs.append(np.where(A_ij_s.sum(axis=0))[0].tolist())
    found_motifs += A_ij_s.sum(axis=0)
    found_genus_ids.append(genus_id)

    # Count the viable genera
    nViable += 1
    nUnviable -= 1

    # Get the overlap for the subset
    o_ij_s = compute_relative_overlaps(A_ij_s, Bs, ind='all')

    # Get the distribution of number of RM systems
    nB_per_RM.append(np.sum(A_ij_s > 0, axis=0))
    nRM_per_B.append(np.sum(A_ij_s > 0, axis=1))


    # Get the metrics from the simulation
    _, RMs_per_B, _, _, _, _ = get_metrics_of_real_network(Bs)
    hist_avg_RM_per_B.append(RMs_per_B.mean())

    max_o_ij.append(np.nanmax(o_ij_s, axis=0))  # Take the marginal nansum to get the max overlap

    # Exclude diagnoal elements
    triu_ind_x, triu_ind_y = np.triu_indices(len(Bs), k=1)
    tril_ind_x, tril_ind_y = np.tril_indices(len(Bs), k=-1)
    ind = (np.concatenate((triu_ind_x, tril_ind_x)), np.concatenate((triu_ind_y, tril_ind_y)))

    df_overlap_s = pd.DataFrame({
                        'genus' : genus,
                        'genus_id' : genus_id,
                        'o_ij' : o_ij_s[ind]})

    df_overlap = df_overlap.append(df_overlap_s)

    # Store average measures
    hist_strain_counts.append(len(Bs))

treshold = inspect.signature(iter_sequence_data).parameters['treshold'].default
print('-------------- Data report -----------------')
print(f'Initial sample contains : {A_ij.shape[0]} strains with {A_ij.shape[1]} RM identifiers')
print(f'Filtering out duplicates within genus yields {nSamples} strains')
print(f'After filtering {nViable} genera contains {treshold} or more samples and {nUnviable} contains less')
print(f'In total we have {nViable} genera with a total of {nViableSamples} samples')
print(f'Out of {A_ij.shape[1]} RM identifiers searched for, {len(np.unique(flatten_list(found_RMs)))} RM identifiers were found in the final sample')

# Store the bacterial strain IDs
data_taxonomy.index[np.isin(data_taxonomy['genus_id'], np.array(found_genus_ids))].to_frame().to_csv(os.path.join(fig_path,  'found_strains.csv'), index = False)

# Store list of found motifs
I = found_motifs > 0
data_motifs.loc[I].sort_values(by='R-M motif').to_csv(os.path.join(fig_path,  'found_motifs.csv'),   index = False)
data_motifs.loc[~I].sort_values(by='R-M motif').to_csv(os.path.join(fig_path, 'missing_motifs.csv'), index = False)



lpad = 5   # Distance between x/y labels and the axes
dy = '   '
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(8, 10))
fig.subplots_adjust(wspace=0.25, hspace=0.3)
axes = axes.flatten()



# Average RM count per genera
data = np.array(list(map(np.mean, nRM_per_B)))
weights = np.ones_like(data)
axes[0].hist(data, bins=np.arange(start=-0.5, stop=np.max(data)+1.5), weights=weights, color=plt.cm.Pastel1(0), linewidth=plt.rcParams['axes.linewidth'], edgecolor='k', density=True)
axes[0].set_xticks(np.arange(0, 20, step=2))
axes[0].set_xlim(-0.5, 18+0.5)
axes[0].set_ylim(0, 0.5)
axes[0].set_xlabel('$\langle $# RM systems$\\rangle$ per genera', labelpad=lpad)
axes[0].set_ylabel('Fraction of genera'+dy,                       labelpad=lpad)


# Number of genera RM system is found in
RMs = np.arange(A_ij.shape[1])
nGenera_per_RM = np.array([np.isin(RMs, genus) for genus in found_RMs]).sum(axis=0)

# Distribution of RM across genera
data = np.clip(nGenera_per_RM, 0, 18)
data = data[data > 0]
weights = np.ones_like(data)
axes[1].hist(data, bins=np.arange(start=-0.5, stop=np.max(data)+1.5), weights=weights, color=plt.cm.Pastel1(1), linewidth=plt.rcParams['axes.linewidth'], edgecolor='k', density=True)
axes[1].set_xticks(np.arange(0, 20, step=2))
axes[1].set_xlim(-0.5, np.max(data)+0.5)
ticks_loc = axes[1].get_xticks().tolist()
ticks_labels = ['{:,.0f}'.format(x) for x in ticks_loc]
ticks_labels[-1] = ' ' + ticks_labels[-1] + '+'
axes[1].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes[1].xaxis.set_ticklabels(ticks_labels)
axes[1].set_ylim(0, 0.4)
axes[1].set_xlabel('# genera where RM system is found', labelpad=lpad)
axes[1].set_ylabel('Fraction of RM systems'+dy,         labelpad=lpad)


# Plot the RM per B distribution
data = np.array(flatten_list(nRM_per_B))
print(f'Average RM abundency across genera is {data.mean():.2f}, with a maximum abundecy of {data.max()}')
weights = flatten_list([np.full(len(n), fill_value = 1/len(n)) for n in nRM_per_B])
axes[2].hist(data, bins=np.arange(start=-0.5, stop=np.max(data)+1.5), weights=weights, color=plt.cm.Pastel1(2), linewidth=plt.rcParams['axes.linewidth'], edgecolor='k', density=True)
axes[2].set_xticks(np.arange(0, 20, step=2))
axes[2].set_xlim(-0.5, np.max(data)+0.5)
axes[2].set_ylim(0, 0.5)
axes[2].set_xlabel('# RM systems per bacterium', labelpad=lpad)
axes[2].set_ylabel('Fraction of bacteria'+dy,    labelpad=lpad)

# Plot the B per RM distribution
data = np.clip(flatten_list(nB_per_RM), 0, 18)
weights = np.array(flatten_list([np.full(len(n), fill_value = 1/len(n)) for n in nB_per_RM]))
weights = weights[data > 0]
data    = data[data > 0]
axes[3].hist(data, bins=np.arange(start=0.5, stop=np.max(data)+1.5), weights=weights, color=plt.cm.Pastel1(3), linewidth=plt.rcParams['axes.linewidth'], edgecolor='k', density=True)
axes[3].set_xticks(np.arange(0, 20, step=2))
axes[3].set_xlim(-0.5, np.max(data)+0.5)
ticks_loc = axes[3].get_xticks().tolist()
ticks_labels = ['{:,.0f}'.format(x) for x in ticks_loc]
ticks_labels[-1] = ' ' + ticks_labels[-1] + '+'
axes[3].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
axes[3].xaxis.set_ticklabels(ticks_labels)
axes[3].set_ylim(0, 0.6)
axes[3].set_xlabel('# bacteria with RM system', labelpad=lpad)
axes[3].set_ylabel('Fraction of RM systems'+dy, labelpad=lpad)


# Plot the histogram of overlap
data = df_overlap.o_ij
axes[4].hist(data, bins=np.linspace(0, 1, 11), weights=np.ones_like(data)/len(data), color=plt.cm.Pastel1(4))
axes[4].set_xlim(0, 1)
axes[4].set_ylim(0, 0.8)
axes[4].set_xlabel('# intersection / # union',    labelpad=lpad)
axes[4].set_ylabel('Fraction of strain pairs'+dy, labelpad=lpad)

# Plot the histogram of genera samples
data = hist_strain_counts
axes[5].hist(data, bins=np.arange(15, 110, 10), color=plt.cm.Pastel1(5))
axes[5].set_xticks(np.arange(15, 110, step=15))
axes[5].set_xlim(15, 105)
axes[5].set_ylim(0,  25)
axes[5].set_xlabel('# sampled strains'+dy, labelpad=lpad)
axes[5].set_ylabel('# genera',             labelpad=lpad)


# Add the figure labels
add_figure_labels(['A', 'B', 'C', 'D', 'E', 'F'], axes, dx=-0.067, dy=0.015)

fig.savefig(os.path.join(fig_path, 'fig2.png'), bbox_inches = 'tight')