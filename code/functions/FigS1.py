import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

from helpers import make_dirs, get_path, set_rc_params, add_figure_labels, get_metrics_of_real_network
from helpers import load_sequence_data, iter_sequence_data, compute_network
import os

set_rc_params()

n_bootstraps = 100
n_subsets = 9

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_S1')
make_dirs(fig_path)


# Prepare figure
fig, axes = plt.subplots(nrows=1, ncols=3, sharex=True, figsize=(8, 3))
plt.subplots_adjust(wspace=0.5)
axes[0].set_xlim(0, 100)

A_ij, s, data_taxonomy, _, genus_id_to_name = load_sequence_data()

# Determine the genera to compare
average_RM_abundence = np.array([])
for strain_id, A_ij_s, Bs, genus in iter_sequence_data(A_ij, s, data_taxonomy, tqdm_disable=True) :
    average_RM_abundence = np.append(average_RM_abundence, np.mean(A_ij_s.sum(axis=1)))

# Choose examples based on quantiles of average RM abundence
sort_ind = np.argsort(average_RM_abundence)
q_cuts = np.quantile(average_RM_abundence[sort_ind], (0.25, 0.50, 0.75))
q_inds = np.searchsorted(average_RM_abundence[sort_ind], q_cuts)


percentages = (np.arange(n_subsets) + 1) / n_subsets
np.random.seed(0)

# Loop over subsets to generate figures
data_id = 0
with tqdm(total=len(q_inds)) as pbar:
    for d, (strain_id, A_ij_s, Bs, genus) in enumerate(iter_sequence_data(A_ij, s, data_taxonomy, tqdm_disable=True)) :

        # Only plot if index is in q_inds
        if d not in sort_ind[q_inds] :
            continue
        pbar.update(1)

        # Allocate array for bootstrapping
        nRMs     = np.full((n_subsets, n_bootstraps), fill_value=np.nan)
        overlaps = np.full((n_subsets, n_bootstraps), fill_value=np.nan)
        uniques  = np.full((n_subsets, n_bootstraps), fill_value=np.nan)

        n_samples = len(A_ij_s)

        # Draw n_subsets subsets of increasing size
        for percentage_id, percentage in tqdm(enumerate(percentages), total=n_subsets, leave=False) :

            for bootstrap_id in range(n_bootstraps) :

                # Draw bootstrap subset
                inds = np.random.choice(np.arange(n_samples), size=int(np.round(n_samples * percentage)), replace=False)
                A_ij_ss = A_ij_s[inds]
                Bss     = [b for i, b in enumerate(Bs) if i in inds]

                # Compute the network of RM systems
                _, A_RMs_ss, A_Bs_ss = compute_network(Bss)

                # Skip last (only 1 way to sample 100 percent)
                if percentage_id == n_subsets and bootstrap_id > 0 :
                    continue

                # Get the metrics from the subset data
                _, RMs_per_B, _, _, f_unique_ss, avg_o_ij_ss = get_metrics_of_real_network(Bss)

                nRMs[percentage_id, bootstrap_id]     = RMs_per_B.mean()
                overlaps[percentage_id, bootstrap_id] = avg_o_ij_ss
                uniques[percentage_id, bootstrap_id]  = f_unique_ss


        color = plt.cm.tab10(data_id)
        data_id += 1

        axes[0].errorbar(100*percentages[:-1], nRMs[:-1].mean(axis=1),     yerr=nRMs[:-1].std(axis=1)     / np.sqrt(n_bootstraps), fmt = '.', capsize = 5, markeredgewidth = 2, color = color, label=genus)
        axes[1].errorbar(100*percentages[:-1], overlaps[:-1].mean(axis=1), yerr=overlaps[:-1].std(axis=1) / np.sqrt(n_bootstraps), fmt = '.', capsize = 5, markeredgewidth = 2, color = color, label=genus)
        axes[2].errorbar(100*percentages[:-1], uniques[:-1].mean(axis=1),  yerr=overlaps[:-1].std(axis=1) / np.sqrt(n_bootstraps), fmt = '.', capsize = 5, markeredgewidth = 2, color = color, label=genus)

        axes[0].plot([0, 100], [nRMs[-1].mean()]     * 2, ':', color = color)
        axes[1].plot([0, 100], [overlaps[-1].mean()] * 2, ':', color = color)
        axes[2].plot([0, 100], [uniques[-1].mean()]  * 2, ':', color = color)



axes[0].set_ylim(1.2, 2.4)
axes[0].set_xlabel('Percent sampled')
axes[0].set_ylabel('$\langle $#$RM\\rangle$')

axes[1].set_ylim(0, 0.1)
axes[1].set_xlabel('Percent sampled')
axes[1].set_ylabel('$\langle I\;/\;U \\rangle$')

axes[2].set_ylim(0, 1)
axes[2].set_xlabel('Percent sampled')
axes[2].set_ylabel('$f^{ u}$')


axes[2].legend(bbox_to_anchor=(1.05, 1.03), loc='upper left').get_frame().set_edgecolor('k')

# Add the figure labels
add_figure_labels(['A'], [axes[0]], dx=-0.07, dy=0.05)
add_figure_labels(['B'], [axes[1]], dx=-0.08, dy=0.05)
add_figure_labels(['C'], [axes[2]], dx=-0.07, dy=0.05)

fig.savefig(os.path.join(fig_path, 'figS1.png'), bbox_inches = 'tight')