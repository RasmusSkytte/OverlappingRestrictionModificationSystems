import matplotlib.pyplot as plt

from helpers import make_dirs, get_path, set_rc_params
from helpers import load_reference_data, iter_reference_data
from helpers import compute_network, plot_bipartite_network, plot_network_histogram

import os

set_rc_params()

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_S4')
make_dirs(os.path.join(fig_path))

lpad = 5   # Distance between x/y labels and the axes
A_ij, s, genus_id_to_name = load_reference_data()


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
plt.subplots_adjust(wspace=0.05)

# Loop over subsets to generate figures
for strain_id, A_ij_s, Bs in iter_reference_data(A_ij, s) :

    genus = genus_id_to_name[strain_id]

    # Get the matrices for the networks
    _, A_RMs_s, A_Bs_s = compute_network(Bs)

    # Add the plots
    plot_bipartite_network(Bs, A_ij_s, A_RMs_s, genus, scaling=0.5, ax=ax[strain_id])

# Finalize figures
fig.savefig(os.path.join(fig_path, 'figS4.png'), bbox_inches = 'tight')