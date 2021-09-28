import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

from tqdm import tqdm

import os

from helpers import make_dirs, get_path, set_rc_params, pickle_read, compute_gamma_and_omega, add_figure_labels

set_rc_params()

iterations = 1e4

data_path = get_path('data')
fig_path  = get_path()
make_dirs(os.path.join(fig_path, 'Figure_S12'))


fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 6))
plt.subplots_adjust(wspace=0.3)

# Loop over simulations
for simulation_id, (xlim, ylim) in tqdm(enumerate(zip([0.97, 0.94], [-10.0, -10.0])), total=2) :

    # Set settings for the run
    if simulation_id == 0 :
        RMs = np.arange(800)

    elif simulation_id == 1 :
        RMs = np.arange(50)

    # Load the data
    lname = os.path.join(data_path, 'FigS12', f'iterations_1e{np.log10(iterations):.1f}', f'RM_{len(RMs)}_seed_0.pkl')
    if os.path.exists(lname) :
        B, B_end, P_end, bacteria, phages, diversity, _, cost, omega_0, _, _, _, B_all, alpha, _, eta, _, C, _, _, _, _, _, _, _ = pickle_read(lname)
    else :
        raise ValueError(f'{lname} is missing!')

    # Define the extinction function
    extinction_function = lambda g, o : g * (1 - bacteria[-1] / C)  - alpha - eta * o * phages[-1] > 0

    # Select strains in the plotting window
    B_possible = [b for b in B_all if np.prod(cost[b]) > xlim]

    # Add strains that did survive despite being above the extinction curve
    B_possible.extend([b for b in B if b not in B_possible])

    # Determine the indicies of the posible strains that went extinct
    I = np.array([b not in B for b in B_possible])
    I = np.argwhere(I).flatten() - np.arange(I.sum())

    # Get the biomass and phage counts for the possible strains
    B_end_possible = np.insert(B_end,       I, np.zeros_like(I))    # extinct strains have zero biomass
    P_end_possible = np.insert(P_end[0, :], I, np.zeros_like(I))    # extinct strains have phage count of zero

    # Determine the naive parameters of the strains
    gammas_naive = np.array([np.prod(cost[b])    for b in B_possible])
    omegas_naive = np.array([np.prod(omega_0[b]) for b in B_possible])

    # Get the omegas matrix based on the possible strains
    _, omega = compute_gamma_and_omega(B_possible, [[]], cost, omega_0, model='Extended')

    # Weight omega by phage counts
    omega_masked = np.ma.array(omega[0], mask=False)                                        # Convert omega to masked array
    ww = np.ma.array(np.repeat(P_end_possible[np.newaxis, :], len(omega_masked), axis=0))   # Create similarily masked array with the phage counts
    omega_masked.mask = ww.mask = np.eye(len(omega_masked)) > 0                             # Mask the diagonal elements of the arrays
    avg_effective_omega = (np.sum(omega_masked * ww, axis=1) / np.sum(ww, axis=1)).data     # Compute the effective omegas

    # Color the nodes based on survival state
    viridis = cm.get_cmap('viridis')
    colors = [viridis(biomass / B_end.max()) for biomass in B_end]

    # Create plots with naive omegas and with the weighted omegas
    for plot_id, (ax, omegas, ylabel, suffix) in enumerate(zip([axes[0][simulation_id], axes[1][simulation_id]], [omegas_naive, avg_effective_omega], ['$\omega$', '$\omega$ (effective)'], ['naiive', 'weighted'])) :

        # Plot the data
        ax.scatter([g for g, b in zip(gammas_naive, B_end_possible) if b != 0.0],
                   [o for o, b in zip(omegas,       B_end_possible) if b != 0.0],
                   s=14, c=colors, marker='o', linewidths=0.5, edgecolors='w')

        ax.scatter([g for g, b in zip(gammas_naive, B_end_possible) if b == 0.0],
                   [o for o, b in zip(omegas,       B_end_possible) if b == 0.0],
                   s=10, c='k',    marker='x', linewidths=0.5)

        # Plot the extinciton line
        oo = np.logspace(ylim, 0, 200)
        gg = (alpha + oo * eta * phages[-1]) / (1 - bacteria[-1] / C)
        ax.plot(gg, oo, 'r:')

        ax.set_xlim(xlim, 1)

        ax.set_ylim(np.power(10, ylim), 1)
        ax.set_yscale('log')

        if simulation_id == 0 :
            ax.set_ylabel(ylabel)

        if plot_id == 1 :
            ax.set_xlabel('$\gamma$')


# Add the figure labels
add_figure_labels(['A', 'B', 'C', 'D'], axes.flatten(), dx=-0.08, dy=0.024)


fig.savefig(os.path.join(fig_path, 'Figure_S12', 'figS12.png'), bbox_inches = 'tight')