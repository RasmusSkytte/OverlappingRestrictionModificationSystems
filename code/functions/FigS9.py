import numpy as np
import matplotlib.pyplot as plt
import os
from p_tqdm import p_umap

from helpers import default_parameters, get_RMs, dynamical_system, make_dirs, get_path, get_solver, set_rc_params, add_figure_labels

set_rc_params()

# Load parameters
C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations = default_parameters()
Beta = 25
n = 250 # number of runs

fig_path = os.path.join(get_path(), 'Figure_S9')
make_dirs(fig_path)

# Define the solver function
def run_param_set(i) :

    # Define random numbe rgenerator
    rng = np.random.default_rng(i)

    # Determine which case to simulate
    B = [[0], [1], [0, 1]]
    P = [[], [1]]

    # Compute number of RM systems
    RMs = get_RMs(B, P)

    # Number of bacteria and phages
    nB = len(B)
    nP = len(P)

    # Starting condiiton
    deb =  Delta / (Eta * (Beta - 1))
    x0 = deb * np.concatenate((1/nB*np.ones(nB), 10/(nB*nP)*np.ones(nB*nP)))  # Everyone starts with a population of delta/(eta*beta)

    # Define omega and cost
    omega_0 = np.power(10, lb + (ub-lb) * rng.uniform(size = len(RMs)))
    cost    = 1 - f * rng.uniform(size = len(RMs))

    # Run dynamics
    _, P_end, _, _, _ = dynamical_system('Extended', B, P, omega_0=omega_0, cost=cost, params=(Alpha, Beta, Eta, Delta, C, T), x0 = x0, solver=get_solver())

    return P_end.sum(axis=1)


if __name__ == '__main__':

    # Get the data
    PPs = np.array(list(p_umap(run_param_set, np.arange(n), num_cpus=4)))

    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(4, 6))

    for ax in axes :
        pos = ax.get_position()
        pos.x0 += 0.2
        ax.set_position(pos)


    bins = np.linspace(0, 1, 11)
    axes[0].hist(PPs[:, 0] / PPs.sum(axis=1), bins=bins, weights=np.ones(n) / n, color='gainsboro')
    axes[1].hist(PPs[:, 1] / PPs.sum(axis=1), bins=bins, weights=np.ones(n) / n, color='gainsboro')

    axes[0].set_ylim(0, 1)
    axes[0].set_xlim(0, 1)
    axes[1].set_xlabel('Fraction of population')

    for ax in axes :
        ax.set_ylabel('pmf.')

    # Add the figure labels
    add_figure_labels(['A'],      [axes[0]], dx=-1.25, dy=0.025)
    add_figure_labels(['B', 'C'], axes,      dx=-0.13, dy=0.025)

    fig.savefig(os.path.join(fig_path, f'figS9.png'), bbox_inches = 'tight')