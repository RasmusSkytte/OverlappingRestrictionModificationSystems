import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import os
from tqdm import tqdm, trange

from helpers import default_parameters, get_RMs, dynamical_system, make_dirs
from helpers import get_path, get_solver, set_rc_params, add_figure_labels, pickle_read

set_rc_params()

# Load parameters
C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations = default_parameters()
Beta = 25

fig_path = os.path.join(get_path(), 'Figure_5')
make_dirs(fig_path)

fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(8, 10))
axes = axes.flatten()
for ax in axes[np.r_[:2, 3:5, 6:8]] :
    ax.remove()
axes = axes[np.r_[2, 5, 8:12]]

plt.subplots_adjust(hspace=0.5, wspace=0.35)

for ax in axes[:3] :
    pos = ax.get_position()
    pos.x0 = 0.5
    ax.set_position(pos)


# Loop over cases to plot
for ax_id in trange(3) :

    # Set seed
    np.random.seed(8)

    # Determine which case to simulate
    if ax_id == 0 :
        B = [[0], [1]]
        P = [[]]
    elif ax_id == 1 :
        B = [[0], [1], [0, 1]]
        P = [[]]
    elif ax_id == 2:
        B = [[0], [1], [0, 1]]
        P = [[1]]

    # Compute number of RM systems
    RMs = get_RMs(B, P)

    # Number of bacteria and phages
    nB = len(B)
    nP = len(P)

    # Starting condiiton
    deb =  Delta / (Eta * (Beta - 1))
    x0 = deb * np.concatenate((1/nB*np.ones(nB), 10/(nB*nP)*np.ones(nB*nP)))  # Everyone starts with a population of delta/(eta*beta)

    # Define omega and cost
    omega_0 = np.power(10, lb + (ub-lb) * np.random.uniform(size = len(RMs)))
    cost    = 1 - f * np.random.uniform(size = len(RMs))

    # Run dynamics
    _, _, y, t, _ = dynamical_system('Extended', B, P, omega_0=omega_0, cost=cost, params=(Alpha, Beta, Eta, Delta, C, T), x0 = x0, solver=get_solver())

    t[0] = np.finfo(float).eps

    c = ['dodgerblue', 'gold']
    for i in range(nB) :
        axes[ax_id].plot(t, y[i, :], color='k', linewidth=4)
        if i < 2 :
            axes[ax_id].plot(t, y[i, :],       color=c[i], linewidth=2.5)
        else :
            axes[ax_id].plot(t, y[i, :],       color=c[0], linewidth=2.5)
            axes[ax_id].plot(t, y[i, :], '--', color=c[1], linewidth=2.5)

    axes[ax_id].set_xlim(1, t[-1])
    axes[ax_id].set_ylim(1, C)
    axes[ax_id].set_xscale('log')
    axes[ax_id].set_yscale('log')
    axes[ax_id].set_xlabel('Time (generations)')
    axes[ax_id].set_ylabel('Population Size')

# Add the figure labels
add_figure_labels(['A', 'B', 'C'], axes, dx=-0.43, dy=0.025)












# Load parameters
C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations = default_parameters()
Beta = 25

# Number of repititions
N = 17

# Set up test space
omega, gamma = np.meshgrid(np.power(10, np.linspace(0, -2, N)), np.linspace(0, 1, N));

# Get the path
data_path = get_path('data')


axes = axes[3:]
for ax in axes :
    pos = ax.get_position()
    pos.y0 -= 0.03
    ax.set_position(pos)
    ax.set_xlabel('$\gamma$')

axes[0].set_ylabel('$\\omega$')


# Prepare data arrays
ps         = np.zeros(np.shape(omega))
biomass    = np.zeros(np.shape(omega))
biomass_AB = np.zeros(np.shape(omega))

# Loop over omega, gamma
for i, (o, g) in tqdm(enumerate(zip(omega.flatten(), gamma.flatten())), total=N**2) :
    ind = np.unravel_index(i, np.shape(omega))

    # Load data
    fname = os.path.join(data_path, f'Fig5')
    labels = ['A', 'B', 'AB']

    fname  = os.path.join(fname, f'omega_1e{np.log10(o):.3f}_gamma_{g:.4f}.pkl')


    if os.path.exists(fname) :
        coordinates, BBs, PPs = pickle_read(fname)
    else :
        raise ValueError(f'Data not loaded! ({fname}))')

    # Set booleans
    b_A      = False
    b_A_B    = False
    b_AB     = False
    b_B_AB   = False
    b_A_B_AB = False

    # Set up starting coordinates
    deb =      Delta / (Eta * (Beta - 1))
    deb2 = 1.9*Delta / (Eta * (Beta * (1 + o) - 2))  # Use 1.9 instead of 2 to break degeneracy when omega  = 1
    B0 = np.array([deb, deb2, 1e1, 1e3, 1e5, 1e7])

    b0_1, b0_2 = np.meshgrid(B0, B0)

    # Loop over starting conditions
    for j, (b1, b2) in enumerate(zip(b0_1.flatten(), b0_2.flatten())) :

        # Retrieve B_end
        I = np.logical_and(coordinates[:, 0]==b1, coordinates[:, 1]==b2)
        B_end = BBs[I][0]
        biomass[ind]    += np.sum(B_end) / np.size(b0_1)
        biomass_AB[ind] += B_end[-1]     / np.size(b0_1)

        # Determine posible soluitions
        thres = 1
        outcome = tuple((B_end < thres).astype(int))

        if outcome == (1, 0, 0) :
            if not b_A :
                b_A = True
                ps[ind] += 2 ** 0

        elif outcome == (1, 1, 0) :
            if not b_A_B :
                b_A_B = True
                ps[ind] += 2 ** 1

        elif outcome == (1, 1, 1) :
            if not b_A_B_AB :
                b_A_B_AB = True
                ps[ind] += 2 ** 2

        elif outcome == (0, 1, 1) :
            if not b_B_AB :
                b_B_AB = True
                ps[ind] += 2 ** 3

        elif outcome == (0, 0, 1) :
            if not b_AB :
                b_AB = True
                ps[ind] += 2 ** 4

        elif not outcome == (0, 0, 0) : # other non trivial solution exists
            raise ValueError(f'Outcome : {outcome} not defined')

# Plot solution space
dx = 1/(2 * (len(gamma)-1))
dy = 2/(2 * (len(omega)-1))
axes[0].imshow(ps.T, extent=[-dx, 1+dx, -2-dy, dy], aspect=0.5, cmap='tab10')
yticks = np.linspace(0, -2, 3)

# Add annotation
fontsize = 10
axes[0].text(0.10, -1.0, 'No Solution',                          color = 'w', fontweight ='bold', horizontalalignment = 'center', verticalalignment = 'center', fontsize = fontsize,  rotation = 90)
axes[0].text(0.40, -0.3, f'{labels[0]}/{labels[1]}',             color = 'w', fontweight ='bold', horizontalalignment = 'center', verticalalignment = 'center', fontsize = fontsize)
axes[0].text(0.75, -1.5, f'{labels[0]}/{labels[1]}/{labels[2]}', color = 'w', fontweight ='bold', horizontalalignment = 'center', verticalalignment = 'center', fontsize = fontsize)
axes[0].text(1.00, -1.0,    labels[2],                           color = 'w', fontweight ='bold', horizontalalignment = 'center', verticalalignment = 'center', fontsize = fontsize,  rotation = 90)

# Plot change in biomass
B_ref = 2 * Delta / (Eta * (Beta * (1 + omega) - 2))
DeltaB = (biomass - B_ref) / B_ref
h1 = axes[1].imshow(DeltaB[4:, :].T, extent=[gamma[4, 0]-dx, 1+dx, -2-dy, dy], aspect=0.5, cmap='coolwarm', vmin=-1.0, vmax=1.0)
axes[1].set_yticks(yticks)

# Plot relative biomass of AB solution
RelB = biomass_AB / biomass
h2 = axes[2].imshow(RelB[4:, :].T,   extent=[gamma[4, 0]-dx, 1+dx, -2-dy, dy], aspect=0.5, cmap='coolwarm', vmin= -1.0, vmax=1.0)

for ax in axes :
    ax.set_facecolor((0.6, 0.6, 0.6))
    ax.set_xlim(-dx, 1+dx)
    ax.set_ylim(-2-dy, dy)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f'$10^{{{y:.0f}}}$' for y in yticks])


for ax, h in zip(axes[1:], [h1, h2]) :
    axins = inset_axes(ax, width='2.5%', height='90%', loc='center left', bbox_to_anchor=(-0.02, 0, 1, 1), bbox_transform=ax.transAxes)
    vmin, vmax = h.get_clim()
    cb = fig.colorbar(h, ax=ax, cax=axins, ticks=np.linspace(vmin, vmax, 3))
    cb.ax.tick_params(labelsize=7)
    cb.ax.yaxis.set_tick_params(color='w')
    cb.outline.set_edgecolor('w')
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='w')


fig.patch.set_facecolor('w')


# Guide lines
k = lambda o : 1 - 2 * Delta / (C * Eta * (Beta * (1 + o) - 2))
D = lambda o : o**2 * k(o)**2 + k(o)**2 * Alpha * (1 + o) * (1 - o)
f = lambda o : (o * k(o) + np.sqrt( D(o) )) / ((1 + o) * k(o))
f = np.vectorize(f)
ww = np.logspace(-2.1, 0.1)
gg = f(ww)


D = lambda o : o**4 * k(o)**2 + k(o)**2 * Alpha * (1 + o) * (1 + o - 2 * o**2)
f = lambda o : (o**2 * k(o) + np.sqrt( D(o) )) / ((1 + o) * k(o))
f = np.vectorize(f)
gg_n = f(ww)

for ax in axes :
    ax.plot([Alpha, Alpha], ax.get_ylim(), 'w-.', label='$\gamma=\\alpha$')
    ax.plot(gg,             np.log10(ww),  'w-',  label='$\gamma=f(\omega)$')


# Add the figure labels
add_figure_labels(['D', 'E', 'F'], axes, dx=-0.063, dy=0.01)


fig.savefig(os.path.join(fig_path, f'fig5.png'), bbox_inches = 'tight')