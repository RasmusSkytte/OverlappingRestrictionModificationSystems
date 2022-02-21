import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from tqdm import tqdm

import os
import re

from helpers import default_parameters, make_dirs, get_path, set_rc_params, pickle_read, analyze_simulation_data, pickle_write, conserve_RM_degree, add_figure_labels

set_rc_params()
plt.rcParams.update({'font.size': 10.0})

# Start with default parameters
C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations = default_parameters()
iterations = 1e4

# Define sensitivity analysis space
Ks = [50, 100, 200, 400, 800]
Ts = np.logspace(1, 4, 4)

data_path = get_path('data')
fig_path  = get_path()
make_dirs(os.path.join(data_path, 'FigS15_S16_S17'))
make_dirs(os.path.join(fig_path,  'Figure_S17'))


fig, axes = plt.subplots(nrows=4, ncols=len(Ks), figsize=(10, 10))
plt.subplots_adjust(wspace=0.25, hspace=0.3)

# Create plotting function
def add_plot(data, ax, c, T_id, jitter, label = None) :
    #add_boxplot(data, ax, c, T_id, jitter)
    add_errorbarplot(data, ax, c, T_id, jitter, label = label)

def add_boxplot(data, ax, c, T_id, jitter, label = None) :
    ax.boxplot(data, patch_artist=True,
    boxprops     = dict(facecolor = c, color = c),
    capprops     = dict(color=c),
    whiskerprops = dict(color=c),
    flierprops   = dict(color=c, markeredgecolor=c, markersize = 1),
    medianprops  = dict(color=c),
    positions = np.arange(len(data)) - jitter + 2*jitter * T_id / (len(Ts)-1),
    widths = jitter/3,
    label = label)

def add_errorbarplot(data, ax, c, T_id, jitter, label = None) :
    m = list(map(np.mean, data))
    s = list(map(np.std,  data))
    x = np.arange(len(data)) - jitter + 2*jitter * T_id / (len(Ts)-1)

    ax.errorbar(x, m, s,
        fmt = '.', markeredgewidth = 1, elinewidth = 1, markersize = 15*jitter, capsize = 10*jitter, color = c, label = label)



# Connect the axes across tests
for ax_set in axes :
    for ax in ax_set[1:] :
        ax_set[0].get_shared_x_axes().join(ax_set[0], ax)
        ax_set[0].get_shared_y_axes().join(ax_set[0], ax)

# Loop over simulations
for K_id, K in tqdm(enumerate(Ks), total = len(Ks)) :
    for T_id, T in tqdm(enumerate(Ts), total = len(Ts), leave=False) :

        # Load the data
        folder =  f'iterations_1e{np.log10(iterations):.0f}_T_1e{np.log10(T):.1f}_C_1e{np.log10(C):.1f}'
        sub_path = os.path.join('FigS15_S16_S17', folder)
        lname = os.path.join(data_path, sub_path, f'RM_{K}_seed_0.pkl')
        if os.path.exists(lname) :
            _, _, _, _, _, _, mRM, _, _, B_samples, _, _, _, _, _, _, _, _, _, _, _, _, _, _, st = pickle_read(lname)
        else :
            raise ValueError(f'{lname} is missing!')

        # Check if data has been generated
        lpath = os.path.join(data_path, sub_path, 'data_simulations.pkl')
        if not os.path.exists(lpath) :

            # Generate and store the data
            pickle_write(lpath, *analyze_simulation_data(n_its = int(np.log10(iterations)), n_seeds = 1, conserve_RM_degree = conserve_RM_degree(), data_path = sub_path))

        # Load the data
        diff_random_unique, diff_random_o_ij, _, _, names, iters = pickle_read(lpath)

        # Filter out the excess data
        keep = [int(re.findall(r'\d+', name)[0]) == K for name in names]
        diff_random_o_ij     = [x for x, keep in zip(diff_random_o_ij,     keep) if keep]
        diff_random_unique   = [x for x, keep in zip(diff_random_unique,   keep) if keep]


        # Add label to the column
        axes[0][K_id].text(2, 3.5, f'{K=}', horizontalalignment = 'left')

        # Get the current color
        c = plt.cm.tab10(T_id)

        # Plot the average RM count
        if T == 1e4 :
            fmt = '--'
        else :
            fmt = '-'

        t = np.arange(iterations)
        I = np.unique(np.round(np.logspace(0, np.log10(t[-1]), 100)-1).round().astype(int))  # Logartihmically sample points to plot -- fixes dotted line issue
        axes[0][K_id].semilogx(t[I], mRM[I], fmt, color= c, label = f'$T = 10^{int(np.log10(T))}$')

        # Plot the abundency distribution
        jitter = 0.25
        add_plot([list(map(len, B_sample)) for B_sample in B_samples],
            axes[1][K_id], c, T_id, jitter, label = f'$T = 10^{int(np.log10(T))}$')

        # Plot the difference in overlap
        add_plot(diff_random_o_ij,
            axes[2][K_id], c, T_id, jitter)
        axes[2][K_id].plot([-10, 10], [0, 0], 'k:', linewidth = 1, zorder = -10)

        # Plot the difference in uniqueness
        add_plot(diff_random_unique,
            axes[3][K_id], c, T_id, jitter)


axes[0][0].set_xlim(1, iterations)
axes[0][0].set_ylim(0, 4)
axes[0][0].set_ylabel('$\langle $#$RM\\rangle$')

for ax_id, ax in enumerate(axes[0]) :
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins = 10))
    ax.set_xticks(np.power(10, np.arange(np.log10(iterations)+1, step=1)))

    if ax_id > 0 :
        ax.set_yticklabels([])




for ax_id, ax in enumerate(axes[1]) :
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins = 5))

axes[1][0].set_ylim(0, 4)
axes[1][0].set_ylim(0, )
axes[1][0].set_ylabel('#RM', labelpad = 4)



axes[2][0].ticklabel_format(style = 'sci', scilimits=(0,0), axis='y')
axes[2][0].yaxis.major.formatter._useMathText = True

axes[2][0].set_ylim(-0.02, 0.08)
axes[2][0].set_ylabel('$\langle I\;/\;U \\rangle$ - $\langle I\;/\;U\\rangle_{rand}$', labelpad = 0)


for ax_id, ax in enumerate(axes[3]) :
    ax.set_xlabel('Time (# Additions)')

axes[3][0].set_ylim(-0.02, 0.82)
axes[3][0].set_ylabel('$f^{u} - f^{u}_{rand}$')


# Set options for the last 3 rows
ticks = np.arange(len(B_samples))
ticklabels = [f'$10^{n:.0f}$' for n in ticks]

for row in axes[1:] :
    for ax_id, ax in enumerate(row) :

        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels)

        ax.set_xlim(2-2*jitter, len(ticks)-1+2*jitter)

        for x in np.linspace(2.5, 3.5, 2) :
            ax.axvline(x, color = 'k', linewidth = 1, alpha = 0.4)

        if ax_id > 0 :
            ax.set_yticklabels([])



# Add the figure labels
add_figure_labels(['A', 'B', 'C'],  axes[:-1, 0], dx=-0.045, dy=0.012)
add_figure_labels(['D'],           [axes[-1, 0]], dx=-0.045, dy=0.009)

axes[0][-1].legend(bbox_to_anchor=(1.05, 1.04), loc='upper left').get_frame().set_edgecolor('k')
axes[1][-1].legend(bbox_to_anchor=(1.05, 1.04), loc='upper left').get_frame().set_edgecolor('k')

fig.savefig(os.path.join(fig_path, 'Figure_S17', 'figS17.png'), bbox_inches = 'tight')