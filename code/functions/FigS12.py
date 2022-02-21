import numpy as np

import matplotlib.pyplot as plt
import scipy.io

from tqdm import tqdm

import os

from helpers import default_parameters, make_dirs, get_path, set_rc_params


set_rc_params()
plt.rcParams.update({'font.size': 10.0}) # These figures need a little smaller font

Ks = [-1, 800, 400, 200, 100, 50]
seeds = [0, 1, 2, 3, 4, 5]

# Load parameters
C, _, _, beta, _, _, _, _, _, _, iterations = default_parameters()

data_path = os.path.join(get_path('data'))


# Prepare figure
x0 = 0.18   # X coord of first col
w0 = 0.18   # Width of first col
dx = 0.05   # Distance between axes


# Loop over seeds to generate panel
for fig_id, K in tqdm(enumerate(Ks), total = len(Ks)) :

    fig = plt.figure(figsize=(8, 4))
    ax_dynamics     = []
    ax_abundency    = []
    ax_distribution = []

    x0 = 0
    for _ in range(len(seeds)) :
        ax_dynamics.append(    fig.add_axes([x0, 0.375, w0, 1.0 - 0.5]))
        ax_abundency.append(   fig.add_axes([x0, 0.23,  w0, 0.12], sharex=ax_dynamics[-1]))
        ax_distribution.append(fig.add_axes([x0, 0.02,  w0, 0.10]))

        x0 += w0 + dx



    ax_annotate = fig.add_subplot(111, facecolor='none')

    xa = -0.23  # X coord of the annotaiton
    ax_annotate.text(xa,  0.54, 'P / (10*C)', verticalalignment = 'center', horizontalalignment = 'left', rotation = 90, color='r')
    ax_annotate.text(xa,  0.72, ', B / C,',   verticalalignment = 'center', horizontalalignment = 'left', rotation = 90, color='k')
    ax_annotate.text(xa,  0.86, 'D / β',      verticalalignment = 'center', horizontalalignment = 'left', rotation = 90, color='b')

    ax_annotate.text(xa,  0.23, '$\langle $#$RM\\rangle$', verticalalignment = 'center', horizontalalignment = 'left', rotation = 90)

    ax_annotate.text(xa, -0.05, 'pmf.', verticalalignment = 'center', horizontalalignment = 'left', rotation = 90)
    ax_annotate.set_axis_off()


    for ax_id, seed in tqdm(enumerate(seeds), total = len(seeds), leave = False) :

        # Add the simulation dynamics to the plot
        if K > 0 :
            RMs = np.arange(K)

            lname = os.path.join(data_path, 'Fig6', f'RM_{K}_seed_{seed}.mat')
            label = f'${K=}$'

            if os.path.exists(lname) :
                data = scipy.io.loadmat(lname)
                data['nRM'] = [len(lst[0, :].tolist()) for lst in data['B_samples'][0, :][-1].flatten().tolist()] # wtf scipy?

            else :
                raise ValueError(f'{lname} is missing!')

        else :
            lname = os.path.join(data_path, 'Fig6', f'RM_inf_seed_{seed}.mat')
            label = '$K =\infty$'

            if os.path.exists(lname) :
                data = scipy.io.loadmat(lname)
                data['nRM'] = data['nRM'].flatten()
            else :
                raise ValueError(f'{lname} is missing!')

        # Plot dynamics
        ax_dynamics[ax_id].plot(np.arange(iterations), data['phages'].flatten()    / (10 * C), 'r', label='P / (10 C)')
        ax_dynamics[ax_id].plot(np.arange(iterations), data['bacteria'].flatten()  / C,        'k', label='C')
        ax_dynamics[ax_id].plot(np.arange(iterations), data['diversity'].flatten() / beta,     'b', label='D / $\\beta$')

        # Label first column only
        if ax_id == 0 :
            ax_dynamics[ax_id].text(6e5, 2.7, label, horizontalalignment = 'right')

        # Plot the average RM count
        ax_abundency[ax_id].plot(np.arange(iterations), data['mRM'].flatten(), color=plt.cm.tab10(4))

        # Plot the histogram
        ax_distribution[ax_id].hist(data['nRM'], bins=np.arange(0.5, 6.5), color=plt.cm.Pastel1(2), density=True)


        # Dynamics axes
        ax_dynamics[ax_id].set_yticks(range(4))

        ax_dynamics[ax_id].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
        ax_dynamics[ax_id].set_xlim(1, iterations)
        ax_dynamics[ax_id].set_xticks(np.power(10, np.arange(np.log10(iterations)+1)))
        ax_dynamics[ax_id].set_xscale('log')

        # Average RM axes
        ax_abundency[ax_id].set_yticks(range(0, 5, 2))
        ax_abundency[ax_id].set_ylim(0, 4.9)

        ax_abundency[ax_id].set_xlim(1, iterations)
        ax_abundency[ax_id].set_xscale('log')
        ax_abundency[ax_id].set_xticks(np.power(10, np.arange(np.log10(iterations)+1, step=2)))
        ax_abundency[ax_id].set_xlabel('Time (# Additions)', labelpad=0)


        #  Histogram axes
        ax_distribution[ax_id].set_xticks(range(6))
        ax_distribution[ax_id].set_xlim(-0.5, 5.5)
        ax_distribution[ax_id].set_xlabel('# RM', labelpad=3)

        ax_distribution[ax_id].set_yticks(range(2))
        ax_distribution[ax_id].set_ylim(0, 1)


    fig_path = os.path.join(get_path(), 'Figure_S12')
    make_dirs(fig_path)
    fig.savefig(os.path.join(fig_path, f'figS12_{fig_id}.png'), bbox_inches = 'tight')