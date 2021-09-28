import enum
import numpy as np

import matplotlib.pyplot as plt
import scipy.io

import os

from helpers import default_parameters, make_dirs, get_path, set_rc_params, add_figure_labels


set_rc_params()

K = 50
seed = 5


# Load parameters
C, _, _, beta, _, _, _, _, _, _, iterations = default_parameters()

data_path = os.path.join(get_path('data'))


fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))


# Add the simulation dynamics to the plot
RMs = np.arange(K)

lname = os.path.join(data_path, 'Fig6', f'RM_{K}_seed_{seed}.mat')
label = f'${K=}$'

if os.path.exists(lname) :
    data = scipy.io.loadmat(lname)
    data['nRM'] = [len(lst[0, :].tolist()) for lst in data['B_samples'][0, :][-1].flatten().tolist()] # wtf scipy?

else :
    raise ValueError(f'{lname} is missing!')




ax_annotate = fig.add_subplot(111, facecolor='none')

xa = -0.1  # X coord of the annotaiton
y0 = -0.35
ax_annotate.text(xa,  y0+0.67,  'P / C',   verticalalignment = 'center', horizontalalignment = 'left', rotation = 90, color='r')
ax_annotate.text(xa,  y0+0.82, ', B / C,', verticalalignment = 'center', horizontalalignment = 'left', rotation = 90, color='k')
ax_annotate.text(xa,  y0+1.0,  'D / B',    verticalalignment = 'center', horizontalalignment = 'left', rotation = 90, color='b')
ax_annotate.set_axis_off()


for ax, T0 in zip(axes, [4_000, iterations]) :

    T = np.arange(T0 - 1000, T0)
    P = data['phages'].flatten() / C
    B = data['bacteria'].flatten() / C
    D = data['diversity'].flatten() / beta

    ax.plot(T, P[T], 'r')
    ax.plot(T, B[T], 'k')
    ax.plot(T, D[T], 'b')
    ax.set_xlabel('Time (# Additions)')

    ax.set_yticks(range(0, 25, 5))
    ax.set_ylim(0, 17)


# Add the figure labels
add_figure_labels(['A', 'B'], axes, dy=0.025)


fig_path = os.path.join(get_path(), 'Figure_S11')
make_dirs(fig_path)
fig.savefig(os.path.join(fig_path, 'figS11.png'), bbox_inches = 'tight')