import numpy as np

import matplotlib.pyplot as plt

import os

from helpers import default_parameters, make_dirs, get_path, set_rc_params
from helpers import add_figure_labels

set_rc_params()


# Load parameters
C, eta, alpha, beta, delta, T, lb, ub, S, f, iterations = default_parameters()

# Get the path
data_path = get_path('data')

# Prepare figure
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 8/2))
plt.subplots_adjust(wspace=0.35)
axes = axes.flatten()

for ax in axes :
    ax.set_xlabel('$\gamma$')

axes[0].set_ylabel('$\\omega$')


D = 100
T = 10


# Hierarchical triplets -- Helper functions
k1_t = lambda o, D, T :  1 + o * (D - T - 1)
k2_t = lambda o, D, T : (2 + o * (D - T - 2)) * o
k3_t = lambda o, D, T : 1 - (D - T) * delta / C / (eta * (beta * k1_t(o, D, T) - (D - T)))

a_t = lambda o, D, T :   k1_t(o, D, T) * k3_t(o, D, T)
b_t = lambda o, D, T : - k2_t(o, D, T) * k3_t(o, D, T)
c_t = lambda o, D, T : alpha * (k2_t(o, D, T) - k1_t(o, D, T))

s_t = lambda o, D, T : b_t(o, D, T) ** 2 - 4 * a_t(o, D, T) * c_t(o, D, T)

f_t = lambda o, D, T : (-b_t(o, D, T) + np.sqrt(s_t(o, D, T))) / (2 * a_t(o, D, T))
f_t = np.vectorize(f_t)



# Looped triplets -- Helper functions
k1_l = lambda o, D, T :  1 + o * (D - 3*T - 1)
k2_l = lambda o, D, T :   o**2 * (D - 3*T)
k3_l = lambda o, D, T : 1 - (D - 3*T) * delta / C / (eta * (beta * k1_l(o, D, T) - (D - 3*T)))

a_l = lambda o, D, T :   k1_l(o, D, T) * k3_l(o, D, T)
b_l = lambda o, D, T : - k2_l(o, D, T) * k3_l(o, D, T)
c_l = lambda o, D, T : alpha * (k2_l(o, D, T) - k1_l(o, D, T))

s_l = lambda o, D, T : b_l(o, D, T) ** 2 - 4 * a_l(o, D, T) * c_l(o, D, T)

f_l = lambda o, D, T : (-b_l(o, D, T) + np.sqrt(s_l(o, D, T))) / (2 * a_l(o, D, T))
f_l = np.vectorize(f_l)


ww = np.logspace(-4.1, 0.1, 50)


with np.errstate(invalid='ignore') :
    for i, T in enumerate([1, 3, 10, 33]) :
        axes[0].plot(f_t(ww, D, T), ww, color = plt.cm.tab10(i))
        axes[1].plot(f_l(ww, D, T), ww, color = plt.cm.tab10(i), label=f'T = {T}')

for ax in axes :
    ax.plot(f_t(ww, 3, 1), ww, 'k--', zorder = 0)
    ax.set_yscale('log')
    ax.set_box_aspect(1)
    ax.set_xlim(0,    1)
    ax.set_ylim(1e-4, 1)


# Add the figure labels
add_figure_labels(['A', 'B'], axes, dx=-0.07, dy=0.036)

axes[1].legend(bbox_to_anchor=(1.05, 1.03), loc='upper left').get_frame().set_edgecolor('k')

fig_path = os.path.join(get_path(), 'Figure_S8')
make_dirs(fig_path)
fig.savefig(os.path.join(fig_path, 'figS8.png'), bbox_inches = 'tight')