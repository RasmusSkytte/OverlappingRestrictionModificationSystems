import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

from helpers import default_parameters, get_RMs, dynamical_system, make_dirs, get_path, get_solver, set_rc_params, add_figure_labels

set_rc_params()

# Set seed
np.random.seed(60)

# Number of repititions
M = 30

# Load parameters
C, eta, alpha, beta, delta, T, lb, ub, S, f, iterations = default_parameters()
beta = 25
T = 1e4

# Compute SS value
deb = delta / (eta * (beta - 1))

# Determine which case to simulate
B = [[0], [1], [0, 1]]
P = [[]]


# Define test omegas
omegas = np.power(10, [-0.25, -0.5, -1, -2])

# Compute number of RM systems
RMs = get_RMs(B)

# Number of bacteria and phages
nB = len(B)
nP = len(P)


# Prepare figure
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 6))
axes = axes.flatten()
axes[2].remove()
axes = axes[[0, 1, 3]]


# Adjust position of axes
tmp0 = axes[0].get_position()
tmp1 = axes[1].get_position()
tmp2 = axes[2].get_position()

dx = 0.1
axes[1].set_position([tmp1.xmin+dx, tmp1.ymin, tmp1.xmax-tmp1.xmin-dx, tmp1.ymax-tmp1.ymin])
axes[2].set_position([tmp2.xmin+dx, tmp2.ymin, tmp2.xmax-tmp2.xmin-dx, tmp2.ymax-tmp2.ymin])
axes[0].set_position([tmp0.xmin,    tmp2.ymin, tmp0.xmax-tmp0.xmin+dx, tmp0.ymax-tmp2.ymin])

# Add the figure labels
add_figure_labels(['A'],     [axes[0]], dx=-0.05, dy=0.025)
add_figure_labels(['B', 'C'], axes[1:], dx=-0.03, dy=0.025)

# Set the labels
for ax in axes :
    ax.set_xlabel('$\gamma$')

axes[0].set_ylabel('$B$')
axes[1].set_ylabel('$b_A, b_B$')
axes[2].set_ylabel('$b_{AB}$')

for ax in axes :
    ax.set_xlim(0, 1)

ylims = [np.round((3*deb) / 10**np.floor(np.log10(3 * deb)), 2) * 10**np.floor(np.log10(3 * deb)),
         np.ceil( (  deb) / 10**np.floor(np.log10(    deb)))    * 10**np.floor(np.log10(    deb)),
         np.ceil( (  deb) / 10**np.floor(np.log10(    deb)))    * 10**np.floor(np.log10(    deb))]

for ax, ymax in zip(axes, ylims) :
    ax.set_ylim(0, ymax)

for ax in axes[1:] :
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

symbols = ['v', 's', 'o', '^']


# Guide lines
axes[0].plot([alpha, alpha], axes[0].get_ylim(), 'k-.', linewidth = 1.5)
axes[1].plot([alpha, alpha], axes[1].get_ylim(), 'k-.', linewidth = 1.5, label = '$\gamma = \\alpha$')
axes[2].plot([alpha, alpha], axes[2].get_ylim(), 'k-.', linewidth = 1.5)

f = lambda g, o : (g**2 * (1 + o) - 2 * o * g) * (1 - 2 * delta / (beta * (1 + o) - 2)) - alpha * (1 - o)
g = lambda o : scipy.optimize.root(lambda g : f(g, o), 1)
gg = np.array([g(o).x[0] for o in omegas])



label = '$\gamma = f(\omega)$'
for g in gg :
    axes[0].plot([g, g], axes[0].get_ylim(), 'k--', linewidth = 1.5)
    axes[1].plot([g, g], axes[1].get_ylim(), 'k--', linewidth = 1.5, label = label)
    axes[2].plot([g, g], axes[2].get_ylim(), 'k--', linewidth = 1.5)
    label = '_nolegend_'


axes[0].plot(axes[0].get_xlim(),    [deb,    deb], 'k', linewidth = 1.5)
axes[1].plot(axes[1].get_xlim(), [10*deb, 10*deb], 'k', linewidth = 1.5, label = '$\\frac{\delta}{\eta(\\beta-1)}$')
axes[2].plot(axes[2].get_xlim(),    [deb,    deb], 'k', linewidth = 1.5)



label = '$\\frac{2\delta}{\eta(\\beta(1+\\omega)-2)}$'
for o in omegas :
    tmp = delta / (eta * (beta * (1 + o) - 2))
    axes[0].plot(axes[0].get_xlim(), [2*tmp, 2*tmp], 'k:', linewidth = 1.5)
    axes[1].plot(axes[1].get_xlim(),   [tmp,   tmp], 'k:', linewidth = 1.5, label = label)
    label = '_nolegend_'



# Loop over omegas
for j, o in tqdm(enumerate(omegas), total=len(omegas)) :
    label = f'$\omega = 10^{{{np.log10(o)}}}$'

    # Loop over gammas
    for r in tqdm(np.linspace(0, 1, M), total=M, leave=False) :

        # delta eta beta limit
        x0 = deb * np.concatenate((np.ones(nB) / nB, 10 / (nB * nP) * np.ones(nB*nP)))


        omega_0 = np.ones_like(RMs) * o
        cost    = np.ones_like(RMs) * r

        # Run dynamics
        _, _, y, t, res = dynamical_system('Extended', B, P, omega_0=omega_0, cost=cost, params=(alpha, beta, eta, delta, C, T), x0 = x0, solver=get_solver(), rtol = 1e-4, atol = 1e-7)

        if t[-1] < T :
            B_end = y[:nB, -1]
        else :
            B_end = np.median(res.sol(np.linspace(T-10, 10, 100))[:nB, :], axis=1)

        axes[0].scatter(cost[0],  B_end.sum(), marker = symbols[j], color = plt.cm.tab10(j), zorder = 3, label = label)
        axes[1].scatter(cost[:2], B_end[:2],   marker = symbols[j], color = plt.cm.tab10(j), zorder = 3)
        axes[2].scatter(cost[0],  B_end[-1],   marker = symbols[j], color = plt.cm.tab10(j), zorder = 3)

        label = '_nolegend_'

leg = axes[0].legend().get_frame()
leg.set_edgecolor('k')
leg.set_alpha(1.0)
leg.set_linewidth(1.0)

leg = axes[1].legend(bbox_to_anchor=(-0.28, -0.5)).get_frame()
leg.set_edgecolor('k')
leg.set_alpha(1.0)
leg.set_linewidth(1.0)

fig_path = os.path.join(get_path(), 'Figure_S7')
make_dirs(fig_path)
fig.savefig(os.path.join(fig_path, 'figS7.png'), bbox_inches = 'tight')