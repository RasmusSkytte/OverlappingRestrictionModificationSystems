import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.io

import os

from helpers import default_parameters, make_dirs, get_path, set_rc_params, add_figure_labels, plot_summery_figure
from helpers import analyze_sequence_data, analyze_simulation_data, pickle_write, pickle_read, conserve_RM_degree


set_rc_params()

K = 50



# Load parameters
C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations = default_parameters()

data_path = os.path.join(get_path('data'))


# Prepare figure
fig = plt.figure(figsize=(8, 9))

dx = 0.12
dy = 0.05

x0 = 0.0    # X coord of first group
y0 = 0.45   # Y coord of the first group
w0 = 0.4    # width of first group
h0 = 1-y0

x1 = x0 + w0 + dx  # X coord of second group
y1 = y0-0.035      # Y coord of the second group
w1 = 1-x1          # width of second group
h1 = 1-y1

x2 = x0   # X coord of third group
y2 = 0    # Y coord of the third group
w2 = 0.8 * (1-x0)

dx_c = 0.01     # Distance between cut-axes


k = 1 / (1-w2) # Width ratios of the cut-axes
s = 0.4        # scale of insert
x_off = 0.03   # x offset of the inset
y_off = 0.06

ds = 0.05
h = np.array([0.2, 0.2, 0.6])
assert(h.sum() <= 1)
h *= (h0-ds)

ax_distribution = fig.add_axes([x0, y0,                 w0, h[0]])
ax_abundency    = fig.add_axes([x0, y0+ds+h[0],         w0, h[1]])
ax_dynamics     = fig.add_axes([x0, y0+ds+h[:2].sum(),  w0, h[2]], sharex=ax_abundency)

ax_annotate     = fig.add_axes([x0, y0, 1, h0], facecolor='none')


ds = 0.02
h = np.array([0.14, 0.43, 0.43])
assert(h.sum() <= 1)
h *= (h1-2*ds)

axes_RMs       = fig.add_axes([x1, y1,                  w1, h[0]])
axes_unique    = fig.add_axes([x1, y1+ds+h[0],          w1, h[1]])
axes_o_ij      = fig.add_axes([x1, y1+2*ds+h[:2].sum(), w1, h[2]])


ax_data_1      = fig.add_axes([x2,                         y2,                      w2-dx_c/2,         y0-dy])
ax_data_2      = fig.add_axes([x2+w2+dx_c/2,               y2,                      1/k-dx_c/2,        y0-dy])
ax_inset_1     = fig.add_axes([x_off+x2+(w2-dx_c/2)*(1-s), y0-dy-s*(y0-dy)-y_off,   (w2-dx_c/2)*s,  s*(y0-dy)])
ax_inset_2     = fig.add_axes([x_off+x2+(w2+dx_c/2),       y0-dy-s*(y0-dy)-y_off,   (1/k-dx_c/2)*s, s*(y0-dy)])



xa = -0.05 # X coord of the annotaiton
ax_annotate.text(xa, 0.64, 'P / (10*C)',  verticalalignment = 'center', horizontalalignment = 'right', rotation = 90, color='r')
ax_annotate.text(xa, 0.79, ', B / C,',    verticalalignment = 'center', horizontalalignment = 'right', rotation = 90, color='k')
ax_annotate.text(xa, 0.91,  'D / β',      verticalalignment = 'center', horizontalalignment = 'right', rotation = 90, color='b')

ax_annotate.text(xa, 0.36, '$\langle $#$RM\\rangle$', verticalalignment = 'center', horizontalalignment = 'right', rotation = 90)

ax_annotate.text(xa, 0.09, 'pmf.', verticalalignment = 'center', horizontalalignment = 'right', rotation = 90)
ax_annotate.set_axis_off()



# Add the figure labels
add_figure_labels(['A', 'C'], [ax_dynamics, ax_distribution],           dx=-0.07, dy=0.014)
add_figure_labels(['B'], [ax_abundency],                                dx=-0.07, dy=0.005)
add_figure_labels(['D', 'E', 'F'], [axes_o_ij, axes_unique, axes_RMs],  dx=-0.08, dy=0.015)
add_figure_labels(['G'], [ax_data_1], dx=-0.07, dy=0.015)


# Add the simulation dynamics to the plot
RMs = np.arange(K)
seed = 5

lname = os.path.join(data_path, 'Fig6', f'RM_{K}_seed_{seed}.mat')

if os.path.exists(lname) :
    data = scipy.io.loadmat(lname)
    data['nRM'] = [len(lst[0, :].tolist()) for lst in data['B_samples'][0, :][-1].flatten().tolist()] # wtf scipy?

else :
    raise ValueError(f'{lname} is missing!')

# Plot dynamics
ax_dynamics.plot(np.arange(iterations), data['phages'].flatten()    / (10 * C), 'r', label='P / (10 C)')
ax_dynamics.plot(np.arange(iterations), data['bacteria'].flatten()  / C,        'k', label='C')
ax_dynamics.plot(np.arange(iterations), data['diversity'].flatten() / Beta,     'b', label='D / $\\beta$')

#ax_dynamics.text(3e4, 2.7, f'{K=}')

# Plot the average RM count
ax_abundency.plot(np.arange(iterations), data['mRM'].flatten(), color=plt.cm.tab10(4))

# Plot the histogram
ax_distribution.hist(data['nRM'], bins=np.arange(0.5, 6.5), color=plt.cm.Pastel1(2), density=True)



# Dynamics axes
ax_dynamics.set_yticks(range(4))

ax_dynamics.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax_dynamics.set_xlim(1, iterations)
ax_dynamics.set_xticks(np.power(10, np.arange(np.log10(iterations)+1)))
ax_dynamics.set_xscale('log')


# Average RM axes
ax_abundency.set_yticks(range(0, 5, 2))
ax_abundency.set_ylim(0, 4.5)

ax_abundency.set_xlim(1, iterations)
ax_abundency.set_xticks(np.power(10, np.arange(np.log10(iterations)+1)))
ax_abundency.set_xscale('log')
ax_abundency.set_xlabel('Time (# Additions)', labelpad=0)


#  Histogram axes
ax_distribution.set_xticks(range(6))
ax_distribution.set_xlim(-0.5, 5.5)
ax_distribution.set_xlabel('# RM', labelpad=0)

ax_distribution.set_yticks(range(2))
ax_distribution.set_ylim(0, 1)




# Check if data has been generated
lpath = os.path.join(data_path, 'Fig4_6', 'data_simulations.pkl')
if not os.path.exists(lpath) :

    # Generate and store the data
    pickle_write(lpath, *analyze_simulation_data(conserve_RM_degree = conserve_RM_degree()))

# Load the data
diff_random_unique, diff_random_o_ij, RMs_abundence, average_RM_abundence, names, iters = pickle_read(lpath)

# Exclude earlier iterations
its = [3, 4, 5, 6]
average_RM_abundence = [x for x, it in zip(average_RM_abundence, iters) if it in its]
diff_random_o_ij     = [x for x, it in zip(diff_random_o_ij,     iters) if it in its]
diff_random_unique   = [x for x, it in zip(diff_random_unique,   iters) if it in its]
RMs_abundence        = [x for x, it in zip(RMs_abundence,        iters) if it in its]
names                = [x for x, it in zip(names,                iters) if it in its]

# Compute metrics for later
# Plot the correlation significant summery metrics
x2  = np.array([np.mean(s) for s in diff_random_o_ij])
dx2 = np.array([ np.std(s) for s in diff_random_o_ij])
y2  = np.array([np.mean(s) for s in diff_random_unique])
dy2 = np.array([ np.std(s) for s in diff_random_unique])

# Color by the iterations
c2 = np.array(average_RM_abundence)


# Filter out data so only 1 seed remains
average_RM_abundence = [x for x, name in zip(average_RM_abundence, names) if int(name.split(':')[-1]) == seed]
diff_random_o_ij     = [x for x, name in zip(diff_random_o_ij,     names) if int(name.split(':')[-1]) == seed]
diff_random_unique   = [x for x, name in zip(diff_random_unique,   names) if int(name.split(':')[-1]) == seed]
RMs_abundence        = [x for x, name in zip(RMs_abundence,        names) if int(name.split(':')[-1]) == seed]
names                = [name.replace('RM', 'K') for name in names         if int(name.split(':')[-1]) == seed]


# Reoder the data
diff_random_o_ij   =  [x for _, x in sorted(zip(average_RM_abundence, diff_random_o_ij),    key=lambda x: x[0])]
diff_random_unique =  [x for _, x in sorted(zip(average_RM_abundence, diff_random_unique),  key=lambda x: x[0])]
RMs_abundence      =  [x for _, x in sorted(zip(average_RM_abundence, RMs_abundence),       key=lambda x: x[0])]
names              =  [x for _, x in sorted(zip(average_RM_abundence, names),               key=lambda x: x[0])]


plot_summery_figure([axes_o_ij, axes_unique, axes_RMs], names, zip(diff_random_o_ij, diff_random_unique, RMs_abundence), disable_significance=[False, False, True], plot_names=False, sort=False)
lpad = -6

axes_o_ij.set_yticks(np.arange(-1, 1, 0.1))
axes_o_ij.set_ylim(-0.15, 0.3)
axes_o_ij.set_ylabel('$\langle I\;/\;U \\rangle$ - $\langle I\;/\;U\\rangle_{rand}$', labelpad=lpad+6, y = 0.5)

axes_unique.set_yticks(np.arange(-1, 1, 0.1))
axes_unique.set_ylim(-0.15, 0.4)
axes_unique.set_ylabel('$f^{ u} - f^{ u}_{rand}$', labelpad=lpad)


axes_RMs.set_yticks(np.arange(0, 20, 2.5))
axes_RMs.set_ylim(0, 5)
axes_RMs.set_ylabel('#RM', labelpad=lpad+12, y = 0.45)




# Check if data has been generated
lpath = os.path.join(data_path, 'Fig4_6', 'data_sequences.pkl')
if not os.path.exists(lpath) :

    # Generate and store the data
    pickle_write(lpath, *analyze_sequence_data(conserve_RM_degree = conserve_RM_degree()))

# Load the data
diff_random_unique, diff_random_o_ij, _, _, _, _, _, average_RM_abundence, _ = pickle_read(lpath)

# Plot the correlation significant summery metrics
x1  = np.array([np.mean(s) for s in diff_random_o_ij])
dx1 = np.array([ np.std(s) for s in diff_random_o_ij])
y1  = np.array([np.mean(s) for s in diff_random_unique])
dy1 = np.array([ np.std(s) for s in diff_random_unique])

# Color by the average RM abundence
c1 = np.array(average_RM_abundence)


# Set the cutoff
split = 0.103

for plot_id, (axes, (x, dx, y, dy), coloring) in enumerate(zip([(ax_data_1, ax_data_2), (ax_inset_1, ax_inset_2)], [(x1, dx1, y1, dy1), (x2, dx2, y2, dy2)], [c1, c2])) :
    for ax_id, ax in enumerate(axes) :

        if ax_id == 0 :
            s = x <= split
        else :
            s = x > split

        ax.errorbar(x[s], y[s], xerr=dx[s], yerr=dy[s], ls='none', zorder=0, ecolor='k', elinewidth=0.5)

        cm_p = 5
        cm_m = 1
        cmap = 'gist_rainbow'

        h = ax.scatter(x[s], y[s], c=coloring[s], vmin=cm_m, vmax=cm_p, s=20, marker='o', linewidths=0.5, edgecolors='k', cmap=cmap)


        ax.plot([0, 0],      [-0.1, 0.4], ':k', zorder=0, linewidth=2)
        ax.plot([-0.1, 0.3], [0, 0],      ':k', zorder=0, linewidth=2)

        ax.set_yticks(np.arange(-1, 1, 0.1))

        ax.set_ylim(-0.1, 0.4)


    cax = inset_axes(axes[-1], width='350%', height='3%', loc='upper right')
    cbar = fig.colorbar(h, cax=cax, orientation='horizontal', ticks = np.arange(cm_m, cm_p+1))




# Adjust the axes
for ax_left, ax_right in zip([ax_data_1, ax_inset_1], [ax_data_2, ax_inset_2]) :

    ax_left.set_xticks(np.arange(-1, 1, 0.05))
    ax_right.set_xticks(np.arange(-1, 1, 0.1))

    ax_left.set_xlim(-0.035, split)
    ax_right.set_xlim( split,   0.3)

    ax_left.spines['right'].set_visible(False)
    ax_right.spines['left'].set_visible(False)


    d = 0.015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax_left.transAxes, color='k', clip_on=False, linewidth=1)
    ax_left.plot((1-d/2, 1+d/2),  (-d,  +d), **kwargs)
    ax_left.plot((1-d/2, 1+d/2), (1-d, 1+d), **kwargs)

    kwargs.update(transform=ax_right.transAxes)  # switch to the bottom axes
    ax_right.plot((-k*d/2, +k*d/2),  (-d, +d),  **kwargs)
    ax_right.plot((-k*d/2, +k*d/2), (1-d, 1+d), **kwargs)

    ax_left.tick_params(right=False)
    ax_right.tick_params(labelleft=False, left=False)

# Manually adjust right inset axes
ax_inset_2.set_xticks([0.3])


ax_data_1.set_xlabel('$\langle I\;/\;U \\rangle  - \langle I\;/\;U \\rangle_{rand}$ ', x = 0.65)
ax_data_1.set_ylabel('$f^{\;u} - f^{\;u}_{rand}$', labelpad = -6)



fig_path = os.path.join(get_path(), 'Figure_6')
make_dirs(fig_path)
fig.savefig(os.path.join(fig_path, 'fig6.png'), bbox_inches = 'tight')