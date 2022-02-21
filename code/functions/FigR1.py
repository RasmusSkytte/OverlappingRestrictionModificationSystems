import numpy as np
import matplotlib.pyplot as plt

from helpers import make_dirs, get_path, set_rc_params
from helpers import analyze_sequence_data, pickle_write, pickle_read
from helpers import plot_summery_figure, add_figure_labels

import os

set_rc_params()

method = ''
#method = 'hit_miss'
#method = 'linkswapping'
#method = 'configuration_model'

conserve_RM_degree = True

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_R1')
data_path = os.path.join(get_path('data'), 'Fig4_6')
make_dirs(fig_path)
make_dirs(data_path)

# Check if data has been generated
fname   = 'data_sequences'
figname = 'figR1'
if method != '' :
    fname   += f'_{method}'
    figname += f'_{method}'
else :
    method = 'hit_miss'
    conserve_RM_degree = False
fname   += '.pkl'
figname += '.png'

lpath = os.path.join(data_path, fname)
if not os.path.exists(lpath) :

    # Generate and store the data
    pickle_write(lpath, *analyze_sequence_data(conserve_RM_degree = conserve_RM_degree, method = method, max_attempts = 1_000))

# Load the data
_, diff_random_o_ij, diff_random_clustering, diff_random_robins, diff_random_cycles, number_of_cycles, RMs_abundence, average_RM_abundence, names = pickle_read(lpath)

# Prepare figures
fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8, 12.5))
axes_clustering, axes_robins, axes_diff_cycles, axes_cycles, axes_RMs = axes

# Reoder the data
diff_random_o_ij       = [x for _, x in sorted(zip(average_RM_abundence, diff_random_o_ij))]
diff_random_clustering = [x for _, x in sorted(zip(average_RM_abundence, diff_random_clustering))]
diff_random_robins     = [x for _, x in sorted(zip(average_RM_abundence, diff_random_robins))]
diff_random_cycles     = [x for _, x in sorted(zip(average_RM_abundence, diff_random_cycles))]
number_of_cycles       = [x for _, x in sorted(zip(average_RM_abundence, number_of_cycles))]
RMs_abundence          = [x for _, x in sorted(zip(average_RM_abundence, RMs_abundence))]
names                  = [x for _, x in sorted(zip(average_RM_abundence, names))]


# Remove data where there are too few samples
hits = np.array(list(map(len, diff_random_o_ij)))
for i in np.argwhere(hits < hits.max()).flatten() :
    diff_random_o_ij[i]       = np.array([np.nan])
    diff_random_clustering[i] = np.array([np.nan])
    diff_random_robins[i]     = np.array([np.nan])
    diff_random_cycles[i]     = np.array([np.nan])
    number_of_cycles[i]       = np.array([np.nan])

print('-------------- Data report -----------------')
print(f'{names[-1]} strains have on average {np.round(np.mean(RMs_abundence[-1]), 1)} RN systems with standard deviation of +- {np.round(np.std(RMs_abundence[-1]), 1)}')
print(f'{names[-1]} strains have {number_of_cycles[-1]} basic cycles')

plot_summery_figure(axes, names, zip(diff_random_clustering, diff_random_robins, diff_random_cycles, number_of_cycles, RMs_abundence), disable_significance=[False, False, False, True, True], sort=False)
lpad = 4


#axes_clustering.set_yticks(np.arange(-1, 1, 0.1))
axes_clustering.set_ylim(-0.4, 0.6)
axes_clustering.set_ylabel('$C^{Latapy} - C^{Latapy}_{rand}$', labelpad=lpad)

axes_robins.set_ylim(-0.4, 0.6)
axes_robins.set_ylabel('$C^{robins} - C^{robins}_{rand}$', labelpad=lpad)

axes_diff_cycles.set_yticks(np.arange(-5, 20, 5))
axes_diff_cycles.set_ylim(-8, 16)
axes_diff_cycles.set_ylabel('$N^{Cyc} - N^{Cyc}_{rand}$', labelpad=lpad)

#axes_cycles.set_yticks(np.arange(-5, 20, 5))
axes_cycles.set_ylim(-5, 200)
axes_cycles.set_ylabel('$N^{Cyc}$', labelpad=lpad)

pos = axes_RMs.get_position()
pos.y0 += 0.0675
axes_RMs.set_position(pos)

axes_RMs.set_yticks(np.arange(0, 20, 2.5))
axes_RMs.set_ylim(0, 5)
axes_RMs.set_ylabel('#RM', labelpad=lpad+2)


# Add the figure labels
add_figure_labels(['A', 'B', 'C', 'D', 'E'], axes, dx=-0.075, dy=0.01)

fig.savefig(os.path.join(fig_path, figname), bbox_inches = 'tight')