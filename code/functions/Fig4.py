import numpy as np
import matplotlib.pyplot as plt

from helpers import make_dirs, get_path, set_rc_params
from helpers import analyze_sequence_data, pickle_write, pickle_read
from helpers import plot_summery_figure, conserve_RM_degree, add_figure_labels

import os

set_rc_params()

# Prepare figure folders
fig_path = os.path.join(get_path(), 'Figure_4')
data_path = os.path.join(get_path('data'), 'Fig4_6')
make_dirs(fig_path)
make_dirs(data_path)

# Check if data has been generated
lpath = os.path.join(data_path, 'data_sequences.pkl')
if not os.path.exists(lpath) :

    # Generate and store the data
    pickle_write(lpath, *analyze_sequence_data(conserve_RM_degree = conserve_RM_degree()))

# Load the data
diff_random_unique, diff_random_o_ij, _, _, _, _, RMs_abundence, average_RM_abundence, names = pickle_read(lpath)

# Prepare figures
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 8.5))
axes_o_ij, axes_unique, axes_RMs = axes

# Reoder the data
diff_random_o_ij   =  [x for _, x in sorted(zip(average_RM_abundence, diff_random_o_ij))]
diff_random_unique =  [x for _, x in sorted(zip(average_RM_abundence, diff_random_unique))]
RMs_abundence      =  [x for _, x in sorted(zip(average_RM_abundence, RMs_abundence))]
names              =  [x for _, x in sorted(zip(average_RM_abundence, names))]


print('-------------- Data report -----------------')
print(f'{names[-1]} strains have on average {np.round(np.mean(RMs_abundence[-1]), 1)} RN systems with standard deviation of +- {np.round(np.std(RMs_abundence[-1]), 1)}')

plot_summery_figure(axes, names, zip(diff_random_o_ij, diff_random_unique, RMs_abundence), disable_significance=[False, False, True], sort=False)
lpad = -3


axes_o_ij.set_yticks(np.arange(-1, 1, 0.1))
axes_o_ij.set_ylim(-0.1, 0.4)
axes_o_ij.set_ylabel('$\langle I\;/\;U \\rangle$ - $\langle I\;/\;U\\rangle_{rand}$', labelpad=lpad+1, y = 0.45)

axes_unique.set_yticks(np.arange(-1, 1, 0.1))
axes_unique.set_ylim(-0.1, 0.5)
axes_unique.set_ylabel('$f^{ u} - f^{ u}_{rand}$', labelpad=lpad-2)


pos = axes_RMs.get_position()
pos.y0 += 0.125
axes_RMs.set_position(pos)

axes_RMs.set_yticks(np.arange(0, 20, 2.5))
axes_RMs.set_ylim(0, 5)
axes_RMs.set_ylabel('#RM', labelpad=lpad+12)


# Add the figure labels
add_figure_labels(['A', 'B', 'C'], axes, dx=-0.075, dy=0.017)

fig.savefig(os.path.join(fig_path, 'fig4.png'), bbox_inches = 'tight')