import numpy as np
import os
from tqdm import tqdm
from p_tqdm import p_umap
from functools import partial

from helpers import default_parameters, make_dirs, get_path, pickle_write, simulate_original_model, simulate_extended_model, dynamical_system, get_RMs, get_solver

num_cores = 3
tqdm_disable = False

if num_cores > 1 :
    tqdm_disable = True


def generateData(k, nFigS11=None, nFig5=None) :

    # Start with default parameters
    C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations = default_parameters()

    data_path = get_path('data')

    make_dirs(os.path.join(data_path, f'FigS11', f'iterations_1e{np.log10(iterations):.1f}'))
    make_dirs(os.path.join(data_path,  'Fig5'))
    make_dirs(os.path.join(data_path,  'Fig7'))


    # Check which scenario to generate data for
    if k < np.prod(nFigS11) :

        r = k % nFigS11[0]      # Determine run
        s = int(k / nFigS11[0]) # Determine end point

        iterations = 1e4   # Overwrite the iterations

        # Set seed
        np.random.seed(s)

        # Define run parameters
        sampleTimes = np.power(10, np.arange(np.log10(iterations)+1));

        if r == 0 :
            fname = os.path.join(data_path, f'FigS11', f'iterations_1e{np.log10(iterations):.1f}', f'RM_inf_seed_{s}.pkl')
            if not os.path.exists(fname) :
                B_end, P_end, _, bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age, gamma_all, omega_all = simulate_original_model(Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, fname, sampleTimes, disable=tqdm_disable, solver=get_solver())
                pickle_write(fname, B_end, P_end, bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age, gamma_all, omega_all, Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, sampleTimes)

        else :
            if r in [3, 4, 5] :
                return

            RMs = np.arange(25 * 2**(r-1))
            fname = os.path.join(data_path, f'FigS11', f'iterations_1e{np.log10(iterations):.1f}', f'RM_{len(RMs)}_seed_{s}.pkl')
            if not os.path.exists(fname) :
                B, _, B_end, P_end, _, bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age, B_all = simulate_extended_model(Alpha, Beta, Eta, Delta, C, T, S, RMs, f, lb, ub, iterations, fname, sampleTimes, disable=tqdm_disable, solver=get_solver())
                pickle_write(fname, B, B_end, P_end, bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age, B_all, Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, sampleTimes)




    elif k >= np.prod(nFigS11) and k < (np.prod(nFigS11) + np.prod(nFig5)) :

        k -= np.prod(nFigS11)


        # Define lists to sweep over
        omega, gamma = np.meshgrid(np.power(10, np.linspace(0, -2, nFig5[0])), np.linspace(0, 1, nFig5[1]));
        omega = omega.flatten()
        gamma = gamma.flatten()

        # Define omega and gamma
        omega = omega[k]
        gamma = gamma[k]

        # Overwrite the beta value
        Beta = 25

        # Create output name
        fname = os.path.join(data_path, f'Fig5', f'omega_1e{np.log10(omega):.3f}_gamma_{gamma:.4f}.pkl')

        if not os.path.exists(fname) :

            # Define test case
            B = [[0], [1], [0, 1]]

            P = [[]]

            # Compute number of RM systems
            RMs = get_RMs(B, P)

            # Number of bacteria and phages
            nB = len(B)
            nP = len(P)

            # Set up test space
            deb = Delta / (Eta * (Beta - 1))
            deb2 = 1.9 * Delta / (Eta * (Beta * (1 + omega) - 2))    # Use 1.9 instead of 2 to break degeneracy when omega  = 1
            B0 = np.array([deb, deb2, 1e1, 1e3, 1e5, 1e7])

            b0_1, b0_2 = np.meshgrid(B0, B0)
            b0_1 = b0_1.flatten()
            b0_2 = b0_2.flatten()

            # Allocate arrays
            coordinates = np.full((len(b0_1), 2), fill_value=np.nan)
            BBs         = np.full((len(b0_1), 3), fill_value=np.nan)
            PPs         = np.full( len(b0_1),     fill_value=np.nan)

            # Loop over starting conditions
            for j, (b1, b2) in tqdm(enumerate(zip(b0_1, b0_2)), total=len(b0_1), leave=False, disable=tqdm_disable) :

                # Set starting conditions
                x0 = np.array([b1, b1, b2, 10*b1, 10*b1, 10*b2])

                # Define omega and cost
                omega_0 = np.ones_like(RMs) * omega
                cost    = np.ones_like(RMs) * gamma

                # Run dynamics
                _, _, y, _, _ = dynamical_system('Extended', B, P, omega_0=omega_0, cost=cost, params=(Alpha, Beta, Eta, Delta, C, 10*T), x0 = x0, solver=get_solver(), rtol=1e-6, atol=1e-9)

                B_end =            y[:nB, -1]
                P_end = np.reshape(y[nB:, -1],(nP, nB))

                coordinates[j, :] = [b1, b2]
                BBs[j, :] = B_end
                PPs[j]    = P_end.sum(axis=1)

                # Save results
                pickle_write(fname, coordinates, BBs, PPs)




# Define the number of runs for each generator
nFigS11  = (7,  1)  # 7 different K, 1 repeats
nFig5 = (17, 17)  # 17 x 17 phasespace

params = np.arange(np.prod(nFigS11) + np.prod(nFig5))

# Perform maximization locally
if num_cores == 1 :
    for param in tqdm(params) :
        generateData(param, nFigS11=nFigS11, nFig5=nFig5)

elif __name__ == '__main__':
    p_umap(partial(generateData, nFigS11=nFigS11, nFig5=nFig5), params, num_cpus=num_cores)

