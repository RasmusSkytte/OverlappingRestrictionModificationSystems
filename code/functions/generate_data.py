import numpy as np
import os
from tqdm import tqdm
from p_tqdm import p_umap
from functools import partial


from helpers import default_parameters, make_dirs, get_path, pickle_write, simulate_original_model, simulate_extended_model, dynamical_system, get_RMs, get_solver

num_cores = 3
n_batches = 20

tqdm_disable = False

if num_cores > 1 :
    tqdm_disable = True




def generateData(k, nFigS14_S16=None, nFigS15=None, nFig5=None) :

    # Start with default parameters
    C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations = default_parameters()

    data_path = get_path('data')

    # Check which scenario to generate data for
    if k < np.prod(nFigS14_S16) + np.prod(nFigS15) :

        # Overwrite T
        if k < np.prod(nFigS14_S16) :
            r_div = nFigS14_S16[0]
            t = int(k / r_div)

            Ts = np.logspace(1, 4, nFigS14_S16[1])
            T = Ts[t]

            iterations = int(1e4)   # Overwrite the iterations

        # Overwrite C
        else :
            r_div = nFigS15[0]
            c = int((k - np.prod(nFigS14_S16)) / r_div) # Determine the C

            Cs = np.logspace(7, 9, nFigS15[1])
            C = Cs[c]

            iterations = int(1e5)   # Overwrite the iterations

        r = k % r_div      # Determine run
        s = 0              # Determine the seed

        if r in [0, 2, 3, 4] :
            return


        sdir = f'iterations_1e{np.log10(iterations):.0f}_T_1e{np.log10(T):.1f}_C_1e{np.log10(C):.1f}'
        make_dirs(os.path.join(data_path, f'FigS14_S15_S16', sdir))

        # Set seed
        np.random.seed(s)

        # Define run parameters
        sampleTimes = np.power(10, np.arange(np.log10(iterations)+1))

        if r == 0 :
            fname = os.path.join(data_path, f'FigS14_S15_S16', sdir, f'RM_inf_seed_{s}.pkl')
            if not os.path.exists(fname) :
                B_end, P_end, _, bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age, gamma_all, omega_all = simulate_original_model(Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, fname, sampleTimes, disable=tqdm_disable, solver=get_solver())
                pickle_write(fname, B_end, P_end, bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age, gamma_all, omega_all, Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, sampleTimes)

        else :
            RMs = np.arange(50 * 2**(r-1))
            fname = os.path.join(data_path, f'FigS14_S15_S16', sdir, f'RM_{len(RMs)}_seed_{s}.pkl')
            if not os.path.exists(fname) :
                B, _, B_end, P_end, _, bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age, B_all = simulate_extended_model(Alpha, Beta, Eta, Delta, C, T, S, RMs, f, lb, ub, iterations, fname, sampleTimes, disable=tqdm_disable, solver=get_solver())
                pickle_write(fname, B, B_end, P_end, bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age, B_all, Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, sampleTimes)




    elif k >= np.prod(nFigS14_S16) + np.prod(nFigS15) and k < (np.prod(nFigS14_S16) + np.prod(nFigS15) + np.prod(nFig5)) :

        make_dirs(os.path.join(data_path,  'Fig5'))

        k -= np.prod(nFigS14_S16) + np.prod(nFigS15)

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
nFigS14_S16 = (6,  4)  # 6 different K, 4 values of T
nFigS15  = (6,  5)  # 6 different K, 5 values of C
nFig5   = (17, 17) # 17 x 17 phasespace

param_space = np.arange(np.prod(nFigS14_S16) + np.prod(nFigS15) + np.prod(nFig5))

# Perform maximization locally
if num_cores == 1 :
    for k in tqdm(param_space) :
        generateData(k, nFigS14_S16=nFigS14_S16, nFigS15=nFigS15, nFig5=nFig5)

elif __name__ == '__main__':
    p_umap(partial(generateData, nFigS14_S16=nFigS14_S16, nFigS15=nFigS15, nFig5=nFig5), param_space, num_cpus=num_cores)
