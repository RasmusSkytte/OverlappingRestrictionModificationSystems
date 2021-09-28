import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import pickle
import os
import operator
import warnings
from tqdm import tqdm, trange
from numba import jit

from itertools import product
from functools import reduce

import scipy
import scipy.io
from scipy.integrate import solve_ivp

# Set random seed
np.random.seed(0)


def default_parameters() :
    #  Define default parameters
    C       = 1e8   # Carrying capacity
    Eta     = 1e-8  # Adsorption rate
    Alpha   = 0.2   # Dilution rate of bacteria
    Beta    = 100   # Burst size of phages
    Delta   = 0.2   # Decay rate of phages
    T       = 1e3   # Time between addition of species
    lb      = -4    # Lower bound for omega
    ub      = 0     # Upper bound for omega
    S       = 5     # Min number of species

    # Compute cost factor
    f = 0.1

    iterations = int(1e6) 	# Number of iterations

    return C, Eta, Alpha, Beta, Delta, T, lb, ub, S, f, iterations


def get_solver() :
    #return 'LSODA'
    return 'RK45'
    #return 'RK23'
    #return 'BDF'
    #return 'Radau'

def conserve_RM_degree() :
    return False

def set_rc_params() :
    plt.rcParams.update({
        'axes.linewidth': 1.0,
        'figure.dpi': 300.0,
        'font.size': 12.0,
        'lines.linewidth' : 2,
        'patch.facecolor' : plt.cm.Pastel1(3),
        'patch.linewidth' : 1.5,
        'patch.edgecolor' : 'k',
        'patch.force_edgecolor' : True})

def label_args() :
    return {'fontweight' : 'bold',
            'horizontalalignment' : 'right',
            'verticalalignment' : 'top'}

def compute_gamma_and_omega(B, P, cost, omega_0, model='Original') :

    if model == 'Original' :
        nB = len(cost)

        gamma = cost

        # Convert omega to matrix
        omega = np.array([omega_0, ]*nB).T
        omega = omega - np.diag(np.diag(omega)) + np.eye(nB)
        omega = omega[np.newaxis, :, :]


    elif model == 'Extended' :

        nB = len(B)
        nP = max(1, len(P))

        # Compute the growth rates
        gamma = np.nan * np.zeros(nB)
        for i, b in enumerate(B) :
            gamma[i] = np.prod(cost[b])

        # Generate a matrix of how the phages are methylated
        meth = []
        # 2D list: first dimension tells which bacterium methylated it
        #          second dimension tells which phage it is
        for b in B_iter(B) :
            meth_b = []
            for p in P :
                meth_b.append(b[~np.isin(b, p)])

            meth.append(meth_b)

        # Generate a matrix of functional RM systems based on which bacteria
        # the phage has been methylated by.
        RMs_eff = []
        # 3D list: first dimenson tells which bacteria is defending
        #          second dimension tells which bacterium methylated the phage
        #          third dimension tells which phage it is
        for bi in B_iter(B) :
            RMs_eff_bi = []
            # Methylated phages
            for bj in B_iter(B) :
                RMs_eff_bj = []
                for p in P :
                    RMs_eff_bj.append(bi[~np.isin(bi, np.union1d(bj, p))])

                RMs_eff_bi.append(RMs_eff_bj)
            RMs_eff.append(RMs_eff_bi)

        # Compute the RM strength against the phage phenotypes
        omega = np.nan * np.zeros((nP, nB, nB))
        # 3D list: first dimenson tells which bacteria is defeinding
        #          second dimension tells which bacterium methylated the phage
        #          third dimension tells whcih phage is is

        for i in range(nB) :
            for j in range(nB) :
                for k in range(nP) :
                    omega[k, i, j] = np.prod(omega_0[RMs_eff[i][j][k]])

    else :
        raise ValueError('Model not known!')

    return gamma, omega

@jit
def df_dt(t, x, nB, nP, gamma, omega, C, alpha, beta, eta, delta) :

    # Enforce nonzero
    x[x<0] = 1e-10

    # Allocate y
    y = np.zeros((1+nP) * nB)

    # Compute B
    B = np.sum(x[:nB])

    # Compute eta * omega * p
    eop = np.zeros((nB, nP))
    for j in range(nP) :
        eop[:, j] = eta * np.dot(omega[j], x[(j+1)*nB + np.arange(nB)])

    # Equations for bacteria
    y[:nB] = (gamma * (1 - B / C) - alpha - np.sum(eop, axis=1) ) * x[:nB]

    # Phage profilation, adsorption and decay
    for j in range(nP) :
        y[(j+1)*nB + np.arange(nB)] = beta * eop[:, j] * x[:nB] -  x[(j+1)*nB + np.arange(nB)] * (eta * B + delta)

    return y

#@jit
def jac(t, x, nB, nP, gamma, omega, C, alpha, beta, eta, delta) :

    # Enforce nonzero
    x[x<0] = 1e-10

    #jac = [df1_dt / d_x1, df1_dt / d_x2, ----
    #      [df2_dt / d_x1, df2_dt / d_x2, ----  ]

    #jac = [db_dt / db,   db_dt / dp]
    #      [dp_dt / db,   dp_dt / dp]

    b = x[:nB]

    B = b.sum()

    # Compute eta * omega * p
    eop = np.zeros((nB, nP))
    for j in range(nP) :
        eop[:, j] = eta * np.dot(omega[j], x[(j+1)*nB + np.arange(nB)])

    db_dt_db = np.repeat(gamma[:,np.newaxis] / C, nB, axis=1)
    np.fill_diagonal(db_dt_db, gamma*(1-alpha-(B + b) / C) - eop.sum(axis=1))

    db_dt_dp = eta*np.concatenate([o.T for o in omega])

    dp_dt_db = np.concatenate([np.diag(m) for m in beta*eop.T], axis=1)

    dp_dt_dp = np.zeros((nP*nB, nP*nB))
    for j, o in enumerate(omega) :
        dp_dt_dp[j*nB:(j+1)*nB, j*nB:(j+1)*nB] = beta * eta * np.repeat(b[np.newaxis, :], nB, axis=0) * o
    dp_dt_dp -= np.eye(nB*nP)*(eta*B+delta)

    return np.concatenate((np.concatenate((db_dt_db, db_dt_dp)), np.concatenate((dp_dt_db, dp_dt_dp))), axis=1)

def dynamical_system(model, B, P, omega_0=None, cost=None, params=None, x0=None, solver='RK45', rtol=1e-4, atol=1e-7) :

    originalModel = False
    extendedModel = False

    if model == 'Original' :
        originalModel = True
    elif model == 'Extended' :
        extendedModel = True
    else :
        raise ValueError('Model not known!')

    userDefinedOmegas = False
    if omega_0 is not None :
            userDefinedOmegas = True

    userDefinedCost = False
    if cost is not None :
        userDefinedCost = True

    userDefinedParams = False
    if params is not None :
        alpha, beta, eta, delta, C, T_end = params
        userDefinedParams = True

    userDefinedStartingPopulations = False
    if x0 is not None :
        userDefinedStartingPopulations = True


    # Define variables
    if not userDefinedParams :
        C      = 1e8
        alpha  = 0.2
        beta   = 25
        eta    = 1 / C
        delta  = 0.2
        T_end  = 1e4

        RMs = get_RMs(B, P)


    # Draw random omegas for the RMs
    if not userDefinedOmegas :
        lb = -4
        ub = 0
        omega_0 = np.power(10, (lb + (ub-lb) * np.random.uniform(size = len(RMs))))


    # Set the cost of the RMs
    if not userDefinedCost :
        cost = 1 - np.random.uniform(size = len(RMs))


    # Get the omegas and gammas to use in the model
    if originalModel :
        gamma, omega = compute_gamma_and_omega(B, P, cost, omega_0, model='Original')
        nB = len(cost)
        nP = max(1, len(P))
    elif extendedModel :
        gamma, omega = compute_gamma_and_omega(B, P, cost, omega_0, model='Extended')
        nB = len(B)
        nP = len(P)


    # Define the dynamical populations
    if not userDefinedStartingPopulations :
        b = np.ones(nB)   / nB     * delta/(eta*beta)     # Each bacterium has equal amounts of biomass initially
        p = np.ones(nB*nP)/(nB*nP) * np.mean(gamma)/eta   # Each phage has equal amounts of biomass initially
        x0 = np.concatenate((b, p))


    # Run the ode solver
    if solver in ['BDF', 'Radau', 'LSODA'] :
        res = solve_ivp(df_dt, [0, T_end], x0, jac=jac, args=(nB, nP, gamma, omega, C, alpha, beta, eta, delta), rtol=rtol, atol=atol, method=solver, dense_output=True)
    else :
        res = solve_ivp(df_dt, [0, T_end], x0,          args=(nB, nP, gamma, omega, C, alpha, beta, eta, delta), rtol=rtol, atol=atol, method=solver, dense_output=True)

    t = res.t
    y = res.y

    B_end = y[:nB, -1]
    P_end = np.reshape(y[nB:, -1],(nP, nB))

    return B_end, P_end, y, t, res

def generate_species(nRM, lb, ub, f, S) :

    # Allocate
    gamma = []
    omega = []
    n     = []

    # Count the current number of species
    nB = len(nRM)

    k = 1
    while len(nRM) < max(S, nB + 1) :

        # If there are less than S species, add completely new species
        if len(nRM) < S :
            m = 1

        else : # If there are more than S species, mutate from existing

            # Draw from poisson distribution
            m = np.random.poisson(np.median(nRM))
            # Ensure we are within bounds
            m = max(1, m)

        # Sample the distribution
        gamma.append(np.prod(1 - f * np.random.uniform(size=m)))
        omega.append(np.prod(np.power(10, (lb + (ub-lb) * np.random.uniform(size=m)))))
        n.append(m)
        nRM = np.append(nRM, n)

        k += 1

    return gamma, omega, n

def generate_species_extended(B, S, RMs) :

    # Count the current number of species
    nB = len(B)

    k = 0
    while len(B) < max(S, nB + 1) :

        # If there are less than S species, add completely new species
        if len(B) < S :
            b = [RMs[np.random.randint(low=0, high=len(RMs))]]

        else : # If there are more than S species, mutate from existing

            # Draw from poisson distribution
            b = np.random.poisson(np.median(list(map(len, B))))

            # Ensure we are within bounds
            b = min(len(RMs), max(1, b))

            # Choose random RM
            b = np.sort(RMs[np.random.permutation(len(RMs))[:b]]).tolist()

        # Check whether it exists
        exist = False
        for i in range(len(B)) :
            if b == B[i] :
                exist = True
        if exist :
            k += 1
            if k == 10_000 :
                raise ValueError('Recursion detected')
            continue

        # Store the bacterium
        B.append(b)

    return B

def simulate_original_model(Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, fname, sampleTimes=None, solver='RK45', disable=False) :

    # Initialize empty lists
    gamma = np.array([])
    omega = np.array([])
    nRM   = np.array([])
    age   = np.array([])

    gamma_all = np.array([])
    omega_all = np.array([])

    bacteria  = np.full(iterations, fill_value=np.nan)
    phages    = np.full(iterations, fill_value=np.nan)
    diversity = np.full(iterations, fill_value=np.nan)
    mRM       = np.full(iterations, fill_value=np.nan)

    B_samples = []

    # Define filenames
    f_checkpoint = '{0}_{2}.{1}'.format(*fname.rsplit('.', 1) + ['checkpoint'])
    f_log        = '{0}_{2}.{1}'.format(*fname.rsplit('.', 1) + ['status']).replace('.pkl','.txt')

    # Main loop
    for i in tqdm(range(iterations), leave=False, disable=disable) :

        # Generate new gamma, omega value
        g, o, n = generate_species(nRM, lb, ub, f, S)

        gamma = np.append(gamma, g)
        omega = np.append(omega, o)
        nRM   = np.append(nRM,   n)

        gamma_all = np.append(gamma_all, g)
        omega_all = np.append(omega_all, o)

        # Store the age
        age = np.append(age, np.zeros(len(gamma)-len(age)))

        # Count the number of species
        nB = len(gamma);

        # Set starting conditions
        if i == 0 :
            x0 = np.ones(2*nB) # Everyone starts with a population of 1
        else :
            b = np.ones(nB-np.sum(I)) # New bacteria and phages start with a population of 1
            x0 = np.concatenate((y[:np.sum(I)], b, np.abs(y[np.sum(I):]), b))

        # Set tolerances
        rtols = [1e-3, 1e-4]
        atols = [1e-6, 1e-7]

        for rtol, atol in zip(rtols, atols) :

            # Run dynamics
            B_end, P_end, y, _, _ = dynamical_system('Original', [], [], omega_0=omega, cost=gamma, params=(Alpha, Beta, Eta, Delta, C, T), x0=x0, solver=solver, rtol=rtol, atol=atol)

            if np.all(B_end < C) :
                break

            print('Running again with smaller tolerances')

        # Apply threshold of survival
        I = B_end > 1

        # Ignore species below threshold
        B_end = B_end[I]

        # Store data
        diversity[i] = np.sum(I)
        bacteria[i]  = np.sum(B_end)
        phages[i]    = np.sum(P_end)
        mRM[i]       = np.mean(nRM)

        # Remove the dead species
        y = y[np.tile(I, 2), -1]
        gamma = gamma[I]
        omega = omega[I]
        nRM   = nRM[I]
        age   = age[I]

        # Increase the age
        age += 1

        # Sample B if required
        if np.isin(i, sampleTimes) :
            B_samples.append(nRM)

        # Store a checkpoint
        if (i+1) % 100 == 0 :
            pickle_write(f_checkpoint, nRM, bacteria, phages, diversity, [], mRM, gamma, omega, age, B_samples, C, Alpha, Beta, Delta, T, lb, ub, S, f, iterations, sampleTimes)

            with open(f_log, 'w') as log_file :
                log_file.write(f'status: {100*(i+1)/iterations:.1f} %')


    # Cleanup
    os.remove(f_checkpoint)
    os.remove(f_log)

    return B_end, P_end, diversity[-1], bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age, gamma_all, omega_all

def simulate_extended_model(Alpha, Beta, Eta, Delta, C, T, S, RMs, f, lb, ub, iterations, fname, sampleTimes=None, solver='RK45', disable=False) :

    # Draw random base omegas and costs
    omega_0 = np.power(10, lb + (ub-lb) * np.random.uniform(size = len(RMs)))
    cost    = 1 - f * np.random.uniform(size = len(RMs))

    # Prepare empty lists
    B = []      # List of bacteria
    B_all = []

    age   = np.array([])

    bacteria  = np.full(iterations, fill_value=np.nan)
    diversity = np.full(iterations, fill_value=np.nan)
    overlap   = np.full(iterations, fill_value=np.nan)
    mRM       = np.full(iterations, fill_value=np.nan)

    B_samples = []


    # Check for pertubation
    P = [[]]
    phages = np.full(iterations, fill_value=np.nan)


    # Define filenames
    f_checkpoint = '{0}_{2}.{1}'.format(*fname.rsplit('.', 1) + ['checkpoint'])
    f_log        = '{0}_{2}.{1}'.format(*fname.rsplit('.', 1) + ['status']).replace('.pkl','.txt')

    # Main loop
    for i in tqdm(range(iterations), leave=False, disable=disable) :

        # Draw one more species
        l = len(B)
        B = generate_species_extended(B, S, RMs)
        B_all.extend([b for b in B[l:] if b not in B_all])

        # Store the age
        age = np.append(age, np.zeros(len(B) - len(age)))

        # Count the number of species
        nB = len(B)
        nP = len(P)

        # Set starting conditions
        if i == 0 :
            x0 = np.ones((1+nP)*nB) # Everyone starts with a population of 1
        else :
            b = np.ones(nB-np.sum(I)) # New bacteria and phages start with a population of 1
            x0 = np.concatenate((y[:np.sum(I)], b))
            for p in range(nP) :
                x0 = np.concatenate((x0, np.abs(y[(p+1)*np.sum(I):(p+2)*np.sum(I)]), b))

        # Set tolerances
        rtols = [1e-3, 1e-4]
        atols = [1e-6, 1e-7]

        for rtol, atol in zip(rtols, atols) :

            # Run dynamics
            B_end, P_end, y, _, _ = dynamical_system('Extended', B, P, omega_0=omega_0, cost=cost, params=(Alpha, Beta, Eta, Delta, C, T), x0=x0, solver=solver, rtol=rtol, atol=atol)

            if np.all(B_end < C) :
                break

            print('Running again with smaller tolerances')

        # Apply threshold of survival
        I = B_end > 1

        # Ignore species below threshold
        B_end = B_end[I]
        P_end = P_end[:, I]

        # Store data
        diversity[i] = np.sum(I)
        bacteria[i]  = np.sum(B_end)
        phages[i]    = np.sum(P_end)
        mRM[i]       = np.mean(list(map(len, B)))

        # Remove the dead species
        y = y[np.tile(I, 1+nP), -1]
        B = [b for i, b in zip(I, B) if i]
        age = age[I]

        # Increase the age
        age += 1


        # Compute the overlap measure
        Bf = flatten_list(B)
        if len(B) > 1 :
            _, RMs_counts = np.unique(Bf, return_counts=True)
            overlap[i] = np.mean((RMs_counts - 1) / (diversity[i]-1))


        # Sample B if required
        if np.isin(i, sampleTimes) :
            B_samples.append(B)

        # Store a checkpoint
        if (i+1) % 100 == 0 :
            pickle_write(f_checkpoint, B, bacteria, phages, diversity, overlap, mRM, cost, omega_0, age, B_samples, C, Alpha, Beta, Delta, T, lb, ub, S, f, iterations, sampleTimes)

            with open(f_log, 'w') as log_file :
                log_file.write(f'status: {100*(i+1)/iterations:.1f} %')

    # Cleanup
    os.remove(f_checkpoint)
    os.remove(f_log)

    return B, P, B_end, P_end, diversity[-1], bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age, B_all

def compute_network(B) :

    # Compute the span of the RM systems
    RMs = get_RMs(B)

    A_ij  = np.zeros((len(B),   len(RMs)))
    A_RMs = np.zeros((len(RMs), len(RMs)))
    A_Bs  = np.zeros((len(B),   len(B)))

    # Presence absence matrix
    for i in range(len(B)) :
        A_ij[i, B[i]] = 1

    # Overlap matrix
    for b in B :
        for ri, rj in product(b, b) :
            A_RMs[ri, rj] += 1


    # Alternative overlap matrix
    for i, bi in enumerate(B) :
        for j, bj in enumerate(B) :
            A_Bs[i, j] = len(np.intersect1d(bi, bj))


    if len(A_RMs) > 0 :

        # Prune RM network
        I = np.sum(A_RMs, axis=0) > 0
        A_RMs = A_RMs[I, :]
        A_RMs = A_RMs[:, I]

        # Prune adjacency matrix network
        A_ij = A_ij[:, np.sum(A_ij, axis=0) > 0]


    return A_ij, A_RMs, A_Bs

def compute_relative_overlaps(A_ij, B, ind = None) :

    # Compute relative overlap between species
    o_ij = np.full((len(B), len(B)), fill_value = np.nan)
    for i in range(len(B)) :
        for j in range(len(B)) :
            if i == j :
                continue

            o_ij[i, j] = np.logical_and(A_ij[i], A_ij[j]).sum() / np.logical_or(A_ij[i], A_ij[j]).sum()

    if ind == 'all' :
        return o_ij
    elif ind is not None :
        return o_ij[ind]
    else :
        return o_ij[np.triu_indices(len(B), k=1)]  # Only upper triangle of matrix (excluding diagonal)

def compute_overlap_measure(A_RMs, A_ij) :
    c_RMs = np.diag(A_RMs)
    return c_RMs / (np.size(A_ij, 0) - 1)  # p_i = # occurances / # number of other strains

def compute_unique_RM_fraction(B) :

    B_tuples = set([tuple(b) for b in B])  # Cunvert to tuples so we can use set
    return np.mean([0 if np.isin(b, flatten_list(list(B_tuples - set([tuple(b)])))).any() else 1 for b in B_tuples])

def generate_random_sample(A_RMs, A_ij, conserve_RM_degree = False) :

    # Generate the available RM systems
    RMs = np.arange(len(A_RMs))

    # Construct the random strain generator
    B_generator = lambda n, RMs : np.random.choice(RMs, size=n, replace=False)

    j = 0
    while j < 10000 :

        # Reset the B degree
        nRM_tmp = A_ij.sum(axis=1).astype(int).tolist()

        # Reset the RM degree
        RMs_degree = A_ij.sum(axis=0)
        RMs_degree = RMs_degree[RMs_degree > 0]

        # Generate strains for the run
        B_rand = []

        k = 0
        while len(B_rand) < len(A_ij) :

            I = np.random.randint(len(nRM_tmp))

            if nRM_tmp[I] > np.sum(RMs_degree > 0) :    # The strain cannot be generated
                break

            new_B = B_generator(nRM_tmp[I], RMs[RMs_degree > 0]).tolist()

            k += 1
            if new_B not in B_rand :
                B_rand.append(new_B)
                nRM_tmp.pop(I)

                if conserve_RM_degree :
                    for r in new_B :
                        RMs_degree[r] -= 1

            if k > 1000 :
                break

        if len(B_rand) == len(A_ij) :
            return B_rand
        j += 1

    raise ValueError('Random network not generated')

def get_metrics_of_real_network(Bs) :

    A_ij_s, A_RMs_s, _ = compute_network(Bs)

    # Get metrics of real network
    avg_overlap  = np.mean(compute_overlap_measure(A_RMs_s, A_ij_s))
    RMs_per_B    = A_ij_s.sum(axis=1)
    avg_degree   = A_RMs_s.sum(axis=1).mean()
    avg_conn     = np.clip(A_RMs_s, 0, 1).sum(axis=1).mean()
    f_unique     = compute_unique_RM_fraction(Bs)
    avg_o_ij     = np.mean(compute_relative_overlaps(A_ij_s, Bs))

    return avg_overlap, RMs_per_B, avg_degree, avg_conn, f_unique, avg_o_ij

def get_metrics_of_random_network(Bs, n_random_sample_size = 1, conserve_RM_degree = False) :

    A_ij_s, A_RMs_s, _ = compute_network(Bs)

    # Prepare output containers for random networks
    avg_overlap_rand  = []
    avg_RM_per_B_rand = []
    avg_degree_rand   = []
    avg_conn_rand     = []
    f_unique_rand     = []
    avg_o_ij_rand     = []

    # Loop over random runs
    for _ in trange(n_random_sample_size, leave=False, disable=True) :

        B_rand = generate_random_sample(A_RMs_s, A_ij_s, conserve_RM_degree = conserve_RM_degree)

        # Get the associated matrices from the random network
        A_ij_rand, A_RMs_rand, _ = compute_network(B_rand)   # Diagonal counts occurances of RM system. Off diagnoal counts combinations

        # Get metrics of random network
        avg_overlap_rand.append(np.mean(compute_overlap_measure(A_RMs_rand, A_ij_rand)))
        avg_RM_per_B_rand.append(A_ij_rand.sum(axis=1).mean())
        avg_degree_rand.append(A_RMs_rand.sum(axis=1).mean())
        avg_conn_rand.append(np.clip(A_RMs_rand, 0, 1).sum(axis=1).mean())
        f_unique_rand.append(compute_unique_RM_fraction(B_rand))
        avg_o_ij_rand.append(np.mean(compute_relative_overlaps(A_ij_rand, B_rand)))


    return avg_overlap_rand, avg_RM_per_B_rand, avg_degree_rand, avg_conn_rand, f_unique_rand, avg_o_ij_rand

def plot_network(A_ij, A_RMs, A_Bs, title, network_type='RM', pos=None, color=None, no_labels=False, base_width=1, scaling=1, ax=None, **kwargs) :

    if network_type == 'RM' :
        A = A_RMs
        labels = np.argwhere(A_ij.sum(axis=0) > 0).flatten()
        if color is None :
            color = 'r'
    elif network_type == 'B' :
        A = A_Bs
        labels = np.argwhere(A.sum(axis=0) > 0).flatten()
        if color is None :
            color = 'b'

    mar = 0.4
    scale = 1.2

    if ax is None :
        fig = plt.figure(figsize=(4, (2+mar)*scale / (2*scale) * 4))
        ax = plt.gca()
    else :
        fig = ax.figure

    # Adjacency matrix to graph
    G = nx.convert_matrix.from_numpy_matrix(A)
    degree = np.array([d[1] for d in G.degree])

    # Adjust the edge weights
    weights = [G.edges[e]['weight'] for e in G.edges]
    max_weight = np.max(weights)

    # Set the labels
    labels = [str(lb).center(3) for lb in labels]

    # Remove labels with degree less than trehshold
    labels = [lb if d > np.quantile(degree, 0.7) else ' ' for lb, d in zip(labels, degree)]

    # Generate node locations
    if pos is None :
        pos_nodes=nx.circular_layout(G)
        node_size = degree
    else :
        pos_nodes = pos
        node_size = 5

    # Generate label locations
    if pos is None :
        pos_labels = pos_nodes.copy()
        for key, val in pos_labels.items() :
            pos_labels[key] = [1.08 * v for v in val]
    else :
        pos_labels = pos_nodes

    # Draw labels
    if not no_labels :
        nx.draw_networkx_labels(G, pos=pos_labels, labels=dict(zip(np.arange(len(labels)), labels)), font_size=10*scaling, ax=ax)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos=pos_nodes, node_size=node_size*scaling**2, node_color=color, ax=ax)

    # Draw edges
    for weight in np.unique(weights) :
        edges = [e for e in G.edges if G.edges[e]['weight'] == weight]
        nx.draw_networkx_edges(G, pos=pos_nodes, edgelist=edges, width=base_width *  weight / np.sqrt(max_weight)*scaling, ax=ax)


    # Add meta data
    c_RMs = np.diag(A_RMs)                      # Compute overlap
    overlap_s = c_RMs / (np.size(A_ij, 0) - 1)  # p_i = # occurances / # number of other strains

    if pos is None :
        fontsize = 18 * scaling
        ax.text(-scale*0.9, 1.5, title, fontsize=1.5*fontsize, fontweight='bold')

        ax.text(-scale*0.9, 1.3, f'# RM = {len(A_RMs)}',     fontsize=fontsize)
        ax.text(-scale*0.9, 1.2, f'# strains = {len(A_ij)}', fontsize=fontsize)

        ax.text(0, 1.3, f'$\langle$# RM$\\rangle$ = {A_ij.sum(axis=1).mean().round(2)}', fontsize=fontsize)
        ax.text(0, 1.2, f'$\langle$p$\\rangle$ = {overlap_s.mean().round(3)}',           fontsize=fontsize)

        ax.set_xlim(-scale, scale)
        ax.set_ylim(-scale, (1+mar)*scale)

        ax.set_xlim(-scale, scale)
        ax.set_ylim(-scale, (1+mar)*scale)
    else :
        ax.set(**kwargs)


    return fig

def plot_bipartite_network(B, A_ij, A_RMs, title, scaling=1, ax=None) :

    mar = 0.4
    scale = 1.2

    if ax is None :
        fig = plt.figure(figsize=(4, (2+mar)*scale / (2*scale) * 4))
        ax = plt.gca()
    else :
        fig = ax.figure

    # Prepare the network
    G = nx.Graph()
    strain_ids = np.arange(len(B))
    RMs_ids    = [str(r) for r in np.unique(flatten_list(B))]

    # Add nodes with the node attribute "bipartite"
    G.add_nodes_from(strain_ids, bipartite=0)  # Stain nodes
    G.add_nodes_from(RMs_ids,    bipartite=1)  # RM nodes

    for strain_id, b in enumerate(B) :
        for r in b :
            G.add_edge(strain_id, str(r))


    # Plot the graph
    pos = nx.spring_layout(G, iterations=100)

    # Draw strain nodes
    colors = ['b' if G.nodes[n]['bipartite'] == 0 else 'r' for n in G.nodes()]
    nx.draw_networkx_nodes(G, ax=ax, pos=pos, node_color=colors, node_size=30*scaling**2)

    # Draw edges
    nx.draw_networkx_edges(G, ax=ax, pos=pos, width=scaling)

    fontsize = 18 * scaling
    dy = 0.1*(1 + scaling)
    ax.text(-scale*0.9, 1+3.3*dy, title, fontsize=1.5*fontsize, fontweight='bold')
    ax.text(-scale*0.9, 1+2*dy, f'# RM = {len(A_RMs)}',     fontsize=fontsize)
    ax.text(0,          1+2*dy, f'# strains = {len(A_ij)}', fontsize=fontsize)

    ax.text(-scale*0.9, 1+1*dy, f'$\langle$# RM$\\rangle$ = {A_ij.sum(axis=1).mean().round(2)}', fontsize=fontsize)
    ax.text(0,          1+1*dy, f'$f^u$ = {compute_unique_RM_fraction(B):.2f}', fontsize=fontsize)


    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, (1+mar)*scale)

    return fig

def plot_network_histogram(data, bins=None, title=None, normalized=False, clip=None, ax=None, xlabel='# RM per bacterium', ylabel='Count', labelpad = None, color=None, label='_nolegend_', **kwargs) :

    mar = 0.4
    scale = 1.2

    if ax is None :
        fig = plt.figure(figsize=(4, (2+mar)*scale / (2*scale) * 4))
        ax = plt.gca()
    else :
        fig = ax.figure

    if clip is not None :
        data = np.clip(data, 0, clip)

    if bins is None :
        bin_edges = np.arange(-0.5, np.max(data)+1.5)
    else :
        bin_edges = bins

    if normalized :
        weights = np.full(len(data), fill_value=1/(len(bin_edges)-1)/len(data))
    else :
        weights = None

    fill = True
    if color is None :
        color = plt.rcParams['patch.facecolor']
    elif color == 'none' :
        color = None
        fill = False

    ax.hist(data, bin_edges, label=label, color=color, linewidth=plt.rcParams['axes.linewidth'], edgecolor='k', weights=weights, fill=fill, **kwargs)

    if title is not None :
        ax.set_title(title, fontweight='bold')

    if xlabel is not None :
        ax.set_xlabel(xlabel)

    if ylabel is not None :
        ax.set_ylabel(ylabel, labelpad = labelpad)

    if bin_edges[-1] < 15 :
        s = 1
    else :
        s = 2

    ax.set_xticks(np.arange(bin_edges[-1]+1, step = s))
    ax.set_xlim(bin_edges[0], bin_edges[-1])

    return fig

def plot_summery_figure(axes, names, data, disable_significance=None, plot_names=True, sort=True) :

    if sort :
        I = np.argsort(names)
    else :
        I = np.arange(len(names))

    if disable_significance is None :
        disable_significance = [False]*len(axes)

    for strain_id, samples in enumerate(data) :

        color = plt.cm.tab10(I[strain_id] % 10)

        p_vals = [1.0, 0.05, 0.01, 0.001]
        #p_vals = [1.0, 0.317, 0.046, 0.003]
        for ax_id, (sample, ax) in enumerate(zip(samples, axes)) :

            err = np.std(sample)

            if err > 0 :

                # Compute z scores
                z = np.mean(sample) / err

                # Plot the intervals
                ax.errorbar(I[strain_id], np.mean(sample),
                            xerr = 0,
                            yerr = err,
                            marker = '.',
                            capsize = 5, markeredgewidth=2,
                            color = color)

                # Add significance annotation
                if not disable_significance[ax_id] :

                    if err > 0 :
                        # Create significance label
                        p = scipy.special.erfc(abs(z) / np.sqrt(2))
                        sig = '∗' * (np.argwhere(p < p_vals)[-1][-1])
                    else :
                        sig = '+'

                    if z > 0 :
                        sig = '  ' + sig
                        ha = 'center'
                        va = 'bottom'

                    else :
                        sig += '  '
                        ha = 'center'
                        va = 'top'

                    ax.text(I[strain_id], np.mean(sample) + np.sign(z)*np.std(sample), sig, verticalalignment=va, horizontalalignment=ha, fontstretch='condensed', fontsize='x-small', rotation='vertical')

            else :

                # Plot the central value intervals
                ax.scatter(I[strain_id], np.mean(sample),
                            marker = 'x',
                            color = color)

    for ax_id, ax in enumerate(axes) :

        if not disable_significance[ax_id] :
            ax.plot([-1, len(names)+1], [0, 0], 'k:', linewidth=2, zorder=-1)
        ax.set_xlim(-0.5, len(names)-0.5)
        ax.set_xticks(np.arange(len(names)))

        if plot_names and ax_id == 0 :
            ax.xaxis.tick_top()
            ax.set_xticklabels([names[i] for i  in I], rotation='vertical', ha='center', fontsize='small')

        else :
            ax.tick_params(labelbottom=False, bottom=False)

def plot_significant_metrics(x, y, dx, dy, colorings, suffixes, xlabel, ylabel, path, prefix = 'summery_significant_metrics') :

    # Set the cutoff
    split = 0.103

    print(x.min())
    print(x.max())
    print(y.min())
    print(y.max())

    k = 5
    for suffix, coloring in zip(suffixes, colorings) :
        fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, facecolor='w', gridspec_kw={'width_ratios': [k, 1]}, figsize=(10, 4))
        # https://stackoverflow.com/questions/32185411/break-in-x-axis-of-matplotlib

        fig.subplots_adjust(wspace=0.05)

        for ax_id, ax in enumerate(axes) :
            if ax_id == 0 :
                s = x <= split
            else :
                s = x > split

            ax.errorbar(x[s], y[s], xerr=dx[s], yerr=dy[s], ls='none', zorder=0, ecolor='k', elinewidth=0.5)

            I = np.isfinite(coloring)
            if np.any(coloring[I] > 0) and np.any(coloring[I] < 0) :    # Use symmetric coloring
                cm_p = np.max(np.abs(coloring[I]))
                cm_m = - cm_p
                cmap = 'coolwarm'
            elif np.any(np.array(coloring[I]) < 0) :
                cm_p = np.max(coloring[I])
                cm_m = np.min(coloring[I])
                cmap = 'Reds'
            else :
                cm_p = np.max(coloring[I])
                cm_m = 0
                cmap = 'Reds'

            h = ax.scatter(x[np.logical_and(s,  I)], y[np.logical_and(s,  I)], c=coloring[np.logical_and(s,  I)], vmin=cm_m, vmax=cm_p, s=20, marker='o', linewidths=0.5, edgecolors='k', cmap=cmap)
            ax.scatter(    x[np.logical_and(s, ~I)], y[np.logical_and(s, ~I)], c='w',                             vmin=cm_m, vmax=cm_p, s=20, marker='o', linewidths=0.5, edgecolors='k', cmap=cmap)


            ax.plot([0, 0],      [-0.1, 0.4], ':k', zorder=0, linewidth=1)
            ax.plot([-0.1, 0.3], [0, 0],      ':k', zorder=0, linewidth=1)

            ax.set_yticks(np.arange(-1, 1, 0.1))

            ax.set_ylim(-0.1, 0.4)

            cax = inset_axes(axes[1], width='200%', height='5%', loc='upper right')
            fig.colorbar(h, cax=cax, orientation='horizontal')

        axes[0].set_xticks(np.arange(-1, 1, 0.05))
        axes[1].set_xticks(np.arange(-1, 1, 0.1))

        axes[0].set_xlim(-0.035, split)
        axes[1].set_xlim( split,   0.3)

        axes[0].spines['right'].set_visible(False)
        axes[1].spines['left'].set_visible(False)
        axes[0].tick_params(labelright='off')

        axes[0].yaxis.tick_left()
        axes[1].yaxis.tick_right()


        d = 0.015 # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=axes[0].transAxes, color='k', clip_on=False, linewidth=1)
        axes[0].plot((1-d/2, 1+d/2),  (-d,  +d), **kwargs)
        axes[0].plot((1-d/2, 1+d/2), (1-d, 1+d), **kwargs)

        kwargs.update(transform=axes[1].transAxes)  # switch to the bottom axes
        axes[1].plot((-k*d/2, +k*d/2),  (-d, +d),  **kwargs)
        axes[1].plot((-k*d/2, +k*d/2), (1-d, 1+d), **kwargs)

        axes[0].set_xlabel(xlabel, x = 0.65)
        axes[0].set_ylabel(ylabel)


        fig_name = prefix + suffix

        fig.savefig(os.path.join(path, fig_name+'.png'), bbox_inches = 'tight')

def add_figure_labels(labels, axes, dx=-0.05, dy=0.05) :

    # Add the figure labels
    for ax, label in zip(axes, labels) :
        pos = ax.get_position()
        ax_label = inset_axes(ax, width=0, height=0,
                        bbox_to_anchor=(pos.x0+dx, pos.y1+dy),
                        bbox_transform=ax.figure.transFigure,
                        loc='upper left')

        ax_label.set_axis_off()
        ax_label.text(1, 1, label, label_args())

def load_sequence_data() :

    # Load data
    data_path = os.path.join(get_path('data'), 'Fig2')
    A_ij = np.loadtxt(os.path.join(data_path, 'A_ij.csv'), dtype=int, delimiter=',')
    data_taxonomy = pd.read_csv(os.path.join(data_path, 'data_taxonomy.csv'), index_col=0)
    data_motifs = pd.read_csv(os.path.join(data_path, 'motifs.csv'))

    # Group by 'genus'
    s = data_taxonomy.genus_id.values
    genus_id_to_name = pd.Series(index=s, data=data_taxonomy.genus.values).drop_duplicates()

    # Ensure the data is sorted
    sort_index = np.argsort(s)
    s = s[sort_index]
    A_ij = A_ij[sort_index, :]

    # Perform sanity checks
    for genus_id in data_taxonomy.genus_id.unique() :
        subset = data_taxonomy.loc[data_taxonomy['genus_id'] == genus_id]

        if len(subset.subgroup.unique()) > 1 :
            raise warnings.warn(f'The genus "{subset.genus.unique()}" has several sup-groups')

        if len(subset.group.unique()) > 1 :
            raise warnings.warn(f'The genus "{subset.genus.unique()}" has several groups')

    return A_ij, s, data_taxonomy, data_motifs, genus_id_to_name

def iter_sequence_data(A_ij, s, data_taxonomy, treshold = 15, tqdm_disable=False) :

    # Loop over subsets
    for strain_id in tqdm(np.unique(s), disable=tqdm_disable) :

        # Locate the start and size of the subset
        idx = np.argmax(s == strain_id) # Start
        ns  = np.sum(s == strain_id)    # Size

        # Extract the subset
        A_ij_s = A_ij[idx:idx+(ns-1), :]

        # Filter the subset for unique combinations
        A_ij_s = np.unique(A_ij_s, axis=0)

        # Label the (sub) adjacency matrix
        A_ij_s_tmp = A_ij_s * np.repeat(np.arange(start=1, stop=np.size(A_ij_s, 1)+1)[np.newaxis, :], len(A_ij_s), axis=0)

        # Convert to list of RM format
        Bs = []
        for a_ij in A_ij_s_tmp :
            Bs.append(tuple(np.where(a_ij > 0)[0])) # Here, use tuple instead of lists


        # Get the name of the subset
        genus = data_taxonomy.genus[data_taxonomy.genus_id == strain_id][0]

        # Filter out small sample sets
        if np.size(A_ij_s, 0) < treshold :
            continue


        yield strain_id, A_ij_s, Bs, genus

def iter_simulation_data(n_sims = 5, n_seeds = 6, n_its = 6, tqdm_disable = False) :

    for simulation_id, seed_id, log_iterations in tqdm(product(range(n_sims), range(n_seeds), range(3, n_its+1)), total=n_sims*n_seeds*(n_its-2), disable=tqdm_disable) :
        yield 50*2**simulation_id, seed_id, log_iterations

def load_reference_data() :

    # Load data
    data_path = os.path.join(get_path('data'), 'FigS2_S3')

    lname = os.path.join(data_path, 'data_Roer_Fullmer.mat')
    if os.path.exists(lname) :
        data = scipy.io.loadmat(lname)
    else :
        raise ValueError('Data not found!')

    A_ij = data['A_ij'].astype(int)
    s    = data['s'] - 1
    genus_id_to_name = pd.Series({ 0 : 'Halobacteria (Roer et al.)', 1 : 'Salmonella (Fullmer et al.)'})

    return A_ij, s.flatten(), genus_id_to_name

def iter_reference_data(A_ij, s, tqdm_disable=False) :

    # Loop over subsets
    for strain_id in tqdm(np.unique(s), disable=tqdm_disable) :

        # Locate the start and size of the subset
        idx = np.argmax(s == strain_id) # Start
        ns  = np.sum(s == strain_id)    # Size

        # Extract the subset
        A_ij_s = A_ij[idx:idx+(ns-1), :]

        # Filter the subset for unique combinations
        A_ij_s = np.unique(A_ij_s, axis=0)

        # Label the (sub) adjacency matrix
        A_ij_s_tmp = A_ij_s * np.repeat(np.arange(start=1, stop=np.size(A_ij_s, 1)+1)[np.newaxis, :], len(A_ij_s), axis=0)

        # Convert to list of RM format
        Bs = []
        for a_ij in A_ij_s_tmp :
            Bs.append(tuple(np.where(a_ij > 0)[0])) # Here, use tuple instead of lists

        yield strain_id, A_ij_s, Bs

def analyze_sequence_data(n_random_sample_size = 100, conserve_RM_degree = False) :

    # Load the sequence data
    A_ij, s, data_taxonomy, _, genus_id_to_name = load_sequence_data()

    # Prepare outputs
    diff_random_overlap = []
    diff_random_unique  = []
    diff_random_o_ij    = []

    RMs_abundence       = []

    hist_avg_overlap    = []
    hist_avg_RM_per_B   = []
    hist_imag_eig       = []

    names = []

    # Loop over subsets to data figures
    for strain_id, _, Bs, _ in iter_sequence_data(A_ij, s, data_taxonomy) :

        # Get the matrices for the networks
        _, _, A_Bs_s = compute_network(Bs)

        # Get the metrics from the sequence data
        avg_overlap, RMs_per_B, _, _, f_unique, avg_o_ij = get_metrics_of_real_network(Bs)

        # Get the metrics for the random comparrison
        avg_overlap_rand, _, _, _, f_unique_rand, avg_o_ij_rand = get_metrics_of_random_network(Bs, n_random_sample_size, conserve_RM_degree = conserve_RM_degree)

        # Store metric of difference between random distrubution and the sample
        diff_random_overlap.append(avg_overlap  - np.array(avg_overlap_rand) )
        diff_random_unique.append( f_unique     - np.array(f_unique_rand)  )
        diff_random_o_ij.append(   avg_o_ij     - np.array(avg_o_ij_rand)    )

        RMs_abundence.append(RMs_per_B)

        # Store average measures
        E, _ = scipy.linalg.eig(A_Bs_s)
        E = np.abs(np.imag(E))

        hist_imag_eig.append(np.sum(E[E>0]) / len(A_Bs_s))
        hist_avg_RM_per_B.append(RMs_per_B.mean())
        hist_avg_overlap.append(avg_overlap)

        # Store the names
        names.append(genus_id_to_name[strain_id])

    return diff_random_overlap, diff_random_unique, diff_random_o_ij, RMs_abundence, hist_avg_RM_per_B, hist_avg_overlap, hist_imag_eig, names

def analyze_simulation_data(n_random_sample_size = 10, conserve_RM_degree = False) :

    # Allocate arrays
    names = []
    iters = []

    # Define the nunmber of random samples to generate
    diff_random_overlap = []
    diff_random_unique  = []
    diff_random_o_ij    = []

    RMs_abundence       = []

    hist_avg_overlap    = []
    hist_avg_RM_per_B   = []

    # Loop over subsets to generate figures
    for K, seed_id, log_iterations in tqdm(iter_simulation_data()) :

        # Set settings for the run
        names.append(f'RM : {K}, it : 1e{log_iterations}, seed : {seed_id}')

        # Load the data
        lname = os.path.join(get_path('data'), 'Fig6', f'RM_{K}_seed_{seed_id}.mat')
        if os.path.exists(lname) :
            data = scipy.io.loadmat(lname)

        else :
            raise ValueError(f'{lname} is missing!')


        # Get data from the run
        Bs = [b[0, :].tolist() for b in data['B_samples'][0, :][log_iterations].flatten().tolist()] # wtf scipy?

        # Get the metrics from the simulation
        avg_overlap, RMs_per_B, _, _, f_unique, avg_o_ij = get_metrics_of_real_network(Bs)

        # Get the metrics for the random comparrison
        avg_overlap_rand, _, _, _, f_unique_rand, avg_o_ij_rand = get_metrics_of_random_network(Bs, n_random_sample_size, conserve_RM_degree=conserve_RM_degree)

        # Store metric of difference between random distrubution and the sample
        diff_random_overlap.append(avg_overlap- np.array(avg_overlap_rand) )
        diff_random_unique.append( f_unique   - np.array(f_unique_rand)  )
        diff_random_o_ij.append(   avg_o_ij   - np.array(avg_o_ij_rand)    )

        RMs_abundence.append(RMs_per_B)

        # Store the measures
        hist_avg_RM_per_B.append(RMs_per_B.mean())
        hist_avg_overlap.append(avg_overlap)

        # Store iteration id
        iters.append(log_iterations)

    return diff_random_overlap, diff_random_unique, diff_random_o_ij, RMs_abundence, hist_avg_RM_per_B, hist_avg_overlap, names, iters

def make_dirs(path) :
    # Prepare the figure folder
    if not os.path.exists(path) :
        os.makedirs(path)

def B_iter(B) :
    for b in B :
        yield np.array(b)

def get_RMs(B, P=[]) :

    # Special case without RM systems
    if B == [] :
        return

    RMs = np.max(flatten_list(B+P)) + 1
    return np.arange(RMs)

def flatten_list(lst) :
    return reduce(operator.iconcat, lst, [])

def get_path(path_type='figures') :

    if path_type == 'figures':
        path = 'figures'

    elif path_type == 'data':
        path = os.path.join('code', 'data')

    else :
        raise ValueError('path_type not known!')

    make_dirs(path)

    return path

def pickle_write(fname, *args) :

    with open(fname, 'wb') as outfile:
        for x in args :
            pickle.dump(x,    outfile)

def pickle_read(fname) :

    with open(fname, 'rb') as outfile:

        return list(pickle_reader(outfile))

def pickle_reader(f) :

    while True:
        try:
            yield pickle.load(f)
        except EOFError:
            break
