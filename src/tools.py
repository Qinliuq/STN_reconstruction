from neuron import h
import numpy as np
import matplotlib.pyplot as plt

def set_values_from_file(cell, conductance, modifier=None, soma_scale=1.0, dend_scale=1.0):
    # use `modifier` to replicate some of paper result.
    # Using `soma_scale` and/or `dend_scale` should theoretically have the same effect, although in practice it's not always the case

    modif = f'-{modifier}' if modifier else ''
    f = open(f'cells/sth-data/cell-{conductance}{modif}', 'r')
    totval=0
    numval = 0
    for line in f:
        tree, sec, loc, val = line.strip().split(' ')
        loc = float(loc)
        val = float(val)

        if tree == '-1': # 'soma'
            hObj = cell.secs.soma.hObj
            val *= soma_scale
        else:
            hObj = getattr(cell.secs, f'dend{tree}_{sec}').hObj
            val *= dend_scale
            totval += val
            numval += 1
        setattr(hObj(loc), conductance, val)
    # print(f"Cond {conductance} avg {totval/numval} (n. {numval})")

def apply_CSF_Beurrier():
    print("Setting in vitro parameters based on Beurrier et al (1999)")

    h.nai0_na_ion = 15
    h.nao0_na_ion = 150

    h.ki0_k_ion = 140
    h.ko0_k_ion = 3.6

    h.cai0_ca_ion = 1e-04
    h.cao0_ca_ion = 2.4

    # # This was in original model, but seems to be not relevant
    # h.cli0_cl_ion = 4
    # h.clo0_cl_ion = 135

def apply_CSF_Bevan():
    print("Setting in vitro parameters based on Bevan & Wilson (1999)")

    h.nai0_na_ion = 15
    h.nao0_na_ion = 128.5

    h.ki0_k_ion = 140
    h.ko0_k_ion = 2.5

    h.cai0_ca_ion = 1e-04
    h.cao0_ca_ion = 2.0

    # # This was in original model, but seems to be not relevant
    # cli0_cl_ion = 4
    # clo0_cl_ion = 132.5

def generate_locs(num, pv_percent, seed=0, plot=True):

    from scipy.stats import skewnorm

    randGen = np.random.default_rng(seed=seed)
    #rawLocsX = randGen.uniform(0, netParams.sizeX, int(num))
    #rawLocsY = randGen.uniform(0, netParams.sizeY, int(num))

    # now we can apply any transform to the raw locations (maybe multiply with some matrix, or use any custom rule)
    # ...
    # ...
    # Total number of PV neurons
    total_points = num*pv_percent

    # Define parameters for the skew normal distribution
    n_grids = 12
    peak_grid = 9 # 8
    mean = (peak_grid - 1) / n_grids
    std_dev = 0.15  # Standard deviation
    a = .2  # Skewness parameter: positive for right skew, negative for left skew

    # Generate the skew normal distribution values
    grid_centers = np.linspace(0, 1, n_grids+1)
    skewed_dist = skewnorm.pdf(grid_centers, a, loc=mean, scale=std_dev)
    skewed_dist /= skewed_dist.sum()  # Normalize to make it a probability distribution

    # Allocate the number of data points to each grid based on the skewed distribution
    pvp_points_per_grid = (skewed_dist * total_points).astype(int)
    error = total_points - pvp_points_per_grid.sum()
    pvp_points_per_grid[6]+= error

    temp = num/(n_grids-2)
    points_per_grid = [temp*.3,temp*.7,temp,temp,temp,temp,temp,temp,temp,temp,temp*.7,temp*.3] #a temporary hack!!
    pvn_points_per_grid = np.array([np.int(points_per_grid[i]-pvp_points_per_grid[i]) for i in range(n_grids)])

    pvp_x = np.zeros(pvp_points_per_grid.sum(),)
    pvp_y = np.zeros(pvp_points_per_grid.sum(),)
    pvn_x = np.zeros(pvn_points_per_grid.sum(),)
    pvn_y = np.zeros(pvn_points_per_grid.sum(),)

    start_idx = 0
    for i in range(n_grids):
        end_idx = start_idx + pvp_points_per_grid[i]
        if end_idx > start_idx:  # Avoid empty ranges
            pvp_x[start_idx:end_idx] = randGen.normal(loc=grid_centers[i], scale=0.05, size=(end_idx - start_idx))
            pvp_y[start_idx:end_idx] = randGen.normal(loc=grid_centers[i], scale=0.05, size=(end_idx - start_idx))
        start_idx = end_idx
        
    start_idx = 0
    for i in range(n_grids):
        end_idx = start_idx + pvn_points_per_grid[i]
        if end_idx > start_idx:  # Avoid empty ranges
            pvn_x[start_idx:end_idx] = randGen.normal(loc=grid_centers[i], scale=0.05, size=(end_idx - start_idx))
            pvn_y[start_idx:end_idx] = randGen.normal(loc=grid_centers[i], scale=0.05, size=(end_idx - start_idx))
        start_idx = end_idx

    pvp_z = randGen.uniform(0, 1, len(pvp_x))
    pvn_z = randGen.uniform(0, 1, len(pvn_x))

    if plot:
        plt.figure(dpi=200)
        plt.scatter(pvn_x, pvn_y, marker='.', label=f'PV- ({len(pvn_x)} cells)', color='tab:blue')
        plt.scatter(pvp_x, pvp_y, marker='.', label=f'PV+ ({len(pvp_x)} cells)', color='tab:green')
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'STN ({num} cells, PV+ fraction = {pv_percent}, seed = {seed})')
        plt.savefig(f'net2D-n{num}-frac{pv_percent}-seed{seed}.png')

    return [{'xnorm': x, 'ynorm': y, 'znorm': z} for x, y, z in zip(pvp_x, pvp_y, pvp_z)], \
           [{'xnorm': x, 'ynorm': y, 'znorm': z} for x, y, z in zip(pvn_x, pvn_y, pvn_z)]