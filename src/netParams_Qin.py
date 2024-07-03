from netpyne import specs
import numpy as np

from __main__ import cfg

netParams = specs.NetParams()

# ------------------------------------------------ CELL TYPES ------------------------------------------------
from neuron import h
resultCode = h.load_file('cells/SThprotocell.hoc')
if resultCode:
    protoCell = h.SThproto()
else:
    raise RuntimeError("Couldn't load prototype STN cell from `cells/STHprotocell.hoc")

pvp = netParams.importCellParams(label='PVP_cell', fileName='cells/SThprotocell.hoc', cellName='SThcell', cellArgs=[0, protoCell]) # 0 in cellArgs is just placeholder
pvn = netParams.importCellParams(label='PVN_cell', fileName='cells/SThprotocell.hoc', cellName='SThcell', cellArgs=[0, protoCell])
# so far they are identical. Difference in some param values will be introduced after network is created in init.py

pvp.secs.soma['threshold'] = -30
pvn.secs.soma['threshold'] = -30


# ------------------------------------------------ NETWORK & POPULATIONS ------------------------------------------------
popOrder = 10 # we can start with lower value for speed, and then increase gradually
popScale = 0.2 # the percentage of PVP neurons; 1-popScale is the percetage of PVN neurons; cholinergic interneurons neglected for now

def generate_locs(num,pv_percent):
    randGen = np.random.default_rng(seed=cfg.seeds['loc'])
    #rawLocsX = randGen.uniform(0, netParams.sizeX, int(num))
    #rawLocsY = randGen.uniform(0, netParams.sizeY, int(num))

    # now we can apply any transform to the raw locations (maybe multiply with some matrix, or use any custom rule)
    # ...
    # ...
    # Total number of PV neurons
    total_points = num*pv_percent

    # Define parameters for the skew normal distribution
    n_grids = 10
    peak_grid = 7
    mean = (peak_grid - 1) / n_grids
    std_dev = 0.4  # Standard deviation
    a = .2  # Skewness parameter: positive for right skew, negative for left skew

    # Generate the skew normal distribution values
    grid_centers = np.linspace(0, 1, n_grids+1)
    skewed_dist = skewnorm.pdf(grid_centers, a, loc=mean, scale=std_dev)
    skewed_dist /= skewed_dist.sum()  # Normalize to make it a probability distribution

    # Allocate the number of data points to each grid based on the skewed distribution
    pvp_points_per_grid = (skewed_dist * total_points).astype(int)
    points_per_grid = num/n_grids
    pvn_points_per_grid = np.array([np.int(points_per_grid-i) for i in pvp_points_per_grid])

    pvp_x = np.zeros(pvp_points_per_grid.sum(),)
    pvp_y = np.zeros(pvp_points_per_grid.sum(),)
    pvn_x = np.zeros(pvn_points_per_grid.sum(),)
    pvn_y = np.zeros(pvn_points_per_grid.sum(),)

    start_idx = 0
    for i in range(n_grids):
        end_idx = start_idx + pvp_points_per_grid[i]
        if end_idx > start_idx:  # Avoid empty ranges
            pvp_x[start_idx:end_idx] = randGen(grid_centers[i], grid_centers[i+1], size=(end_idx - start_idx))
            pvp_y[start_idx:end_idx] = randGen(grid_centers[i], grid_centers[i+1], size=(end_idx - start_idx))
        start_idx = end_idx
        
    start_idx = 0
    for i in range(n_grids):
        end_idx = start_idx + pvn_points_per_grid[i]
        if end_idx > start_idx:  # Avoid empty ranges
            pvn_x[start_idx:end_idx] = randGen(grid_centers[i], grid_centers[i+1], size=(end_idx - start_idx))
            pvn_y[start_idx:end_idx] = randGen(grid_centers[i], grid_centers[i+1], size=(end_idx - start_idx))
        start_idx = end_idx

    return [{'x': x, 'y': y} for x,y in zip(pvp_x, pvp_y)],[{'x': x, 'y': y} for x,y in zip(pvn_x, pvn_y)] # alternatively, `x_norm` and/or `y_norm` can be used, presuming that values are within the range [0,1]

#def generate_locs_PVN(num):
    # TODO: jjust using placeholder, replace with proper logic
    #return generate_locs_PVP(num)

pvpLocs,pvnLocs = generate_locs(popOrder,popScale)

netParams.popParams['PVP_pop'] = {'cellType': 'PVP_cell', 'cellsList': pvpLocs}
netParams.popParams['PVN_pop'] = {'cellType': 'PVN_cell','cellsList': pvnLocs}

# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# to be defined..

# ------------------------------------------------ STIMULATION ------------------------------------------------

cfg.IClamp_amp = -0.16
netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': cfg.IClamp_amp}
netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'PVP_cell'}, 'sec': 'soma', 'loc': 0.5}