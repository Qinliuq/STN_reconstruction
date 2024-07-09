from netpyne import specs
import numpy as np
from scipy.stats import skewnorm

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
popOrder = 30 # we can start with lower value for speed, and then increase gradually
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
    n_grids = 12
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
    error = total_points - pvp_points_per_grid.sum()
    pvp_points_per_grid[6]+= error

    temp = num/(n_grids-2)
    points_per_grid = [temp*.3,temp*.7,temp,temp,temp,temp,temp,temp,temp,temp,temp*.7,temp*.3]
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

    return [{'x': x, 'y': y} for x,y in zip(pvp_x, pvp_y)],[{'x': x, 'y': y} for x,y in zip(pvn_x, pvn_y)] # alternatively, `x_norm` and/or `y_norm` can be used, presuming that values are within the range [0,1]

#def generate_locs_PVN(num):
    # TODO: jjust using placeholder, replace with proper logic
    #return generate_locs_PVP(num)

pvpLocs,pvnLocs = generate_locs(popOrder,popScale)

# STN pops
netParams.popParams['PVP_pop'] = {'cellType': 'PVP_cell', 'cellsList': pvpLocs}
netParams.popParams['PVN_pop'] = {'cellType': 'PVN_cell', 'cellsList': pvnLocs}

# External pops:
# Noisy excitation from cortex
netParams.popParams['Ctx_pop'] = {'cellModel': 'NetStim', 'numCells': int(popOrder*2), 'start': '300 + uniform(-2, 2)', 'interval': 'normal(35,5)', 'number': 1e9, 'noise': 0.3}

# rhythmic inhiibiitiion from GPe
spikePat = {'type': 'rhythmic', 'start': -1, 'startMin':340, 'startMax': 360, 'stop': 10000, 'freq': 15, 'freqStd': 1, 'distribution': 'normal', 'eventsPerCycle': 1, 'repeats': 45} # 'freqStd': 50
netParams.popParams['GPe_pop'] = {'cellModel': 'VecStim', 'numCells': popOrder, 'spikePattern': spikePat}

# healthy inhibition from GPe

# ------------------------------------------------ SYNAPTIC MECHS ------------------------------------------------
netParams.synMechParams['glut'] = {'mod': 'Exp2Syn', 'e': 0, 'tau1': 1.0, 'tau2': 5.0}
netParams.synMechParams['GABA'] = {'mod': 'Exp2Syn', 'e': -85, 'tau1': 2, 'tau2': 15.0}


# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# # noisy input (poisson)
# strengthFun = ''
# netParams.connParams['Ctx->PVP'] = {'synMech': 'exc', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': 'PVP_pop'}}
# netParams.connParams['Ctx->PVN'] = {'synMech': 'exc', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'cellType': 'PVN_pop'}}

I_ratio = 0.9 #PVP/PVN
E_ratio = 0.6 #PVP/PVN 

#netParams.connParams['Ctx->STN'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 0.75, 'weight': 0.0002, 'delay': 1}
netParams.connParams['Ctx->PVP'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': 'PVP_pop'}, 'probability':  0.5, 'weight': E_ratio * 0.0002, 'delay': 1}
netParams.connParams['Ctx->PVN'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability':  0.5, 'weight': 0.0002, 'delay': 1}

netParams.connParams['GPe->PVP'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVP_pop'}, 'probability':  0.5, 'weight': I_ratio * 0.00002, 'delay': 1}
netParams.connParams['GPe->PVN'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability':  0.5, 'weight': 0.00002, 'delay': 1}

# # STN->STN - free parameter
netParams.probLengthConst = .1

netParams.connParams['PVP->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVP_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5} # 'connList': [[0,1], [0,2]]

netParams.connParams['PVN->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVN_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5}

# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# to be defined..

# ------------------------------------------------ STIMULATION ------------------------------------------------

cfg.IClamp_amp = -0.16
netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': cfg.IClamp_amp}
netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'PVP_cell'}, 'sec': 'soma', 'loc': 0.5}