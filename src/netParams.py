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
popScales = 1.5, 0.3

def generate_locs_PVP(num):
    randGen = np.random.default_rng(seed=cfg.seeds['loc'])
    rawLocsX = randGen.uniform(0, netParams.sizeX, int(num))
    rawLocsY = randGen.uniform(0, netParams.sizeY, int(num))

    # now we can apply any transform to the raw locations (maybe multiply with some matrix, or use any custom rule)
    # ...
    # ...

    return [{'x': x, 'y': y} for x,y in zip(rawLocsX, rawLocsY)] # alternatively, `x_norm` and/or `y_norm` can be used, presuming that values are within the range [0,1]

def generate_locs_PVN(num):
    # TODO: jjust using placeholder, replace with proper logic
    return generate_locs_PVP(num)

pvpLocs = generate_locs_PVP(popOrder * popScales[0])
pvnLocs = generate_locs_PVN(popOrder * popScales[1])

netParams.popParams['PVP_pop'] = {'cellType': 'PVP_cell', 'cellsList': pvpLocs}
netParams.popParams['PVN_pop'] = {'cellType': 'PVN_cell','cellsList': pvnLocs}

# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# to be defined..

# ------------------------------------------------ STIMULATION ------------------------------------------------

cfg.IClamp_amp = -0.16
netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': cfg.IClamp_amp}
netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'PVP_cell'}, 'sec': 'soma', 'loc': 0.5}