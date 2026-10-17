from netpyne import specs
import numpy as np
from scipy.stats import skewnorm

from __main__ import cfg

netParams = specs.NetParams()
netParams.sizeX = cfg.sizeX
netParams.sizeY = cfg.sizeY
netParams.sizeZ = cfg.sizeZ
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
PVPproportion = 0.2 # the percentage of PVP neurons; 1-popScale is the percetage of PVN neurons; cholinergic interneurons neglected for now

from tools import generate_locs
pvpLocs, pvnLocs = generate_locs(popOrder, PVPproportion, seed=cfg.seeds['loc'])

# STN pops
netParams.popParams['PVP_pop'] = {'cellType': 'PVP_cell', 'cellsList': pvpLocs}
netParams.popParams['PVN_pop'] = {'cellType': 'PVN_cell', 'cellsList': pvnLocs}

# External inputs
# Cortex: arrival events from src/cortical_arrivals.py (settings: cfg.ctx_*), one VecStim per STN cell, connected 1:1
# GPe: closed-loop model from src/gpe_feedback.py (settings: cfg.gpe_*), driven during the simulation from init.py
from stn_inputs import add_cortical_input
add_cortical_input(netParams, cfg, {'PVP_pop': pvpLocs, 'PVN_pop': pvnLocs})

# ------------------------------------------------ SYNAPTIC MECHS ------------------------------------------------
netParams.synMechParams['glut'] = {'mod': 'Exp2Syn', 'e': 0, 'tau1': 1.0, 'tau2': 5.0}
netParams.synMechParams['GABA'] = {'mod': 'Exp2Syn', 'e': -85, 'tau1': 2, 'tau2': 15.0}


# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# Cortex->STN connections are added by add_cortical_input() above; GPe->STN synapses are created in init.py.
# Cortical topography (cfg.ctx_topography) scales each STN cell's cortical arrival rate by its position along the STN axis.

# # STN->STN - free parameter
netParams.probLengthConst = 300

netParams.connParams['PVP->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVP_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5} # 'connList': [[0,1], [0,2]]

netParams.connParams['PVN->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVN_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5} # exp(-dist_3D/probLengthConst) * 0.25'

# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# to be defined..

# ------------------------------------------------ STIMULATION ------------------------------------------------

#cfg.IClamp_amp = -0.16
#netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': cfg.IClamp_amp}
#netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'PVP_cell'}, 'sec': 'soma', 'loc': 0.5}