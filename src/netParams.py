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
PVPproportion = 0.2 # the percentage of PVP neurons; 1-popScale is the percetage of PVN neurons; cholinergic interneurons neglected for now

from netParams_Qin import generate_locs
pvpLocs, pvnLocs = generate_locs(popOrder, PVPproportion)

# STN pops
netParams.popParams['PVP_pop'] = {'cellType': 'PVP_cell', 'cellsList': pvpLocs}
netParams.popParams['PVN_pop'] = {'cellType': 'PVN_cell', 'cellsList': pvnLocs}

# External pops:
# Noisy excitation from cortex
netParams.popParams['Ctx_pop'] = {'cellModel': 'NetStim', 'numCells': int(popOrder*2), 'start': '300 + uniform(-2, 2)', 'interval': 'normal(35,5)', 'number': 1e9, 'noise': 0.3}

# rhythmic inhiibiitiion from GPe
spikePat = {'type': 'rhythmic', 'start': -1, 'startMin':340, 'startMax': 360, 'stop': 10000, 'freq': 3, 'freqStd': 80, 'distribution': 'normal', 'eventsPerCycle': 1, 'repeats': 45} # 'freqStd': 50
netParams.popParams['GPe_pop'] = {'cellModel': 'VecStim', 'numCells': popOrder*2, 'spikePattern': spikePat}

# ------------------------------------------------ SYNAPTIC MECHS ------------------------------------------------
netParams.synMechParams['glut'] = {'mod': 'Exp2Syn', 'e': 0, 'tau1': 1.0, 'tau2': 5.0}
netParams.synMechParams['GABA'] = {'mod': 'Exp2Syn', 'e': -85, 'tau1': 2, 'tau2': 15.0}


# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# # noisy input (poisson)
# strengthFun = ''
# netParams.connParams['Ctx->PVP'] = {'synMech': 'exc', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': 'PVP_pop'}}
# netParams.connParams['Ctx->PVN'] = {'synMech': 'exc', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'cellType': 'PVN_pop'}}

E_I_balance_PVP = 1.3
E_I_balance_PVN = 0.8

netParams.connParams['Ctx->STN'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 0.75, 'weight': 0.0002, 'delay': 1}

netParams.connParams['GPe->PVP'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVP_pop'}, 'probability':  0.75, 'weight': E_I_balance_PVP * 0.0002, 'delay': 1}
netParams.connParams['GPe->PVN'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability':  0.75, 'weight': E_I_balance_PVN * 0.0002, 'delay': 1}

# # STN->STN - free parameter
netParams.probLengthConst = 1

netParams.connParams['PVP->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVP_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5} # 'connList': [[0,1], [0,2]]

netParams.connParams['PVN->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVN_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5}

# ------------------------------------------------ STIMULATION ------------------------------------------------

# cfg.IClamp_amp = -0.16
# netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 300, 'dur': 40000, 'amp': cfg.IClamp_amp}
# netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'PVP_cell'}, 'sec': 'soma', 'loc': 0.5}