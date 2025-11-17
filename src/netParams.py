from netpyne import specs
import numpy as np
from scipy.stats import skewnorm

from __main__ import cfg

netParams = specs.NetParams()
netParams.sizeX = 850 #8.5 * 1e3 # so that diagonal axis is 12
netParams.sizeY = 850 #8.5 * 1e3
netParams.sizeZ = 100 #3 * 1e3
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

# External pops:
# Noisy excitation from cortex

#netParams.popParams['Ctx_pop'] = {'cellModel': 'NetStim', 'numCells': int(popOrder*2), 'start': '300 + uniform(-2, 2)', 'interval': 'normal(35,5)', 'number': 1e9, 'noise': 0.3}
spikePat_c = {'type': 'rhythmic', 'start':-1,'startMin':1, 'startMax': 10, 'stop': 1400, 'freq': 3, 'freqStd': 90, 'distribution': 'normal','eventsPerCycle': 1, 'repeats': 20} # 'freqStd': 50
netParams.popParams['Ctx_pop'] = {'cellModel': 'VecStim', 'numCells': popOrder*2, 'spikePattern': spikePat_c}

# rhythmic inhibiitiion from GPe
spikePat = {'type': 'rhythmic', 'start': -1, 'startMin':1, 'startMax': 10, 'stop': 1300, 'freq': 15, 'freqStd': 10, 'distribution': 'normal', 'eventsPerCycle': 1, 'repeats': 5} # 'freqStd': 50
netParams.popParams['GPe_pop'] = {'cellModel': 'VecStim', 'numCells': popOrder, 'spikePattern': spikePat}

# healthy inhibition from GPe
#netParams.popParams['GPe_H_pop'] = {'cellModel': 'NetStim', 'numCells': int(popOrder*2), 'start': '10 + uniform(-2, 2)', 'interval': 'normal(12,5)', 'number': 1e9, 'noise': 1.}
spikePat_g = {'type': 'rhythmic', 'start':-1,'startMin':1, 'startMax': 10, 'stop': 1400, 'freq': 3, 'freqStd': 90, 'distribution': 'normal','eventsPerCycle': 1, 'repeats': 20} # 'freqStd': 50
netParams.popParams['GPe_H_pop'] = {'cellModel': 'VecStim', 'numCells': popOrder*2, 'spikePattern': spikePat_g}

# ------------------------------------------------ SYNAPTIC MECHS ------------------------------------------------
netParams.synMechParams['glut'] = {'mod': 'Exp2Syn', 'e': 0, 'tau1': 1.0, 'tau2': 5.0}
netParams.synMechParams['GABA'] = {'mod': 'Exp2Syn', 'e': -85, 'tau1': 2, 'tau2': 15.0}


# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# # noisy input (poisson)
# strengthFun = ''
# netParams.connParams['Ctx->PVP'] = {'synMech': 'exc', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': 'PVP_pop'}}
# netParams.connParams['Ctx->PVN'] = {'synMech': 'exc', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'cellType': 'PVN_pop'}}

#I_ratio = 0.9 #PVP/PVN
#E_ratio = 0.6 #PVP/PVN 

target_secs = [sec for sec in netParams.cellParams['PVP_cell'].secs.keys() if sec is not 'soma']

#netParams.connParams['Ctx->STN'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 0.75, 'weight': 0.0002, 'delay': 1}
#netParams.connParams['Ctx->PVP'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': 'PVP_pop'}, 'probability':  0.7, 'weight': E_ratio * 0.0003, 'delay': 1}
#netParams.connParams['Ctx->PVN'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability':  0.5, 'weight': 0.0003, 'delay': 1}

#
# 
#netParams.connParams['GPe->PVP'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVP_pop'}, 'secs': target_secs,'probability':  0.5, 'weight': I_ratio * 0.00015, 'delay': 1}
#netParams.connParams['GPe->PVN'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability':  0.5, 'weight': 0.0001, 'delay': 1}
#netParams.connParams['GPeH->PVP'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_H_pop'}, 'postConds': {'pop': 'PVP_pop'},'secs': target_secs, 'probability':  0.3, 'weight': I_ratio * 0.00015, 'delay': 1}
#netParams.connParams['GPeH->PVN'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_H_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability':  0.5, 'weight': 0.000, 'delay': 1}
#netParams.connParams['GPeH->PVP'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_H_pop'}, 'postConds': {'pop': 'PVP_pop'},'probability':  0.5, 'weight': I_ratio * 0.02, 'delay': 1}
#weight -0.001 or lower,probability 0.5: burst-pause caused by disinhibition
#weight -0.0001 no burst, beta
#weight -0.0004 little rebound
#weight -0.0007 rebound




E_I_balance_PVP = 1 #1.3
E_I_balance_PVN = 1 #0.8

# get coord along STN axis (45 deg.)
deg45 = np.cos(np.pi/4)
proj_matrix = np.array([[deg45,  deg45], # inverse of 45-deg rotation matrix
                        [-deg45,  deg45]])
# project on axis, keep first coord only, and normalize
norm_ax = f'({deg45} * post_x + {deg45} * post_y) / sqrt(2)'

scale = 1
base_prob_E = 0.2/scale # before applying topographic rules, presume that each STN cell receives input with <base_prob_E/I> probability
base_prob_I = 0.2/scale

# Jeon et al, Fig 7H.2 (linear fit)
# Gria (exc): y = 0.55-0.3x
# Gabra (inh): y = 0.51-0.2x

# exaggerated (more dramatically decreasing alon the STN axis excitation, and spatially constant inhibition):
rule_glut =  f'{base_prob_E} * (0.65 - 0.5 * {norm_ax})'
rule_gaba = base_prob_I * 0.55



# target_secs = [sec for sec in netParams.cellParams['PVP_cell'].secs.keys() if sec is not 'soma']
# target_secs = [f'dend0_{i}' for i in range(10, 22)] + [f'dend1_{i}' for i in range(5, 10)]
netParams.connParams['Ctx->STN'] = {'synMech': 'glut', 'preConds': {'pop': 'Ctx_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': rule_glut, 'weight': 0.0003, 'delay': 1, 'sec': 'soma'}

#netParams.connParams['GPe->PVP'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVP_pop'}, 'secs': target_secs,'probability': rule_gaba, 'weight': E_I_balance_PVP * 0.0002, 'delay': 1}
#netParams.connParams['GPe->PVN'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability': rule_gaba, 'weight': E_I_balance_PVN * 0.0002, 'delay': 1}
netParams.connParams['GPeH->PVP'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_H_pop'}, 'postConds': {'pop': 'PVP_pop'},'secs': target_secs, 'probability': rule_gaba, 'weight': E_I_balance_PVP * 0.00015, 'delay': 1}
netParams.connParams['GPeH->PVN'] = {'synMech': 'GABA', 'preConds': {'pop': 'GPe_H_pop'}, 'postConds': {'pop': 'PVN_pop'}, 'probability': rule_gaba, 'weight': E_I_balance_PVN * 0.00015,  'delay': 1}




# # STN->STN - free parameter
netParams.probLengthConst = 300/scale

netParams.connParams['PVP->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVP_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5} # 'connList': [[0,1], [0,2]]

netParams.connParams['PVN->all'] = {'synMech': 'glut', 'preConds': {'pop': 'PVN_pop'}, 'postConds': {'pop': ['PVP_pop', 'PVN_pop']}, 'probability': 'exp(-dist_3D/probLengthConst) * 0.25', 'weight': 0.001, 'delay': 0.5} # exp(-dist_3D/probLengthConst) * 0.25'

# ------------------------------------------------ CONNECTIONS ------------------------------------------------
# to be defined..

# ------------------------------------------------ STIMULATION ------------------------------------------------

#cfg.IClamp_amp = -0.16
#netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': cfg.IClamp_amp}
#netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'PVP_cell'}, 'sec': 'soma', 'loc': 0.5}