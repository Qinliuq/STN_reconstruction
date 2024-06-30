from netpyne import specs

from __main__ import cfg

netParams = specs.NetParams()

from neuron import h
resultCode = h.load_file('cells/SThprotocell.hoc')
if resultCode:
    protoCell = h.SThproto()
else:
    raise RuntimeError("Couldn't load prototype STN cell from `cells/STHprotocell.hoc")

stnParams = netParams.importCellParams(label='STN_cell', fileName='cells/SThprotocell.hoc', cellName='SThcell', cellArgs=[0, protoCell]) # 0 in cellArgs is just placeholder
# # modifing arbitrary properties, e,g:
# for iSec, sec in stnParams.secs.items():
#     sec.mechs.Kv31.gk = 0

stnParams.secs.soma['threshold'] = -30

netParams.popParams['STN_pop'] = {'cellType': 'STN_cell', 'numCells': cfg.num_vals}

# # rebound burst
# cfg.IClamp_amp = -0.35
# netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 600, 'dur': 1000, 'amp': cfg.IClamp_amp}

# slow bursting
cfg.IClamp_amp = -0.25
netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': cfg.IClamp_amp}

# # # fast bursting
# cfg.IClamp_amp = -0.35
# netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': cfg.IClamp_amp}

netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'STN_cell'}, 'sec': 'soma', 'loc': 0.5}

# # alternatively, to search through multiple amps:
# import numpy as np
# iapps = np.linspace(0, 10, cfg.num_vals) / 100
# for i, amp in enumerate(iapps):
#     netParams.stimSourceParams[f'clamp{i}'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': amp}
#     netParams.stimTargetParams[f'clamp{i}->{i}'] = {'source': f'clamp{i}', 'conds': {'cellList': [i]}, 'sec': 'soma', 'loc': 0.5}