

from netpyne import specs

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

netParams.popParams['STN_pop'] = {'cellType': 'STN_cell', 'numCells': 1}

# # rebound burst
# netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 770, 'dur': 500, 'amp': -0.25}
# # slow bursting
# netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': -0.25}
# fast bursting
netParams.stimSourceParams['clamp'] = {'type': 'IClamp', 'del': 0, 'dur': 40000, 'amp': -0.35}

netParams.stimTargetParams['clamp->all'] = {'source': 'clamp', 'conds': {'cellType': 'STN_cell'}, 'sec': 'soma', 'loc': 0.5}