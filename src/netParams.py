

from netpyne import specs

netParams = specs.NetParams()

from neuron import h
resultCode = h.load_file('cells/SThprotocell.hoc')
if resultCode:
    protoCell = h.SThproto()
else:
    raise RuntimeError("Couldn't load prototype STN cell from `cells/STHprotocell.hoc")

stnParams = netParams.importCellParams(label='STN_cell', fileName='cells/SThprotocell.hoc', cellName='SThcell', cellArgs=[0, protoCell]) # 0 in cellArgs is just placeholder
stnParams.secs.soma['threshold'] = -30

netParams.popParams['STN_pop'] = {'cellType': 'STN_cell', 'numCells': 1}

netParams.synMechParams['AMPA'] = {
    'mod': 'Exp2Syn',
    'tau1': 0.1,
    'tau2': 1.0,
    'e': 0
}

netParams.stimSourceParams['stim'] = {'type': 'NetStim', 'interval': 10, 'number': 4, 'start': 10}
netParams.stimTargetParams['stim->STN'] = {
    'source': 'stim', 
    'synMech': 'AMPA',
    'weight': 0.05,
    'sec': 'soma', 
    'loc': 0.5, 
    'conds': {'pop': 'STN_pop'}
}