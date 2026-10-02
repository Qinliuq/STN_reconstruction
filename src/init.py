from netpyne import sim
import sys

# TODO: can't avoid this hack?
isBatchRun = sys.argv[0][-5:] == 'nrniv'

if not isBatchRun:
    from tools import set_conductances, set_ion_styles, apply_CSF_Beurrier, apply_CSF_Bevan # for regular run
else:
    from src.tools import set_conductances, set_ion_styles, apply_CSF_Beurrier, apply_CSF_Bevan # for batch run

cfg, netParams = sim.readCmdLineArgs(simConfigDefault='src/cfg.py', netParamsDefault='src/netParams.py')

sim.create(netParams, cfg)

for cell in sim.net.cells:
    if cell.tags.get('cellType') in ['PVP_cell', 'PVN_cell']:
        set_conductances(cell, cfg)
set_ion_styles()

# apply_CSF_Beurrier()
apply_CSF_Bevan()

sim.simulate()
sim.analyze()

