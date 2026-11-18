from netpyne import sim
import sys

# TODO: can't avoid this hack?
isBatchRun = sys.argv[0][-5:] == 'nrniv'

if not isBatchRun:
    from tools import set_conductances, set_ion_styles, apply_CSF_Beurrier, apply_CSF_Bevan # for regular run
    from stn_inputs import GPeFeedback
    from network_metrics import save_metrics
else:
    from src.tools import set_conductances, set_ion_styles, apply_CSF_Beurrier, apply_CSF_Bevan # for batch run
    from src.stn_inputs import GPeFeedback
    from src.network_metrics import save_metrics

cfg, netParams = sim.readCmdLineArgs(simConfigDefault='src/cfg.py', netParamsDefault='src/netParams.py')

sim.create(netParams, cfg)

for cell in sim.net.cells:
    if cell.tags.get('cellType') in ['PVP_cell', 'PVN_cell']:
        set_conductances(cell, cfg)
set_ion_styles()

# apply_CSF_Beurrier()
apply_CSF_Bevan()

# GPe input: closed loop with the STN spikes, so the simulation is run in steps of cfg.gpe_dt (replaces sim.simulate())
gpe = GPeFeedback(sim, cfg, netParams.synMechParams['GABA'])
gpe.run()
sim.gatherData()
gpe.save()
save_metrics(sim, cfg, gpe)
gpe.plot(f'{sim.cfg.filename}_GPe_rates.png')
sim.analyze()

