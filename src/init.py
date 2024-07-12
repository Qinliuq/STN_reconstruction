from netpyne import sim
import sys

# TODO: can't avoid this hack?
isBatchRun = sys.argv[0][-5:] == 'nrniv'

if not isBatchRun:
    from tools import set_values_from_file, apply_CSF_Beurrier, apply_CSF_Bevan # for regular run
else:
    from src.tools import set_values_from_file, apply_CSF_Beurrier, apply_CSF_Bevan # for batch run

cfg, netParams = sim.readCmdLineArgs(simConfigDefault='src/cfg.py', netParamsDefault='src/netParams.py')

sim.create(netParams, cfg)

import numpy as np

# using PV+ "best candidate" params:
cat, cal, hcn, sk = 1.4, 0.86, 0.65, 0.9 # 28/20, */5, 1.3/2, 0.9

for i, cell in enumerate(sim.net.cells):

    cellType = cell.tags.get('cellType')
    if cellType not in ['PVP_cell', 'PVN_cell']:
        continue

    if cellType == 'PVP_cell': # PV+ params
        gCaT_scale = cat
        gCaL_scale = cal
        gHCN_scale = hcn
        gSK_scale = sk
    else: # PV- params
        gCaT_scale = .8
        gCaL_scale = 1.2
        gHCN_scale = .6
        gSK_scale = 1

    # Na   
    default_gNa_soma = 1.483419823e-02 
    default_gNa_dend = 1.0e-7
    cell.secs.soma.hObj.gna_Na = default_gNa_soma
    # NaL (can int#rchange with the Do & Bean model)
    default_gNaL_soma = 1.108670852e-05
    default_gNaL_dend = 0.81e-5
    cell.secs.soma.hObj.gna_NaL = default_gNaL_soma

    # linear conductances (loaded from files)...
    set_values_from_file(cell, "gk_KDR")  
    set_values_from_file(cell, "gk_Kv31")
    set_values_from_file(cell, "gcaN_HVA")
    set_values_from_file(cell, "gk_Ih", soma_scale=gHCN_scale, dend_scale=gHCN_scale)
    set_values_from_file(cell, "gcaT_CaT", soma_scale=gCaT_scale, dend_scale=gCaT_scale)
    # set_values_from_file(cell, "gcaT_CaT", dend_scale=1.2) # modulate CaT as e.g. in Fig. 6C
    set_values_from_file(cell, "gk_sKCa", soma_scale=gSK_scale, dend_scale=gSK_scale)
    # set_values_from_file(cell, "gk_sKCa", "apamin0.9") # apply apamin
    set_values_from_file(cell, "gcaL_HVA", soma_scale=gCaL_scale, dend_scale=gCaL_scale)
    # set_values_from_file(cell, "gcaL_HVA", "dl0.9") # 10% decrease in dendritic linear CaL (see Figure 6A, 8A)

from neuron import h
for sec in h.allsec():
    sec.push()
    h.ion_style("na_ion",1,2,1,0,1)
    h.ion_style("k_ion",1,2,1,0,1)
    h.ion_style("ca_ion",3,2,1,1,1)   
    h.pop_section()

# apply_CSF_Beurrier()
apply_CSF_Bevan()

sim.simulate()
sim.analyze()

