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

cell = sim.net.cells[0]

# Na   
default_gNa_soma = 1.483419823e-02 
default_gNa_dend = 1.0e-7
cell.secs.soma.hObj.gna_Na = default_gNa_soma
# NaL (can interchange with the Do & Bean model)
default_gNaL_soma = 1.108670852e-05
default_gNaL_dend = 0.81e-5
cell.secs.soma.hObj.gna_NaL = default_gNaL_soma

if isBatchRun:
    # using values from cfg that are variable across batch iterations
    gCaT_scale = cfg.gCaT_scale
    gSK_scale = cfg.gSK_scale
else:
    # using defaults
    gCaT_scale = 1
    gSK_scale = 1

# linear conductances (loaded from files)...
set_values_from_file(cell, "gk_KDR")  
set_values_from_file(cell, "gk_Kv31")
set_values_from_file(cell, "gk_Ih")
set_values_from_file(cell, "gcaN_HVA")

set_values_from_file(cell, "gcaT_CaT", soma_scale=gCaT_scale, dend_scale=gCaT_scale)
# set_values_from_file(cell, "gcaT_CaT", dend_scale=1.2) # modulate CaT as e.g. in Fig. 6C

set_values_from_file(cell, "gcaL_HVA")
# set_values_from_file(cell, "gcaL_HVA", "dl0.9") # 10% decrease in dendritic linear CaL (see Figure 6A, 8A)

set_values_from_file(cell, "gk_sKCa", soma_scale=gSK_scale, dend_scale=gSK_scale)
set_values_from_file(cell, "gk_sKCa", "apamin0.9") # apply apamin

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

