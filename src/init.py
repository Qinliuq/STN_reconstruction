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

saveFolder = 'bursts_37_bevan'

# using PV+ "best candidate" params:
cat, cal, hcn, sk = 1.4, 0.86, 0.65, 0.9 # 28/20, */5, 1.3/2, 0.9
cats, cals, hcns, sks = [cat] * cfg.num_vals, [cal] * cfg.num_vals, [hcn] * cfg.num_vals, [sk] * cfg.num_vals # for consistency with more general way below
cfg.filename = f'{saveFolder}/cat{cat}_cal{cal}_hcn{hcn}_sk{sk}_iclamp_{cfg.IClamp_amp}'

# # or, using search through either of parameters' values (others will be fixed):
# # E.g., search through CaL, as in Fig. 1 in park et al.
# cals = np.linspace(0, 40, cfg.num_vals)
# cals /= 5 # normalize with respect to default CaL value (i.e. 5) from Park et al.
# cfg.filename = f'{saveFolder}/cat{cat}_cals{cals[0]}-{cals[-1]}_hcn{hcn}_sk{sk}_iclamp_{cfg.IClamp_amp}'

for i, cell in enumerate(sim.net.cells):
    # cell = sim.net.cells[0]
    gCaT_scale = cats[i] if cats is not None else 1
    gCaL_scale = cals[i] if cals is not None else 1
    gHCN_scale = hcns[i] if hcns is not None else 1
    gSK_scale = sks[i]

    # Na   
    default_gNa_soma = 1.483419823e-02 
    default_gNa_dend = 1.0e-7
    cell.secs.soma.hObj.gna_Na = default_gNa_soma
    # NaL (can int#rchange with the Do & Bean model)
    default_gNaL_soma = 1.108670852e-05
    default_gNaL_dend = 0.81e-5
    cell.secs.soma.hObj.gna_NaL = default_gNaL_soma

    if isBatchRun:
        # using values from cfg that are variable across batch iterations
        gCaT_scale = cfg.gCaT_scale
        gSK_scale = cfg.gSK_scale

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

