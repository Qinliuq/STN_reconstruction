from netpyne import sim

cfg, netParams = sim.loadFromIndexFile('index.npjson')

sim.create(netParams, cfg)

from tools import set_values_from_file, apply_CSF_Beurrier, apply_CSF_Bevan
cell = sim.net.cells[0]

# Na   
default_gNa_soma = 1.483419823e-02 
default_gNa_dend = 1.0e-7
cell.secs.soma.hObj.gna_Na = default_gNa_soma
# NaL (can interchange with the Do & Bean model)
default_gNaL_soma = 1.108670852e-05
default_gNaL_dend = 0.81e-5
cell.secs.soma.hObj.gna_NaL = 1.108670852e-05
   
# linear conductances (loaded from files)...
set_values_from_file(cell, "gk_KDR")  
set_values_from_file(cell, "gk_Kv31")
set_values_from_file(cell, "gk_Ih")

set_values_from_file(cell, "gcaT_CaT")
# set_values_from_file(cell, "gcaT_CaT", dend_scale=1.2) # (modulate CaT as in Fig. 6C)

set_values_from_file(cell, "gcaN_HVA")
# set_values_from_file(cell, "gcaL_HVA", "d0.9")

# set_values_from_file(cell, "gk_sKCa") # no apamin
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

