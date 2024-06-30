from neuron import h

def set_values_from_file(cell, conductance, modifier=None, soma_scale=1.0, dend_scale=1.0):
    # use `modifier` to replicate some of paper result.
    # Using `soma_scale` and/or `dend_scale` should theoretically have the same effect, although in practice it's not always the case

    modif = f'-{modifier}' if modifier else ''
    f = open(f'cells/sth-data/cell-{conductance}{modif}', 'r')
    totval=0
    numval = 0
    for line in f:
        tree, sec, loc, val = line.strip().split(' ')
        loc = float(loc)
        val = float(val)

        if tree == '-1': # 'soma'
            hObj = cell.secs.soma.hObj
            val *= soma_scale
        else:
            hObj = getattr(cell.secs, f'dend{tree}_{sec}').hObj
            val *= dend_scale
            totval += val
            numval += 1
        setattr(hObj(loc), conductance, val)
    # print(f"Cond {conductance} avg {totval/numval} (n. {numval})")

def apply_CSF_Beurrier():
    print("Setting in vitro parameters based on Beurrier et al (1999)")

    h.nai0_na_ion = 15
    h.nao0_na_ion = 150

    h.ki0_k_ion = 140
    h.ko0_k_ion = 3.6

    h.cai0_ca_ion = 1e-04
    h.cao0_ca_ion = 2.4

    # # This was in original model, but seems to be not relevant
    # h.cli0_cl_ion = 4
    # h.clo0_cl_ion = 135

def apply_CSF_Bevan():
    print("Setting in vitro parameters based on Bevan & Wilson (1999)")

    h.nai0_na_ion = 15
    h.nao0_na_ion = 128.5

    h.ki0_k_ion = 140
    h.ko0_k_ion = 2.5

    h.cai0_ca_ion = 1e-04
    h.cao0_ca_ion = 2.0

    # # This was in original model, but seems to be not relevant
    # cli0_cl_ion = 4
    # clo0_cl_ion = 132.5