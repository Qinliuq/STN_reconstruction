# Calibration of the realistic parkinsonian state (input_version 3, pop_order 60): topographic GPe-STN width
# (topo_sigma) x base GPe->STN conductance x parkinsonian GPe->STN synapse scale; healthy reference per
# (sigma, conductance). Targets: parkinsonian more bursty than healthy (>= 1.5x), beta and synchrony peaking
# at the upper-middle STN axis (~0.65), plausible healthy rates.
duration = 3000
seeds = [1]
x_label = 'parkinsonian GPe->STN synapse scale (healthy at 1)'
conditions = []
for sigma in [0.10, 0.15]:
    for g in [1.0, 2.0]:
        tag = f's{sigma:g}_g{g:g}'
        net = {'gpe_conductance_nS': g}
        conditions.append({'label': f'{tag}_healthy', 'group': tag, 'x': 1.0,
                           'cfg': {**net, 'ctx_state': 'healthy', 'gpe_state': 'healthy', 'gpe_params': {'topo_sigma': sigma}}})
        for scale in [1.5, 2.5]:
            conditions.append({'label': f'{tag}_pd{scale:g}', 'group': tag, 'x': scale,
                               'cfg': {**net, 'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian',
                                       'gpe_params': {'topo_sigma': sigma, 'synapse_scale': scale}}})

# pinned to the cortical conductance these runs used (default later moved to 0.7 nS)
for c in conditions:
    c['cfg'] = {'ctx_overrides': {}, **c['cfg']}
