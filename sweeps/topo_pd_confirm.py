# Confirmation of the best calibration setting (topo_sigma 0.1, base GPe->STN 1 nS, parkinsonian synapse
# scale 1.5 = gpe_feedback_v3 default) with more seeds, and a lower cortical conductance to bring rates down.
duration = 3000
seeds = [1, 2, 3]
x_label = 'condition'
reference_group = 'healthy'
base = {'gpe_conductance_nS': 1.0, 'gpe_params': {'topo_sigma': 0.10}}
conditions = []
for ctx in [1.0, 0.7]:
    net = {**base, 'ctx_overrides': {'conductance_nS': ctx}}
    conditions += [
        {'label': f'ctx{ctx:g}_healthy', 'group': f'ctx{ctx:g}', 'x': 0, 'cfg': {**net, 'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
        {'label': f'ctx{ctx:g}_pd', 'group': f'ctx{ctx:g}', 'x': 1, 'cfg': {**net, 'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}},
    ]

# pinned to the Kv3.1 scales these runs used (defaults later moved to 1.2 / 1.15)
for c in conditions:
    c['cfg'] = {'PVP_gKv31_scale': 1.0, 'PVN_gKv31_scale': 1.0, **c['cfg']}

# pinned to the Kir scales these runs used (defaults later moved to 0.3 in both types)
for c in conditions:
    c['cfg'] = {'PVP_gKir_scale': 0.5, 'PVN_gKir_scale': 0.1, **c['cfg']}

# pinned to the parkinsonian GPe->STN synapse scale these runs used (default later moved to 2.0)
for c in conditions:
    if c['cfg'].get('gpe_state') == 'parkinsonian':
        c['cfg']['gpe_params'] = {'synapse_scale': 1.5, **c['cfg'].get('gpe_params', {})}

# pinned to equal GPe->STN conductance in PV+ and PV- (default later moved to cfg.gpe_gaba_scale 1.13 / 0.87)
for c in conditions:
    c['cfg'] = {'gpe_gaba_scale': {'PVP': 1.0, 'PVN': 1.0}, **c['cfg']}
