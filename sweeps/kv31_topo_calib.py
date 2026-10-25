# Topographic healthy / parkinsonian network with the raised Kv3.1 defaults (PV+ 1.2, PV- 1.15): does the
# parkinsonian state survive, and which cortical conductance brings rates back to the earlier ~30 / 20 Hz?
duration = 3000
seeds = [1, 2, 3]
x_label = 'condition'
reference_group = 'healthy'
conditions = []
for ctx in [0.7, 0.8]:
    net = {'ctx_overrides': {'conductance_nS': ctx}}
    conditions += [
        {'label': f'ctx{ctx:g}_healthy', 'group': f'ctx{ctx:g}', 'x': 0, 'cfg': {**net, 'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
        {'label': f'ctx{ctx:g}_pd', 'group': f'ctx{ctx:g}', 'x': 1, 'cfg': {**net, 'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}},
    ]

# pinned to the Kir scales these runs used (defaults later moved to 0.3 in both types)
for c in conditions:
    c['cfg'] = {'PVP_gKir_scale': 0.5, 'PVN_gKir_scale': 0.1, **c['cfg']}

# pinned to the parkinsonian GPe->STN synapse scale these runs used (default later moved to 2.0)
for c in conditions:
    if c['cfg'].get('gpe_state') == 'parkinsonian':
        c['cfg']['gpe_params'] = {'synapse_scale': 1.5, **c['cfg'].get('gpe_params', {})}
