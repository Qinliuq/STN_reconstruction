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
