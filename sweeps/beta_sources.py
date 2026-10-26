# Which input carries the parkinsonian beta (input_version 2)? Healthy vs parkinsonian, and each source alone:
# pd_gpe_only = parkinsonian GPe + healthy cortex, pd_ctx_only = healthy GPe + parkinsonian cortex.
duration = 3000
seeds = [1, 2, 3]
x_label = 'condition'
conditions = [
    {'label': 'healthy', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
    {'label': 'pd_ctx_only', 'group': 'parkinsonian', 'x': 1, 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'healthy'}},
    {'label': 'pd_gpe_only', 'group': 'parkinsonian', 'x': 2, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'parkinsonian'}},
    {'label': 'pd', 'group': 'parkinsonian', 'x': 3, 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}},
]

# pinned to the model version these runs used (defaults later moved to input_version 3, pop_order 60, Kv3.1 1.2 / 1.15)
for c in conditions:
    c['cfg'] = {'input_version': 2, 'pop_order': 30, 'gpe_conductance_nS': 0.5, 'ctx_overrides': {},
                'PVP_gKv31_scale': 1.0, 'PVN_gKv31_scale': 1.0, **c['cfg']}

# pinned to the Kir scales these runs used (defaults later moved to 0.3 in both types)
for c in conditions:
    c['cfg'] = {'PVP_gKir_scale': 0.5, 'PVN_gKir_scale': 0.1, **c['cfg']}

# pinned to equal GPe->STN conductance in PV+ and PV- (default later moved to cfg.gpe_gaba_scale 1.13 / 0.87)
for c in conditions:
    c['cfg'] = {'gpe_gaba_scale': {'PVP': 1.0, 'PVN': 1.0}, **c['cfg']}
