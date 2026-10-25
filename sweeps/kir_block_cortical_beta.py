# Kir block in the parkinsonian network with input_version 1 (beta from cortex), as run on 2026-10-17
# (results/kir_sweep_cortical_beta/). Healthy reference: same network, healthy cortical input.
duration = 3000
seeds = [1, 2, 3]
x_label = 'Kir factor (1 = default, 0 = full block)'
reference_group = 'healthy'
baseline_label = 'pd_kir1'
base = {'input_version': 1, 'ctx_overrides': {'conductance_nS': 1.0}, 'gpe_syn': {}}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {**base, 'ctx_state': 'healthy'}}]
conditions += [{'label': f'pd_kir{f:g}', 'group': 'parkinsonian', 'x': f,
                'cfg': {**base, 'ctx_state': 'parkinsonian', 'channel_block': {'gKir': f}}} for f in [1.0, 0.75, 0.5, 0.25, 0.0]]

# pinned to the model version these runs used (defaults later moved to input_version 3, pop_order 60, Kv3.1 1.2 / 1.15)
for c in conditions:
    c['cfg'] = {'input_version': 2, 'pop_order': 30, 'gpe_conductance_nS': 0.5, 'ctx_overrides': {},
                'PVP_gKv31_scale': 1.0, 'PVN_gKv31_scale': 1.0, **c['cfg']}

# pinned to the Kir scales these runs used (defaults later moved to 0.3 in both types)
for c in conditions:
    c['cfg'] = {'PVP_gKir_scale': 0.5, 'PVN_gKir_scale': 0.1, **c['cfg']}
