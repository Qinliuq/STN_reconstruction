# HCN and CaT block, alone and combined with Kir block (same factor for both channels), in the
# parkinsonian network with input_version 2 (GPe is the primary beta source). The healthy, untreated
# and Kir-only conditions are the same runs as in sweeps/kir_block_gpe_beta.py.
duration = 3000
seeds = [1, 2, 3]
x_label = 'conductance factor (1 = default, 0 = full block)'
reference_group = 'healthy'
baseline_label = 'untreated'
band_groups = ['healthy', 'untreated']
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
              {'label': 'untreated', 'group': 'untreated', 'x': 1.0, 'cfg': {**pd}}]
blocks = {'Kir': ['gKir'], 'HCN': ['gHCN'], 'CaT': ['gCaT'], 'HCN+Kir': ['gHCN', 'gKir'], 'CaT+Kir': ['gCaT', 'gKir']}
for name, channels in blocks.items():
    for f in [0.75, 0.5, 0.25, 0.0]:
        conditions.append({'label': f"{name.lower().replace('+', '_')}{f:g}", 'group': name, 'x': f,
                           'cfg': {**pd, 'channel_block': {ch: f for ch in channels}}})

# pinned to the model version these runs used (defaults later moved to input_version 3, pop_order 60, Kv3.1 1.2 / 1.15)
for c in conditions:
    c['cfg'] = {'input_version': 2, 'pop_order': 30, 'gpe_conductance_nS': 0.5, 'ctx_overrides': {},
                'PVP_gKv31_scale': 1.0, 'PVN_gKv31_scale': 1.0, **c['cfg']}

# pinned to the Kir scales these runs used (defaults later moved to 0.3 in both types)
for c in conditions:
    c['cfg'] = {'PVP_gKir_scale': 0.5, 'PVN_gKir_scale': 0.1, **c['cfg']}
