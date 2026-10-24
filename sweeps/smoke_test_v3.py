# quick check of input_version 3 at the default network size (pop_order 60)
duration = 600
seeds = [1]
conditions = [
    {'label': 'healthy', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
    {'label': 'pd', 'group': 'pd', 'x': 1, 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}},
]

# pinned to the Kv3.1 scales and cortical conductance these runs used (defaults later moved to 1.2 / 1.15, 0.7 nS)
for c in conditions:
    c['cfg'] = {'PVP_gKv31_scale': 1.0, 'PVN_gKv31_scale': 1.0, 'ctx_overrides': {}, 'gpe_conductance_nS': 0.5, **c['cfg']}
