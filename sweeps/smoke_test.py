# quick check that both input versions run through the sweep tool
duration = 600
seeds = [1]
conditions = [
    {'label': 'v1_pd', 'group': 'v1', 'x': 1, 'cfg': {'input_version': 1, 'ctx_overrides': {'conductance_nS': 1.0}, 'ctx_state': 'parkinsonian'}},
    {'label': 'v2_pd_kir0.5', 'group': 'v2', 'x': 2, 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian', 'channel_block': {'gKir': 0.5}}},
]

# pinned to the model version these runs used (defaults later moved to input_version 3, pop_order 60)
for c in conditions:
    c['cfg'] = {'input_version': 2, 'pop_order': 30, 'gpe_conductance_nS': 0.5, 'ctx_overrides': {}, **c['cfg']}
